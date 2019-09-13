% Based on earlier DFN CasADi implementation by SHP, Bosch Year 1
% By: Zach Gima 2019-9-12

% INPUTS --- 
%%% FOR JUST VOLTAGE SIMULATION
% (1) p: structure containing all of the parameters needed for the model
% (2) data: struct with fields cur (current), time, V0

%%% FOR VOLTAGE + SENSITIVITY SIMULATION
% (3) varargin{1} -- theta_0: initial numerical values for the parameters we'd like to
% identify
% (4) theta_str: string format of the names of the parameters we'd like to
% identify; match the field name associated w/ that param in the p struct

% OUTPUTS --- 
% (1) v_sim: voltage 
% (2) state_info: additional internal state info of the model
% (3) varargout{1} -- sens: sensitivity of parameters (p) passed into the
% function (partial V / partial theta) 

function [v_sim,state_info,varargout] = dfn_casadi(p, data, varargin)
    import casadi.*

    %% Parse Inputs
    I = data.cur/p.Area;
    t = data.time;
    V0 = data.V0;
    T_amb = data.T_amb;

    dt =  t(2)-t(1);
    NT = length(t);
    
    % varargin used if computing sensitivity / identifying parameters in the model
    if nargin > 2
        theta_0 = varargin{1};
        theta_str = varargin{2};
        SensFlag = 1;
        np = length(theta_0);

        % Generate & Assign CasADi Variables to p struct
        [p,theta_sx] = update_casadi_vars(p,theta_str);
        
        % Update parameter dependencies
        p = update_dependencies(p);
                
        % Electrolyte concentration, solid & elyte electric potential, and
        % elyte current matrices (function of other params)
        p = update_mats(p,SensFlag);
        
    else % just simulating the model
        theta_sx = [];
        theta_0 = [];
        SensFlag = 0;
        
        % Update parameter dependencies
        p = update_dependencies(p);
        
        % Electrolyte concentration, solid & elyte electric potential, and
        % elyte current matrices (function of other params)
        p = update_mats(p,SensFlag);
    end

    %% (DFN Code Copy) Initial Conditions, Preallocation
    [x0_nom,z0_nom,SOC,c_avg_n] = init_dfn(p,NT,V0,T_amb,SensFlag);  

    %% Initialize dae system in casadi
    % Define States
    x = SX.sym('x',size(x0_nom));
    
    % Define algebraic variables
    z = SX.sym('z',size(z0_nom));

    % Input
    u = SX.sym('u'); 

    %% DAE Builder
%     [x_dot, g_, L ] = dae_dfn(x,z,u,p); % [SHP Change]
    [x_dot, g_, L, x_outs, z_outs, info_outs, param_outs ] = dae_dfn(x,z,u,p);

    %% Sensitivity Calculation (if turned on)

    % Initial condition
    x0_call = Function('x0_call',{theta_sx},{x0_nom},{'theta_sx'},{'result'});
    x0 = full(x0_call(theta_0));

    if SensFlag == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1. Calculate Jacobian automatically
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Jac_V0_param = jacobian(L,theta_sx);
        V0_call = Function('V0_call',{theta_sx,x,u},{Jac_V0_param},{'theta_sx','x','u'},{'result'});
        S3_0_V0_params = full(V0_call(theta_0,x0,0));

        for j=1:np
            if strcmp(theta_sx(j),'c_e') % c_e requires special handling
                Jac_x0_ce0 = jacobian(x0,p.c_e);
                Jac_x0_ce0_call = Function('Jac_x0_ce0_call',{theta_sx},{Jac_x0_ce0},{'theta_sx'},{'result'});
                S1_0_x0_ce0 = full(Jac_x0_ce0_call(theta_0));
            end
        end

        % Calculate Jacobian matrix (\frac{\partial f}{\partial x}}), A11
        f_x_jac = jacobian(x_dot,x);

        % Calculate Jacobian matrix (\frac{\partial f}{\partial z}}), A12
        f_z_jac = jacobian(x_dot,z);

        % Calculate Jacobian matrix (\frac{\partial f}{\partial \theta}}), B1
        f_theta_jac = jacobian(x_dot,theta_sx);

        % Calculate Jacobian matrix (\frac{\partial g}{\partial x}}), A21
        g_x_jac = jacobian(g_,x);

        % Calculate Jacobian matrix (\frac{\partial g}{\partial z}}), A22
        g_z_jac = jacobian(g_,z);

        % Calculate Jacobian matrix (\frac{\partial g}{\partial \theta}}), B2
        g_theta_jac = jacobian(g_,theta_sx);

        % Calculate Jacobian matrix (\frac{\partial h}{\partial x}}), C
        h_x_jac = jacobian(L,x);

        % Calculate Jacobian matrix (\frac{\partial h}{\partial z}}), D
        h_z_jac = jacobian(L,z);

        % Calculate Jacobian matrix (\frac{\partial h}{\partial \theta}}), E
        h_theta_jac = jacobian(L,theta_sx);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. Build sensitivity equation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        S1_0 = zeros(size(x,1),np);
        for j=1:np
            if strcmp(theta_sx(j),'c_e')% c_e index
                S1_0(:,j) = S1_0_x0_ce0;
            end
        end
        S2_0 = zeros(size(z,1),np);
        S3_0 = S3_0_V0_params;

        % Symbolic Sensitivity

        %%% S1_dot
        S1_blk = SX.sym('S1_blk',size(x,1)*np,1);
        S1 = reshape(S1_blk,size(x,1),np);

        S2_blk = SX.sym('S2_blk',size(z,1)*np,1);
        S2 = reshape(S2_blk,size(z,1),np);

        tmp = repmat({f_x_jac},np,1);
        f_x_blk = blkdiag(tmp{:});

        tmp = repmat({f_z_jac},np,1);
        f_z_blk = blkdiag(tmp{:});

        f_theta_blk = reshape(f_theta_jac,size(x,1)*np,1);

        S1_dot = f_x_blk*S1_blk + f_z_blk*S2_blk + f_theta_blk;

        %%% S2_
        tmp = repmat({g_x_jac},np,1);
        g_x_blk = blkdiag(tmp{:});

        tmp = repmat({g_z_jac},np,1);
        g_z_blk = blkdiag(tmp{:});

        g_theta_blk = reshape(g_theta_jac,size(z,1)*np,1);

        S2_ = g_x_blk*S1_blk + g_z_blk*S2_blk + g_theta_blk;

        %%% S3
        S3 = h_x_jac * S1 + h_z_jac * S2 + h_theta_jac;
        S3 = S3'; % change row-vector to column vector. 

        % IC for S1_blk, S2_blk
        S1_0_blk(:,1) = reshape(S1_0,size(x,1)*np,1);
        S2_0_blk(:,1) = reshape(S2_0,size(z,1)*np,1);

        %%% Explicitly define sensitivity values needed for simulation
        S1_sim_blk(:,1) = S1_0_blk(:,1);
        S2_sim_blk(:,1) = S2_0_blk(:,1);
        
        sens = zeros(NT,np);
        S3_sim = zeros(np,NT); % necessary here regardless of sensflag
        S3_sim(:,1) = S3_0';
    end
    
    %% Integrator

    % Note: Sensitivity wouldnt converge/work for higher tolerance (10^-6) but
    % the DFN did, so separate CasADi integrator functions made for the DFN (F)
    % and Sensitivity DAEs (SS)
    opts1 = struct('tf',dt, 'abstol',1e-6,'reltol',1e-6);
    opts2 = struct('tf',dt, 'abstol', 1e-1, 'reltol', 1e-1);

    dae = struct('x',x, 'z',z, 'p',[theta_sx;u], 'ode',x_dot, 'alg', g_, 'quad', L);
    F = integrator('F','idas', dae, opts1);

    if SensFlag == 1
        dae_s = struct('x',S1_blk, 'z', S2_blk, 'p',[theta_sx;x;z;u], 'ode', S1_dot, 'alg', S2_, 'quad', S3);
        SS = integrator('SS', 'idas', dae_s, opts2);
    end

    %% Build function for simulated states % [SHP Change]
    f_out = Function('f_out',{theta_sx,x,z,u},{x_outs},{'theta_sx','x','z','u'},{'f_out'});
    g_out = Function('g_out',{theta_sx,x,z,u},{z_outs},{'theta_sx','x','z','u'},{'g_out'});
    alg_out = Function('alg_out',{theta_sx,x,z,u},{info_outs},{'theta_sx','x','z','u'},{'alg_out'});
    par_out = Function('par_out',{theta_sx,x,z,u},{param_outs},{'theta_sx','x','z','u'},{'par_out'});
    f0 = full(f_out(theta_0,x0,z0_nom,0));
    g0 = full(g_out(theta_0,x0,z0_nom,0));
    a0 = full(alg_out(theta_0,x0,z0_nom,0));
    p0 = full(par_out(theta_0,x0,z0_nom,0));

    %% Indexing

    % index for x
    out_csn_idx = 1:(p.PadeOrder * (p.Nxn-1)); 
    out_csp_idx = (p.PadeOrder * (p.Nxn-1) + 1) : (p.PadeOrder * (p.Nxn-1) + p.PadeOrder * (p.Nxp-1));
    out_ce_idx = (p.PadeOrder * (p.Nxn-1) + p.PadeOrder * (p.Nxp-1) + 1) : (p.PadeOrder * (p.Nxn-1) + p.PadeOrder * (p.Nxp-1) + p.Nxn-1+p.Nxs-1+p.Nxp-1);
    out_T_idx = (p.PadeOrder * (p.Nxn-1) + p.PadeOrder * (p.Nxp-1) + p.Nxn-1+p.Nxs-1+p.Nxp-1 + 1);
    
    % index for z
    out_phisn_idx = 1:p.Nxn-1; % 1:Nn
    out_phisp_idx = (p.Nxn-1+1) : (p.Nxn-1 + p.Nxp-1); % Nn+1 : Nnp
    out_ien_idx = (p.Nxn-1 + p.Nxp-1 +1) : (p.Nxn-1 + p.Nxp-1 + p.Nxn-1);
    out_iep_idx = (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 +1) : (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1);
    out_phie_idx = (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1 + 1) : (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1 + (p.Nxn-1 + p.Nxs-1 + p.Nxp-1)+2); % Why 2? not 4?? used to plus 4 to consider boundary conditions.
    out_jn_idx = (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1 + (p.Nxn-1 + p.Nxs-1 + p.Nxp-1 + 2) + 1) : (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1 + (p.Nxn-1 + p.Nxs-1 + p.Nxp-1 + 2) + (p.Nxn-1));
    out_jp_idx = (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1 + (p.Nxn-1 + p.Nxs-1 + p.Nxp-1 + 2) + (p.Nxn-1) +1) : (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1 + (p.Nxn-1 + p.Nxs-1 + p.Nxp-1+ 2) + (p.Nxn-1) + (p.Nxp-1));

    % index for information
    out_cssn_idx = 1:p.Nxn-1 ;
    out_cssp_idx = (p.Nxn-1+1) : (p.Nxn-1 + p.Nxp-1);
    out_cex_idx = (p.Nxn-1 + p.Nxp-1 +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4));
    out_theta_avgn_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1));
    out_theta_avgp_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1));
    out_etan_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1));
    out_etap_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) + (p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1));
    out_ce0n_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +1);
    out_ce0p_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +2) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +2);
    out_etas_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +3) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +3);
    out_Volt_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +4) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +4);
    out_nLis_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +5) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +5);
    out_nLie_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +6) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +6);

    % index for state-dep param
    out_Den_idx = 1:(p.Nxn-1) ;
    out_dDen_idx = (p.Nxn-1 +1) : (p.Nxn-1 + p.Nxn-1);
    out_Des_idx = (p.Nxn-1 + p.Nxn-1 +1) : (p.Nxn-1 + p.Nxn-1 + p.Nxs-1);
    out_dDes_idx = (p.Nxn-1 + p.Nxn-1 + p.Nxs-1 +1) : (p.Nxn-1 + p.Nxn-1 + p.Nxs-1 + p.Nxs-1);
    out_Dep_idx = (p.Nxn-1 + p.Nxn-1 + p.Nxs-1 + p.Nxs-1 +1) : (p.Nxn-1 + p.Nxn-1 + p.Nxs-1 + p.Nxs-1 + p.Nxp -1);
    out_dDep_idx = (p.Nxn-1 + p.Nxn-1 + p.Nxs-1 + p.Nxs-1 + p.Nxp-1 +1) : (p.Nxn-1 + p.Nxn-1 + p.Nxs-1 + p.Nxs-1 + p.Nxp-1 + p.Nxp-1);

    %% Simulation
    % pre-allocate
    v_sim = zeros(NT,1);
    T1_sim = zeros(NT,1);
    T2_sim = zeros(NT,1);
    
    % Need to set these to value for error causing script to exit reasons
    if SensFlag == 0
        S3_sim = 0;
        sens = 0;
    end

    % Initialize
    x_sim(:,1) = x0;
    z_sim(:,1) = z0_nom;
    v_sim(1) = V0;

    % [SHP Change]
    % save x
    csn_sim(:,1) = f0(out_csn_idx);
    csp_sim(:,1) = f0(out_csp_idx);
    ce_sim(:,1) = f0(out_ce_idx);
    cen_sim(:,1) = ce_sim(1:(p.Nxn-1),1);
    ces_sim(:,1) = ce_sim((p.Nxn-1 +1):(p.Nxn-1 + p.Nxs-1),1);
    cep_sim(:,1) = ce_sim((p.Nxn-1 + p.Nxs-1 +1):(p.Nxn-1 + p.Nxs-1 + p.Nxp-1),1);
    T1_sim(1) = f0(out_T_idx);
    T2_sim(1) = f0(end);
    
    % save z
    phi_s_n_sim(:,1) = g0(out_phisn_idx);
    phi_s_p_sim(:,1) = g0(out_phisp_idx);
    ien_sim(:,1) = g0(out_ien_idx);
    iep_sim(:,1) = g0(out_iep_idx);
    phie_sim(:,1) = g0(out_phie_idx);
    jn_sim(:,1) = g0(out_jn_idx);
    jp_sim(:,1) = g0(out_jp_idx);

    % save alg. states
    cssn_sim(:,1) = a0(out_cssn_idx);
    cssp_sim(:,1) = a0(out_cssp_idx);
    cex_sim(:,1) = a0(out_cex_idx);
    theta_avgn_sim(:,1) = a0(out_theta_avgn_idx);
    theta_avgp_sim(:,1) = a0(out_theta_avgp_idx);
    etan_sim(:,1) = a0(out_etan_idx);
    etap_sim(:,1) = a0(out_etap_idx);
    ce0n_sim(:,1) = a0(out_ce0n_idx);
    ce0p_sim(:,1) = a0(out_ce0p_idx);
    etas_sim(:,1) = a0(out_etas_idx);
    volt_sim(:,1) = a0(out_Volt_idx);
    nLis_sim(:,1) = a0(out_nLis_idx);
    nLie_sim(:,1) = a0(out_nLie_idx);

    % save state-dependent param
    Den0_sim(:,1) = p0(out_Den_idx);
    dDen0_sim(:,1) = p0(out_dDen_idx);
    Des0_sim(:,1) = p0(out_Des_idx);
    dDes0_sim(:,1) = p0(out_dDes_idx);
    Dep0_sim(:,1) = p0(out_Dep_idx);
    dDep0_sim(:,1) = p0(out_dDep_idx);

%     disp('---------------------------------')
%     disp('Running')
    
    % Simulate DFN & Sensitivity Equations
    try
        for k=1:(NT-1)
            if (mod(k,1000) == 0)
                fprintf('Iter:%d, v_sim: %f\n',k,v_sim(k));
            end

            Cur = I(k+1);

            % DFN
            Fk = F('x0',x_sim(:,k),'z0',z_sim(:,k),'p',[theta_0;Cur]);
            x_sim(:,k+1) = full(Fk.xf);
            z_sim(:,k+1) = full(Fk.zf);
            v_sim(k+1) = full(Fk.qf)/dt;

            f0 = full(f_out(theta_0,x_sim(:,k+1),z_sim(:,k+1),Cur));
            g0 = full(g_out(theta_0,x_sim(:,k+1),z_sim(:,k+1),Cur));
            a0 = full(alg_out(theta_0,x_sim(:,k+1),z_sim(:,k+1),Cur));
            p0 = full(par_out(theta_0,x_sim(:,k+1),z_sim(:,k+1),Cur));

            % save x
            csn_sim(:,k+1) = f0(out_csn_idx);
            csp_sim(:,k+1) = f0(out_csp_idx);
            ce_sim(:,k+1) = f0(out_ce_idx);
            cen_sim(:,k+1) = ce_sim(1:(p.Nxn-1),k+1);
            ces_sim(:,k+1) = ce_sim((p.Nxn-1 +1):(p.Nxn-1 + p.Nxs-1),k+1);
            cep_sim(:,k+1) = ce_sim((p.Nxn-1 + p.Nxs-1 +1):(p.Nxn-1 + p.Nxs-1 + p.Nxp-1),k+1);

            T1_sim(k+1) = f0(out_T_idx);
            T2_sim(k+1) = f0(end);

            % save z
            phi_s_n_sim(:,k+1) = g0(out_phisn_idx);
            phi_s_p_sim(:,k+1) = g0(out_phisp_idx);
            ien_sim(:,k+1) = g0(out_ien_idx);
            iep_sim(:,k+1) = g0(out_iep_idx);
            phie_sim(:,k+1) = g0(out_phie_idx);
            jn_sim(:,k+1) = g0(out_jn_idx);
            jp_sim(:,k+1) = g0(out_jp_idx);

            % save alg.
            cssn_sim(:,k+1) = a0(out_cssn_idx);
            cssp_sim(:,k+1) = a0(out_cssp_idx);
            cex_sim(:,k+1) = a0(out_cex_idx);
            theta_avgn_sim(:,k+1) = a0(out_theta_avgn_idx);
            theta_avgp_sim(:,k+1) = a0(out_theta_avgp_idx);
            etan_sim(:,k+1) = a0(out_etan_idx);
            etap_sim(:,k+1) = a0(out_etap_idx);
            ce0n_sim(:,k+1) = a0(out_ce0n_idx);
            ce0p_sim(:,k+1) = a0(out_ce0p_idx);
            etas_sim(:,k+1) = a0(out_etas_idx);
            volt_sim(:,k+1) = a0(out_Volt_idx);
            nLis_sim(:,k+1) = a0(out_nLis_idx);
            nLie_sim(:,k+1) = a0(out_nLie_idx);

            % save param.
            Den0_sim(:,k+1) = p0(out_Den_idx);
            dDen0_sim(:,k+1) = p0(out_dDen_idx);
            Des0_sim(:,k+1) = p0(out_Des_idx);
            dDes0_sim(:,k+1) = p0(out_dDes_idx);
            Dep0_sim(:,k+1) = p0(out_Dep_idx);
            dDep0_sim(:,k+1) = p0(out_dDep_idx);

            % Step SOC forward
            SOC(k+1) = (mean(c_avg_n(:,k+1)) - p.cn_low) / (p.cn_high - p.cn_low);

            if SensFlag == 1 % Calculate sensitivity flag == true
                % Sensitivity eqns
                Sk = SS('x0',S1_sim_blk(:,k),'z0',S2_sim_blk(:,k),'p',[theta_0;x_sim(:,k);z_sim(:,k);Cur]);
                S1_sim_blk(:,k+1) = full(Sk.xf);
                S2_sim_blk(:,k+1) = full(Sk.zf);
                S3_sim(:,k+1) = full(Sk.qf)/dt;
            end

             % Check voltage constraints
            if v_sim(k+1) <= p.volt_min
                fprintf('Min voltage is reached at %d iteration. Stopping simulation and concatenating data. \n',k);
                v_sim = concatenate_data(v_sim);
                sens(1:length(S3_sim),:) = S3_sim'; 
                varargout{1} = sens;
                state_info = [];
                return
            end

            if v_sim(k+1) >= p.volt_max
                fprintf('Max voltage is reached at %d iteration. Stopping simulation and concatenating data. \n',k);
                v_sim = concatenate_data(v_sim);
                sens(1:length(S3_sim),:) = S3_sim';
                varargout{1} = sens;
                state_info = [];
                return 
            end
        end
    catch e
        errorMessage = sprintf('%s',getReport( e, 'extended', 'hyperlinks', 'on' ))
        fprintf('CasADi error. Stopping simulation and concatenating data. \n');
        v_sim = concatenate_data(v_sim);
        sens(1:length(S3_sim),:) = S3_sim';
        varargout{1} = sens;
        state_info = [];
        return
    end

    %% Collect Outputs
    % save x info
    state_info.csn_sim = csn_sim;
    state_info.csp_sim = csp_sim;
    state_info.ce_sim = ce_sim;
    state_info.T1_sim = T1_sim;
    state_info.T2_sim = T2_sim;
    
    % Internal States
    state_info.cssn_sim = cssn_sim;
    state_info.cssp_sim = cssp_sim;
    state_info.cex_sim = cex_sim;
    state_info.theta_avgn_sim = theta_avgn_sim;
    state_info.theta_avgp_sim = theta_avgp_sim;    
    state_info.etan_sim = etan_sim;
    state_info.etap_sim = etap_sim;
    state_info.ce0n_sim = ce0n_sim;
    state_info.ce0p_sim = ce0p_sim;
    state_info.etas_sim = etas_sim;
%     state_info.phi_s_n_sim = phi_s_n_sim;
%     state_info.phi_s_p_sim = phi_s_p_sim;
%     state_info.ien_sim = ien_sim;
%     state_info.iep_sim = iep_sim;
%     state_info.phie_sim = phie_sim;
%     state_info.jn_sim = jn_sim;
%     state_info.jp_sim = jp_sim;
%     state_info.volt_sim = volt_sim;
%     state_info.nLis_sim = nLis_sim;
%     state_info.nLie_sim = nLie_sim;
%     state_info.SOC = SOC;
%     state_info.p = p;
%     state_info.Den0_sim = Den0_sim;
%     state_info.dDen0_sim = dDen0_sim;
%     state_info.Des0_sim = Des0_sim;
%     state_info.dDes0_sim = dDes0_sim;
%     state_info.Dep0_sim = Dep0_sim;
%     state_info.dDep0_sim = dDep0_sim;
%     state_info.cen_sim = cen_sim;
%     state_info.ces_sim = ces_sim;
%     state_info.cep_sim = cep_sim;

    %% Return
    if SensFlag == 1
        sens = S3_sim'; % Nt-by-Np matrix
        varargout{1} = sens;
    end
end

