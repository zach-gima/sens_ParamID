%%% CasADi implementation of spmet
% By: Zach Gima 2019-4-15

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

function [v_sim,state_info,varargout] = spmet_casadi(p,data,varargin)
    import casadi.*
      
    %% Parse Inputs
    I = data.cur;
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
        p = update_c_e_mats(p,SensFlag);
        
    else % just simulating the model
        theta_sx = [];
        theta_0 = [];
        SensFlag = 0;
        
        % Update parameter dependencies
        p = update_dependencies(p);
                
        % Electrolyte concentration, solid & elyte electric potential, and
        % elyte current matrices (function of other params)
        p = update_c_e_mats(p,SensFlag);
    end
    
    %% SPMeT State Initialization
    x0_nom = init_spmet(p,V0,T_amb,0);

    %% Initialize ode system setup
    % Define Model Variables
    x = SX.sym('x',size(x0_nom));
    
    %%% Input
    u = SX.sym('u');
    
    % Build the SPMeT ode System for the intial states; 
    % ode_spmet_casadi holds all of the odes that define the model
    [x_dot, x_outs, L, alg_states] = ode_spmet_casadi(x,u,p); % ZTG add 2019-6-17
    
    %% Setup Sensitivity Equations

    % Initial condition
    x0_call = Function('x0_call',{theta_sx},{x0_nom},{'params_sx'},{'result'});
    x0 = full(x0_call(theta_0)); % State I.C
    
    if SensFlag == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1. Calculate Jacobian automatically
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Jac_V0_param = jacobian(L,theta_sx);
        V0_call = Function('V0_call',{theta_sx,x,u},{Jac_V0_param},{'params_sx','x','u'},{'result'});
        S3_0_V0_params = full(V0_call(theta_0,x0,0));

        % Calculate Jacobian matrix (\frac{\partial f}{\partial x}}), A11
        f_x_jac = jacobian(x_dot,x);

        % Calculate Jacobian matrix (\frac{\partial f}{\partial \theta}}), B1
        f_theta_jac = jacobian(x_dot,theta_sx);

        % Calculate Jacobian matrix (\frac{\partial h}{\partial x}}), C
        h_x_jac = jacobian(L,x);

        % Calculate Jacobian matrix (\frac{\partial h}{\partial \theta}}), E
        h_theta_jac = jacobian(L,theta_sx);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. Build sensitivity equation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        S1_0 = zeros(size(x,1),np);
        S3_0 = S3_0_V0_params;

        % Symbolic Sensitivity

        %%% S1_dot
        S1_blk = SX.sym('S1_blk',size(x,1)*np,1);
        S1 = reshape(S1_blk,size(x,1),np);

        tmp = repmat({f_x_jac},np,1);
        f_x_blk = blkdiag(tmp{:});

        f_theta_blk = reshape(f_theta_jac,size(x,1)*np,1);

        S1_dot = f_x_blk*S1_blk + f_theta_blk;

        %%% S3
        S3 = h_x_jac * S1 + h_theta_jac;
        S3 = S3'; % change row-vector to column vector. 

        % IC for S1_blk
        S1_0_blk(:,1) = reshape(S1_0,size(x,1)*np,1);

        %%% Explicitly define sensitivity values needed for simulation
        S1_sim_blk(:,1) = S1_0_blk(:,1);
        
        sens = zeros(NT,np);
        S3_sim = zeros(np,NT); % necessary here regardless of sensflag
        S3_sim(:,1) = S3_0';
    end
    
    %% Integrator
    % Build ode System
    ode = struct('x', x, 'p', [theta_sx;u], 'ode', x_dot, 'quad', L);
    
    % Note: Play with values below; for the DFN, sensitivity wouldnt converge/work for higher tolerance (10^-6) but
    % the DFN did, so separate CasADi integrator functions made for the DFN (F)
    % and Sensitivity odes (SS)
    % opts = struct('fsens_err_con',true,'quad_err_con',true,'t0',0,'tf',dt);
    opts1 = struct('tf',dt, 'abstol',1e-6,'reltol',1e-6,'disable_internal_warnings',1);
    opts2 = struct('tf',dt, 'abstol', 1e-4, 'reltol', 1e-4,'disable_internal_warnings',1);

    F = integrator('F','cvodes', ode, opts1);

    if SensFlag == 1
        ode_s = struct('x', S1_blk, 'p', [theta_sx;x;u], 'ode', S1_dot, 'quad', S3);
        SS = integrator('SS', 'cvodes', ode_s, opts2);
    end
    
    %% Build function for simulated states % [SHP Change]
    f_out = Function('f_out',{theta_sx,x,u},{x_outs},{'params_sx','x','u'},{'f_out'});
    alg_out = Function('alg_out',{theta_sx,x,u},{alg_states},{'params_sx','x','u'},{'alg_out'}); % ZTG add 2019-6-17
    C_out = Function('C_out',{},{},{},{'C_out'});
    f0 = full(f_out(theta_0,x0,0));
    a0 = full(alg_out(theta_0,x0,0));

     %% Indexing

    % index for x
    out_csn_idx = 1 : (p.Nr-1); % 1:29
    out_csp_idx = p.Nr : 2*(p.Nr-1); % 30:58
    out_ce_idx = 2*p.Nr-1 :  2*p.Nr-1 + p.Nx-4; % 59:80
%     out_T1_idx = (end-2) of x; % 81
%     out_T2_idx = (end-1) of x;% 82
%     out_dsei_idx = end of x
    
    % index for other internl state outputs
    % alg_states = [c_ss_n; c_ss_p; c_ex; SOC_n; SOC_p; eta_n; eta_p; c_e0n; c_e0p; eta_s];
    out_cssn_idx = 1;
    out_cssp_idx = 2;
    out_cex_idx = (2+1) : (p.Nxn + p.Nxs + p.Nxp-1 + 4); % check this index
    out_SOCn_idx = ((p.Nxn + p.Nxs + p.Nxp-1 + 4) + 1);
    out_SOCp_idx = ((p.Nxn + p.Nxs + p.Nxp-1 + 4) + 1 + 1);
    out_etan_idx = ((p.Nxn + p.Nxs + p.Nxp-1 + 4) + 1 + 1 + 1);
    out_etap_idx = ((p.Nxn + p.Nxs + p.Nxp-1 + 4) + 1 + 1 + 1 + 1);
    out_ce0n_idx = ((p.Nxn + p.Nxs + p.Nxp-1 + 4) + 1 + 1 + 1 + 1 + 1);
    out_ce0p_idx = ((p.Nxn + p.Nxs + p.Nxp-1 + 4) + 1 + 1 + 1 + 1 + 1 + 1);
    out_etas_idx = ((p.Nxn + p.Nxs + p.Nxp-1 + 4) + 1 + 1 + 1 + 1 + 1 + 1 + 1);
%     out_Volt_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +4) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +4);
%     out_nLis_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +5) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +5);
%     out_nLie_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +6) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +6);

    %% Simulation
    % pre-allocate
    v_sim = zeros(NT,1);
    T1_sim = zeros(NT,1);
    T2_sim = zeros(NT,1);
    delta_sei_sim = zeros(NT,1);

    % Need to set these to value for error causing script to exit reasons
    if SensFlag == 0
        S3_sim = 0;
        sens = 0;
    end
    
    % pre-allocate other internal states
    cssn_sim = zeros(NT,1);
    cssp_sim = zeros(NT,1);
    SOCn_sim = zeros(NT,1);
    SOCp_sim = zeros(NT,1);
    etan_sim = zeros(NT,1);
    etap_sim = zeros(NT,1);
    ce0n_sim = zeros(NT,1);
    ce0p_sim = zeros(NT,1);
    etas_sim = zeros(NT,1);

    % initialize
    x_sim(:,1) = x0;
    v_sim(1) = V0;

    % initialize other states
    csn_sim(:,1) = f0(out_csn_idx);
    csp_sim(:,1) = f0(out_csp_idx);
    ce_sim(:,1) = f0(out_ce_idx);
    T1_sim(1) = f0(end-2);
    T2_sim(1) = f0(end-1);
    delta_sei_sim(1) = f0(end);
    % ZTG Note: Below states added
    cssn_sim(1) = a0(out_cssn_idx);
    cssp_sim(1) = a0(out_cssp_idx);
    cex_sim(:,1) = a0(out_cex_idx);
    SOCn_sim(1) = a0(out_SOCn_idx);
    SOCp_sim(1) = a0(out_SOCp_idx);
    etan_sim(1) = a0(out_etan_idx);
    etap_sim(1) = a0(out_etap_idx);
    ce0n_sim(1) = a0(out_ce0n_idx);
    ce0p_sim(1) = a0(out_ce0p_idx);
    etas_sim(1) = a0(out_etas_idx);
%     volt_sim(:,1) = a0(out_Volt_idx);
%     nLis_sim(:,1) = a0(out_nLis_idx);
%     nLie_sim(:,1) = a0(out_nLie_idx);
    
    % Simulate SPMeT & Sensitivity Equations
    try
        for k=1:(NT-1)
            if (mod(k,1000) == 0)
                fprintf('Iter:%d, v_sim: %f\n',k,v_sim(k));
            end

            Cur = I(k+1);

            % DFN
            Fk = F('x0',x_sim(:,k),'p',[theta_0;Cur]);
            x_sim(:,k+1) = full(Fk.xf);
            v_sim(k+1) = full(Fk.qf)/dt;

            f0 = full(f_out(theta_0,x_sim(:,k+1),Cur));
            a0 = full(alg_out(theta_0,x_sim(:,k+1),Cur));
                        
            % save x
            csn_sim(:,k+1) = f0(out_csn_idx);
            csp_sim(:,k+1) = f0(out_csp_idx);
            ce_sim(:,k+1) = f0(out_ce_idx);
            T1_sim(k+1) = f0(end-2);
            T2_sim(k+1) = f0(end-1);
            delta_sei_sim(k+1) = f0(end);
            
            % Save other internal states
            cssn_sim(k+1) = a0(out_cssn_idx);
            cssp_sim(k+1) = a0(out_cssp_idx);
            cex_sim(:,k+1) = a0(out_cex_idx);
            SOCn_sim(k+1) = a0(out_SOCn_idx);
            SOCp_sim(k+1) = a0(out_SOCp_idx);
            etan_sim(k+1) = a0(out_etan_idx);
            etap_sim(k+1) = a0(out_etap_idx);
            ce0n_sim(k+1) = a0(out_ce0n_idx);
            ce0p_sim(k+1) = a0(out_ce0p_idx);
            etas_sim(k+1) = a0(out_etas_idx);
%             volt_sim(:,1) = a0(out_Volt_idx);
%             nLis_sim(:,1) = a0(out_nLis_idx);
%             nLie_sim(:,1) = a0(out_nLie_idx);

            % Step SOC forward
        %     SOC(k+1) = (mean(c_avg_n(:,k+1)) - cn_low) / (cn_high - cn_low);

            if SensFlag == 1
                % Sensitivity eqns
                Sk = SS('x0',S1_sim_blk(:,k),'p',[theta_0;x_sim(:,k);Cur]);
                S1_sim_blk(:,k+1) = full(Sk.xf);
                S3_sim(:,k+1) = full(Sk.qf)/dt;
            end

            % Check voltage constraints
            if v_sim(k+1) <= p.volt_min
                fprintf('Min voltage is reached at %d iteration. Stopping simulation and concatenating data. \n',k);
                v_sim = concatenate_data(v_sim);
                state_info = [];
                sens(1:length(S3_sim),:) = S3_sim'; 
                varargout{1} = sens;
                return
            end

            if v_sim(k+1) >= p.volt_max
                fprintf('Max voltage is reached at %d iteration. Stopping simulation and concatenating data. \n',k);
                v_sim = concatenate_data(v_sim);
                state_info = [];
                sens(1:length(S3_sim),:) = S3_sim';
                varargout{1} = sens;
                return 
            end
        end
    catch e
        errorMessage = sprintf('%s',getReport( e, 'extended', 'hyperlinks', 'on' ))
        fprintf('CasADi error. Stopping simulation and concatenating data. \n');
        v_sim = concatenate_data(v_sim);
        state_info = [];
        sens(1:length(S3_sim),:) = S3_sim';   
        varargout{1} = sens;
        return
    end

    %% Collect Outputs
    if SensFlag == 1
        % Provides sensitivity data / Jacobian to fmincon
        sens = S3_sim';
        varargout{1} = sens;
    end
    
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
    state_info.SOCn_sim = SOCn_sim;
    state_info.SOCp_sim = SOCp_sim;
    state_info.etan_sim = etan_sim;
    state_info.etap_sim = etap_sim;
    state_info.ce0n_sim = ce0n_sim;
    state_info.ce0p_sim = ce0p_sim;
    state_info.etas_sim = etas_sim;
    %     state_info.ien_sim = ien_sim;
    %     state_info.iep_sim = iep_sim;
    %     state_info.phie_sim = phie_sim;
    %     state_info.jn_sim = jn_sim;
    %     state_info.jp_sim = jp_sim;
    %     state_info.cex_sim = cex_sim;
    %     state_info.theta_avgn_sim = theta_avgn_sim;
    %     state_info.theta_avgp_sim = theta_avgp_sim;
    % state_info.etan_sim = etan_sim;
    % state_info.etap_sim = etap_sim;
    %     state_info.ce0n_sim = ce0n_sim;
    %     state_info.ce0p_sim = ce0p_sim;
    %     state_info.etasLn_sim = etasLn_sim;
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
end
    
