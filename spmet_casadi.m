%%% CasADi implementation of spmet
% By: Zach Gima 2019-4-15

% INPUTS --- 
% (1) data: struct with fields cur (current), time, V0
% (2) x0_nom: initial nominal (i.e. not symbolic casadi variables) values associated with the states of the
% model 
% (3) theta_0: initial numerical values for the parameters we'd like to
% identify
% (4) theta_sx: casadi format of the parameters we'd like to identify
% (5) p: parameter vector

% OUTPUTS --- 
% (1) voltage 
% (2) sensitivity of parameters (p) passed into the
% function (partial V / partial theta) 
% (3) algebraic (may be incorrect description) states of the model

function [v_sim,sens,alg_states] = spmet_casadi(data,theta_0,theta_sx,p,SensFlag)
    import casadi.*
    
    %% Parse Parameters
    np = length(theta_0);
    I = data.cur;
    t = data.time;
    V0 = data.V0;
    dt =  t(2)-t(1);
    NT = length(t);
    T_amb = data.T_amb;
    
    %% update parameter dependencies (some parameters are functions of other parameters which may have changed)
    p = update_dependencies(p);
    
    %% SPMeT State Initialization
    x0_nom = init_spmet(p,V0,T_amb,0);

    %% Initialize ode system setup
    % Define Model Variables
    x = SX.sym('x',size(x0_nom));
    
    %%% Input
    u = SX.sym('u');
    
    % Build the SPMeT ode System; 
    % ode_spmet_casadi holds all of the odes that define the model
    [x_dot, x_outs, L, ~] = ode_spmet_casadi(x,u,p);
    
    %% Setup Sensitivity Equations

    % Initial condition
    x0_call = Function('x0_call',{theta_sx},{x0_nom},{'params_sx'},{'result'});
    x0_init = full(x0_call(theta_0));
    
    % State I.C
    x0 = full(x0_call(theta_0));
    
    if SensFlag == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1. Calculate Jacobian automatically
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Jac_V0_param = jacobian(L,theta_sx);
        V0_call = Function('V0_call',{theta_sx,x,u},{Jac_V0_param},{'params_sx','x','u'},{'result'});
        S3_0_V0_params = full(V0_call(theta_0,x0_init,0));

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
    f0 = full(f_out(theta_0,x0,0));

     %% Indexing

    % index for x
    out_csn_idx = 1 : (p.Nr-1); % 1:29
    out_csp_idx = p.Nr : 2*(p.Nr-1); % 30:58
    out_cex_idx = 2*p.Nr-1 :  2*p.Nr-1 + p.Nx-4; % 59:80
%     out_T1_idx = (end-2) of x; % 81
%     out_T2_idx = (end-1) of x;% 82
%     out_dsei_idx = end of x
    
    % index for other outputs
%     out_cssn_idx = 1:p.Nxn-1 ;
%     out_cssp_idx = (p.Nxn-1+1) : (p.Nxn-1 + p.Nxp-1);
%     out_cex_idx = (p.Nxn-1 + p.Nxp-1 +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4));
%     out_theta_avgn_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1));
%     out_theta_avgp_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1));
%     out_etan_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1));
%     out_etap_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) + (p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1));
%     out_ce0n_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +1);
%     out_ce0p_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +2) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +2);
%     out_etasLn_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +3) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +3);
%     out_Volt_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +4) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +4);
%     out_nLis_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +5) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +5);
%     out_nLie_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +6) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +6);

    %% Simulation
    % pre-allocate
    v_sim = zeros(NT,1);
    T1_sim = zeros(NT,1);
    T2_sim = zeros(NT,1);
    delta_sei_sim = zeros(NT,1);
    sens = zeros(NT,np);
    S3_sim = zeros(np,NT);

    % initialize
    x_sim(:,1) = x0;
    v_sim(1) = V0;

    % save x
    csn_sim(:,1) = f0(out_csn_idx);
    csp_sim(:,1) = f0(out_csp_idx);
    cex_sim(:,1) = f0(out_cex_idx);
    T1_sim(1) = f0(end-2);
    T2_sim(1) = f0(end-1);
    delta_sei_sim(1) = f0(end);
    
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

            % save x
            csn_sim(:,k+1) = f0(out_csn_idx);
            csp_sim(:,k+1) = f0(out_csp_idx);
            cex_sim(:,k+1) = f0(out_cex_idx);
            T1_sim(k+1) = f0(end-2);
            T2_sim(k+1) = f0(end-1);
            delta_sei_sim(k+1) = f0(end);

            % Step SOC forward
        %     SOC(k+1) = (mean(c_avg_n(:,k+1)) - cn_low) / (cn_high - cn_low);

            if SensFlag == 1
                % Sensitivity eqns
                Sk = SS('x0',S1_sim_blk(:,k),'p',[theta_0;x_sim(:,k);Cur]);
                S1_sim_blk(:,k+1) = full(Sk.xf);
                S3_sim(:,k+1) = full(Sk.qf)/dt;
            end

            % Check voltage constraints and exit if violated
            if v_sim(k+1) <= p.volt_min
                fprintf('Min voltage is reached at %d iteration. Stopping simulation and concatenating data. \n',k);
                v_sim = concatenate_data(v_sim);
                sens(1:length(S3_sim),:) = S3_sim'; 
                alg_states = [];
                return
            end

            if v_sim(k+1) >= p.volt_max
                fprintf('Max voltage is reached at %d iteration. Stopping simulation and concatenating data. \n',k);
                v_sim = concatenate_data(v_sim);
                sens(1:length(S3_sim),:) = S3_sim';
                alg_states = [];
                return 
            end
        end
    catch e
        errorMessage = sprintf('%s',getReport( e, 'extended', 'hyperlinks', 'on' ))
        fprintf('CasADi error. Stopping simulation and concatenating data. \n');
        v_sim = concatenate_data(v_sim);
        sens(1:length(S3_sim),:) = S3_sim';   
        
        alg_states.csn_sim = csn_sim;
        alg_states.csp_sim = csp_sim;
        alg_states.cex_sim = cex_sim;
        alg_states.T1_sim = T1_sim;
        alg_states.T2_sim = T2_sim;
        
        save('casadi_debug.mat','v_sim','alg_states','sens','I', 't')
        return
    end


    %% Debug
%     Voltage = data.V;
% 
%     %%% Compare DFN voltage to SPMeT
%     figure
%     hold on
%     plot(v_sim(1:3500))
%     plot(Voltage(1:3500))
%     legend('DFN','SPMeT')
%     fs = 25;
%     set(gca,'Fontsize',fs)
%     hold off
% 
%     %%% plot whole v_sim profile
%     figure;plot(v_sim)
%     title('Entire simulated SPMeT voltage profile')
%     fs = 25;
%     set(gca,'Fontsize',fs)  
% 
%     %%% plot current
%     figure;plot(I)
%     title('Current profile')
%     fs = 25;
%     set(gca,'Fontsize',fs)   
% 
%     %%% plot c_s_n
%     figure;plot(csn_sim(end,:))
%     title('c_s^-')
%     fs = 25;
%     set(gca,'Fontsize',fs) 
%     
%     %%% plot c_s_n
%     figure;plot(csp_sim(end,:))
%     title('c_s^+')
%     fs = 25;
%     set(gca,'Fontsize',fs) 
%     
%     %%% plot_c_e
%     figure;plot(cex_sim(end,:))
%     title('c_e')
%     fs = 25;
%     set(gca,'Fontsize',fs) 
% 
%     figure;plot(T1_sim)
%     title('T_1')
%     fs = 25;
%     set(gca,'Fontsize',fs) 
%     
%     figure;plot(T2_sim)
%     title('T_2')
%     fs = 25;
%     set(gca,'Fontsize',fs) 
    
    %% Collect Outputs
    sens = S3_sim';
    
    alg_states.csn_sim = csn_sim;
    alg_states.csp_sim = csp_sim;
    alg_states.cex_sim = cex_sim;
    alg_states.T1_sim = T1_sim;
    alg_states.T2_sim = T2_sim;
    
    %     alg_states.phi_s_n_sim = phi_s_n_sim;
    %     alg_states.phi_s_p_sim = phi_s_p_sim;
    %     alg_states.ien_sim = ien_sim;
    %     alg_states.iep_sim = iep_sim;
    %     alg_states.phie_sim = phie_sim;
    %     alg_states.jn_sim = jn_sim;
    %     alg_states.jp_sim = jp_sim;
    % alg_states.cssn_sim = cssn_sim;
    % alg_states.cssp_sim = cssp_sim;
    %     alg_states.cex_sim = cex_sim;
    %     alg_states.theta_avgn_sim = theta_avgn_sim;
    %     alg_states.theta_avgp_sim = theta_avgp_sim;
    % alg_states.etan_sim = etan_sim;
    % alg_states.etap_sim = etap_sim;
    %     alg_states.ce0n_sim = ce0n_sim;
    %     alg_states.ce0p_sim = ce0p_sim;
    %     alg_states.etasLn_sim = etasLn_sim;
    %     alg_states.volt_sim = volt_sim;
    %     alg_states.nLis_sim = nLis_sim;
    %     alg_states.nLie_sim = nLie_sim;
    %     alg_states.SOC = SOC;
    %     alg_states.p = p;
    %     alg_states.Den0_sim = Den0_sim;
    %     alg_states.dDen0_sim = dDen0_sim;
    %     alg_states.Des0_sim = Des0_sim;
    %     alg_states.dDes0_sim = dDes0_sim;
    %     alg_states.Dep0_sim = Dep0_sim;
    %     alg_states.dDep0_sim = dDep0_sim;
    %     alg_states.cen_sim = cen_sim;
    %     alg_states.ces_sim = ces_sim;
    %     alg_states.cep_sim = cep_sim;

    %% Use below to provide Jacobian to fmincon
%     varargout{1} = S3;
end
    
