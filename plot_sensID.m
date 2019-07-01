%% Plot Sens_ID

function plot_sensID(ID_out,bounds)
    
    fs = 25;
    
    output_folder = 'plots/';
    
    %% Parse Data
    num_events = length(ID_out.data);
    data = ID_out.data;
%     opt_data = data(opt_event_idx);
        
    states_initial = ID_out.states_initial;
    states_final = ID_out.states_final;
    
    theta = ID_out.theta;
    num_params = size(theta.history,1);
    %     params_final_idx = cell2mat(ID_out.params_final_idx);    
    theta_iter = size(theta.history,2);
    theta_iter_vec = 1:theta_iter;
    
    cost_evolution = ID_out.fmincon_history.fval;
    rmse_vec = [cost_evolution(1);cost_evolution(end)];
    fmincon_iter = 1:length(cost_evolution);
    
    wallclock = [0;ID_out.wallclock]/3600;
    
    t_end = 0; % for concatenating the time vector data
    
    cssn_true = [];
    cssp_true = [];
    etas_true =  [];
    ce0n_true =  [];
    ce0p_true =  [];
    
    for ii = 1:num_events
        t_cell{ii,1} = data(ii).time + t_end;
        V_true_cell{ii,1} = data(ii).V_exp;
        states_true{ii,1} = data(ii).states_true;

        
        cssn_true = vertcat(cssn_true,states_true{ii,1}.cssn_sim);
        cssp_true= vertcat(cssp_true,states_true{ii,1}.cssp_sim);
        etas_true = vertcat(etas_true,states_true{ii,1}.etas_sim);
        ce0n_true = vertcat(ce0n_true,states_true{ii,1}.ce0n_sim);
        ce0p_true = vertcat(ce0p_true,states_true{ii,1}.ce0p_sim);
    
        t_end = t_end + data(ii).time(end);
    end
    
    t = cell2mat(t_cell);
    V_true = cell2mat(V_true_cell);

    V_sim_initial = cell2mat(ID_out.V_sim_initial);
    V_sim_final = cell2mat(ID_out.V_sim_final);
    
    %% ParamID Dartboard
    
    %% Internal State RMSE
    
    %%% Pre-allocate states of interest
    % Surface concentration
    cssn_initial = [];
    cssn_final = [];
    cssp_initial = [];
    cssp_final = [];
    % Side rxn overpotential
    etas_initial = [];
    etas_final = [];    
    
    %Electrode overpotential
%     etan_initial = [];
%     etan_final = [];
%     etap_initial  = [];
%     etap_final  = [];
    % elyte conc. at current collectors
    ce0n_initial = [];
    ce0n_final = [];
    ce0p_initial = [];
    ce0p_final = [];

    for jj = 1:num_events
        cssn_initial = vertcat(cssn_initial,states_initial{jj}.cssn_sim);
        cssn_final = vertcat(cssn_final,states_final{jj}.cssn_sim);
        cssp_initial = vertcat(cssp_initial,states_initial{jj}.cssp_sim);
        cssp_final = vertcat(cssp_final,states_final{jj}.cssp_sim);
        
%         etan_initial = vertcat(etan_initial,states_initial{jj}.etan_sim);
%         etan_final = vertcat(etan_final,states_final{jj}.etan_sim);
%         etap_initial = vertcat(etap_initial,states_initial{jj}.etap_sim);
%         etap_final = vertcat(etap_final,states_final{jj}.etap_sim);
        
        etas_initial = vertcat(etas_initial,states_initial{jj}.etas_sim);
        etas_final = vertcat(etas_final,states_final{jj}.etas_sim);
        
        ce0n_initial = vertcat(ce0n_initial,states_initial{jj}.ce0n_sim);
        ce0n_final = vertcat(ce0n_final,states_final{jj}.ce0n_sim);
        ce0p_initial = vertcat(ce0p_initial,states_initial{jj}.ce0p_sim);
        ce0p_final = vertcat(ce0p_final,states_final{jj}.ce0p_sim);
    end
    
    %%% css RMSE and fits
    cssn_rmse(1) = rmse(cssn_true,cssn_initial);
    cssn_rmse(2) = rmse(cssn_true,cssn_final);
    cssp_rmse(1) = rmse(cssp_true,cssp_initial);
    cssp_rmse(2) = rmse(cssp_true,cssp_final);
    
    %rmse
    figure('Position', [100 100 900 700])
    hold on
    plot(theta_iter_vec,cssn_rmse,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k')
    plot(theta_iter_vec,cssp_rmse,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k')
    hold off
    legend('Surface Conc. (Anode)','Surface Conc. (Cathode)')
    xlabel('Iteration','Fontsize',fs)
    ylabel('c_{ss}{\pm} RMSE (mol/m^3)','Fontsize',fs)%,'Interpreter','Latex')
    title('Anodic & Cathodic Surface Concentration RMSE')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(output_folder,'css_rmse.fig'))
    print(strcat(output_folder,'css_rmse'),'-dpng')

    %anode fit
    figure('Position', [100 100 900 700])
    plot(t,cssn_true,t,cssn_initial,t,cssn_final,'LineWidth',3)
    legend('Truth','Before ParamID', 'After ParamID')
    xlabel('Time (s)','Fontsize',fs)
    ylabel('c_{ss}^{-} (mol/m^3)','Fontsize',fs);%,'Interpreter','latex')
    title('Anodic Surface Concentration Evolution')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(output_folder,'cssn_fit.fig'))
    print(strcat(output_folder,'cssn_fit'),'-dpng')
    
    %cathode fit
    figure('Position', [100 100 900 700])
    plot(t,cssp_true,t,cssp_initial,t,cssp_final,'LineWidth',3)
    legend('Truth','Before ParamID', 'After ParamID')
    xlabel('Time (s)','Fontsize',fs)
    ylabel('c_{ss}^{+} (mol/m^3)','Fontsize',fs);%,'Interpreter','latex')
    title('Cathodic Surface Concentration Evolution')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(output_folder,'cssp_fit.fig'))
    print(strcat(output_folder,'cssp_fit'),'-dpng')
    
    %%% electrode overpotential RMSE and fits    
%     etan_rmse(1) = rmse(etan_true,etan_initial);
%     etan_rmse(2) = rmse(etan_true,etan_final);
%     etap_rmse(1) = rmse(etap_true,etap_initial);
%     etap_rmse(2) = rmse(etap_true,etap_final);

    
    %%% sid rxn overpotential RMSE and fits
    etas_rmse(1) = rmse(etas_true,etas_initial);
    etas_rmse(2) = rmse(etas_true,etas_final);
    
    %rmse
    figure('Position', [100 100 900 700])
    plot(theta_iter_vec,etas_rmse,'-o','MarkerSize',10,'LineWidth',3,'MarkerEdgeColor','k')
    xlabel('Iteration','Fontsize',fs)
    ylabel('\eta_s RMSE (V)','Fontsize',fs);%,'Interpreter','latex')
    title('Side Reaction Overpotential RMSE')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(output_folder,'etas_rmse.fig'))
    print(strcat(output_folder,'etas_rmse'),'-dpng')

    % etas fit
    figure('Position', [100 100 900 700])
    plot(t,etas_true,t,etas_initial,t,etas_final,'LineWidth',3)
    legend('Truth','Before ParamID', 'After ParamID')
    xlabel('Time (s)','Fontsize',fs)
    ylabel('\eta_s (V)','Fontsize',fs)%,'Interpreter','latex')
    title('Side Reaction Overpotential Evolution')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(output_folder,'etas_fit.fig'))
    print(strcat(output_folder,'etas_fit'),'-dpng')
    
    %%% electrode concentration @ current collectors
    ce0n_rmse(1) = rmse(ce0n_true,ce0n_initial);
    ce0n_rmse(2) = rmse(ce0n_true,ce0n_final);
    ce0p_rmse(1) = rmse(ce0p_true,ce0p_initial);
    ce0p_rmse(2) = rmse(ce0p_true,ce0p_final);
    
    %rmse
    figure('Position', [100 100 900 700])
    plot(theta_iter_vec,ce0n_rmse,theta_iter_vec,ce0p_rmse,'LineWidth',3)
    legend('Elyte Conc. at CC (Anode)','Elyte Conc. at CC (Cathode)')
    xlabel('Iteration','Fontsize',fs)
    ylabel('c_e(0^{\pm}) RMSE (mol/m^3)','Fontsize',fs);%,'Interpreter','latex')
    title('Elyte Conc. at CC RMSE')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(output_folder,'ce0_rmse.fig'))
    print(strcat(output_folder,'ce0_rmse'),'-dpng')
    
    %anode
    figure('Position', [100 100 900 700])
    plot(t,ce0n_true,t,ce0n_initial,t,ce0n_final,'LineWidth',3)
    legend('Truth','Before ParamID', 'After ParamID')
    xlabel('Time (s)','Fontsize',fs)
    ylabel('c_e(0^{-}) (mol/m^3)','Fontsize',fs);%,'Interpreter','latex')
    title('Anodic Elyte Conc. at CC Evolution')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(output_folder,'ce0n_fit.fig'))
    print(strcat(output_folder,'ce0n_fit'),'-dpng')
    
    %cathode
    figure('Position', [100 100 900 700])
    plot(t,ce0p_true,t,ce0p_initial,t,ce0p_final,'LineWidth',3)
    legend('Truth','Before ParamID', 'After ParamID')
    xlabel('Time (s)','Fontsize',fs)
    ylabel('c_e(0^{+}) (mol/m^3)','Fontsize',fs);%,'Interpreter','latex')
    title('Cathodic Elyte Conc. at CC Evolution')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(output_folder,'ce0p_fit.fig'))
    print(strcat(output_folder,'ce0p_fit'),'-dpng')
    
    %% Plot Voltage Profile before and after ID    
    figure('Position', [100 100 900 700])
    plot(t,V_true,t,V_sim_initial,t,V_sim_final,'LineWidth',3)
    legend('Truth Data','Initial Fit','Final Fit')
    xlabel('Time (s)','Fontsize',fs)
    ylabel('Voltage (V)','Fontsize',fs)
    title('Voltage Fit')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(output_folder,'voltage_fit.fig'))
    print(strcat(output_folder,'voltage_fit'),'-dpng')
    
    %% Cost Function / Voltage RMSE
    figure('Position', [100 100 900 700])
    plot(fmincon_iter,cost_evolution,'-o','LineWidth',3, 'MarkerSize',10,'MarkerEdgeColor','k')
    xlabel('Iteration','Fontsize',fs)
    ylabel('Voltage RMSE (V)','Fontsize',fs)
    title('Voltage RMSE vs. Iter')
    set(gca,'Fontsize',fs)
    box on
    grid on
     
    savefig(strcat(output_folder,'voltage_rmse.fig'))
    print(strcat(output_folder,'voltage_rmse'),'-dpng')
    
    
    %% Plot RMSE vs Computational Time
    figure('Position', [100 100 900 700])
    plot(wallclock, rmse_vec,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k');
    title('RMSE vs. Wall Clock time')
    xlabel('Wall Clock time (Hr)')
    ylabel('RMSE (V)')
    set(gca,'Fontsize',fs)
    box on
    grid on

    savefig(strcat(output_folder,'rmse_evolution.fig'))
    print(strcat(output_folder,'rmse_evolution'),'-dpng')

    %% Plot Normalized Distance between True and Estimated Parameters
    norm_param_dist = zeros(theta_iter,1);
%     norm_truth_param = origin_to_norm(theta.truth,bounds,params_final_idx);
    norm_truth_param = origin_to_norm(theta.truth,bounds);

    for ii = 1:theta_iter
%         norm_ID_param = origin_to_norm(theta.history(:,ii),bounds,params_final_idx);
        norm_ID_param = origin_to_norm(theta.history(:,ii),bounds);
        norm_param_dist(ii) = norm(norm_truth_param-norm_ID_param);
    end
    
    figure('Position', [100 100 900 700])
    plot(theta_iter_vec, norm_param_dist,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k');
    title('Normalized Parameter Guess Error')
    xlabel('Iteration','FontSize',fs)
    ylabel('$|| \theta^* - \hat\theta ||$','Interpreter','Latex','FontSize',fs)  
%     xlim([1,iter(end)]);
    set(gca,'Fontsize',fs)%,'XTick',[1 : 1 : iter(end)])

    box on
    grid on

    savefig(strcat(output_folder,'norm_param_dist.fig'))
    print(strcat(output_folder,'norm_param_dist'),'-dpng')
end