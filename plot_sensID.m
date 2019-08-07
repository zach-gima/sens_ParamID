%% Plot Sens_ID

% function plot_sensID(ID_out)
function plot_sensID(metrics,results_folder,baseline,date)

    fs = 25;
    run param/params_bounds
    
    %% Parse inputs
    cssn_rmse = metrics.cssn_rmse;
    cssp_rmse = metrics.cssp_rmse;
    etas_rmse = metrics.etas_rmse;
    ce0n_rmse = metrics.ce0n_rmse;
    ce0p_rmse = metrics.ce0p_rmse;
    norm_param_dist = metrics.norm_param_dist;
%     fmincon_iter = metrics.fmincon_iter;
%     cost_evolution = metrics.cost_evolution;
    rmse_vec = metrics.rmse_vec;

    theta_iter_vec = metrics.theta_iter_vec;
    per_param_norm_error = metrics.per_param_norm_error;
    
    %% plot
    
    %rmse
    figure('Position', [100 100 900 700])
 
    yyaxis left
    plot(theta_iter_vec,cssn_rmse,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k')
    ylabel('c_{ss}^- RMSE (mol/m^3)','Interpreter','tex')
    yyaxis right
    plot(theta_iter_vec,cssp_rmse,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k')
    ylabel('c_{ss}^+ RMSE (mol/m^3)','Interpreter','tex')
    
    legend('Surface Conc. (Anode)','Surface Conc. (Cathode)')
    xlabel('Batch')
    xticks(theta_iter_vec)
    title('Anodic & Cathodic Surface Concentration RMSE')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(results_folder,date,'_',baseline,'_','css_rmse.fig'))
    print(strcat(results_folder,date,'_',baseline,'_','css_rmse'),'-dpng')

%     %anode fit
%     figure('Position', [100 100 900 700])
%     
%     plot(t,cssn_true,t,cssn_initial,t,cssn_final,'LineWidth',3)
%     legend('Truth','Before ParamID', 'After ParamID')
%     xlabel('Time (s)')
%     ylabel('c_{ss}^- (mol/m^3)','Interpreter','tex')
%     title('Anodic Surface Concentration Evolution')
%     set(gca,'FontSize',fs)
%     box on
%     grid on
% 
%     savefig(strcat(plots_folder,'cssn_fit.fig'))
%     print(strcat(plots_folder,'cssn_fit'),'-dpng')
%     
%     %cathode fit
%     figure('Position', [100 100 900 700])
%     plot(t,cssp_true,t,cssp_initial,t,cssp_final,'LineWidth',3)
%     legend('Truth','Before ParamID', 'After ParamID')
%     xlabel('Time (s)','Fontsize',fs)
%     ylabel('c_{ss}^{+} (mol/m^3)','Fontsize',fs);%,'Interpreter','latex')
%     title('Cathodic Surface Concentration Evolution')
%     set(gca,'FontSize',fs)
%     box on
%     grid on
% 
%     savefig(strcat(plots_folder,'cssp_fit.fig'))
%     print(strcat(plots_folder,'cssp_fit'),'-dpng')
    
    
    %rmse
    figure('Position', [100 100 900 700])
    plot(theta_iter_vec,etas_rmse,'-o','MarkerSize',10,'LineWidth',3,'MarkerEdgeColor','k')
    xlabel('Batch')
    xticks(theta_iter_vec)
    ylabel('\eta_s RMSE (V)','Fontsize',fs);%,'Interpreter','latex')
    title('Side Reaction Overpotential RMSE')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(results_folder,date,'_',baseline,'_','etas_rmse.fig'))
    print(strcat(results_folder,date,'_',baseline,'_','etas_rmse'),'-dpng')

%     % etas fit
%     figure('Position', [100 100 900 700])
%     plot(t,etas_true,t,etas_initial,t,etas_final,'LineWidth',3)
%     legend('Truth','Before ParamID', 'After ParamID')
%     xlabel('Time (s)','Fontsize',fs)
%     ylabel('\eta_s (V)','Fontsize',fs)%,'Interpreter','latex')
%     title('Side Reaction Overpotential Evolution')
%     set(gca,'FontSize',fs)
%     box on
%     grid on
% 
%     savefig(strcat(plots_folder,'etas_fit.fig'))
%     print(strcat(plots_folder,'etas_fit'),'-dpng')
    
    %rmse
    figure('Position', [100 100 900 700])
    plot(theta_iter_vec,ce0n_rmse,'-o',theta_iter_vec,ce0p_rmse,'-o','MarkerSize',10,'LineWidth',3,'MarkerEdgeColor','k')
    legend('Elyte Conc. at CC (Anode)','Elyte Conc. at CC (Cathode)')
    xlabel('Batch')
    xticks(theta_iter_vec)
    ylabel('c_e(0^{\pm}) RMSE (mol/m^3)','Fontsize',fs);%,'Interpreter','latex')
    title('Elyte Conc. at CC RMSE')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(results_folder,date,'_',baseline,'_','ce0_rmse.fig'))
    print(strcat(results_folder,date,'_',baseline,'_','ce0_rmse'),'-dpng')
    
%     %anode
%     figure('Position', [100 100 900 700])
%     plot(t,ce0n_true,t,ce0n_initial,t,ce0n_final,'LineWidth',3)
%     legend('Truth','Before ParamID', 'After ParamID')
%     xlabel('Time (s)','Fontsize',fs)
%     ylabel('c_e(0^{-}) (mol/m^3)','Fontsize',fs);%,'Interpreter','latex')
%     title('Anodic Elyte Conc. at CC Evolution')
%     set(gca,'FontSize',fs)
%     box on
%     grid on
% 
%     savefig(strcat(plots_folder,'ce0n_fit.fig'))
%     print(strcat(plots_folder,'ce0n_fit'),'-dpng')
%     
%     %cathode
%     figure('Position', [100 100 900 700])
%     plot(t,ce0p_true,t,ce0p_initial,t,ce0p_final,'LineWidth',3)
%     legend('Truth','Before ParamID', 'After ParamID')
%     xlabel('Time (s)','Fontsize',fs)
%     ylabel('c_e(0^{+}) (mol/m^3)','Fontsize',fs);%,'Interpreter','latex')
%     title('Cathodic Elyte Conc. at CC Evolution')
%     set(gca,'FontSize',fs)
%     box on
%     grid on
% 
%     savefig(strcat(plots_folder,'ce0p_fit.fig'))
%     print(strcat(plots_folder,'ce0p_fit'),'-dpng')
    
    %% Plot Voltage Profile before and after ID    
%     figure('Position', [100 100 900 700])
%     plot(t,V_true,t,V_sim_initial,t,V_sim_final,'LineWidth',3)
%     legend('Truth Data','Initial Fit','Final Fit')
%     xlabel('Time (s)')
%     ylabel('Voltage (V)')
%     title('Voltage Fit')
%     set(gca,'FontSize',fs)
%     box on
%     grid on
% 
%     savefig(strcat(plots_folder,'voltage_fit.fig'))
%     print(strcat(plots_folder,'voltage_fit'),'-dpng')
    
    %% Cost Function / Voltage RMSE
    figure('Position', [100 100 900 700])
    plot(theta_iter_vec,rmse_vec,'-o','LineWidth',3, 'MarkerSize',10,'MarkerEdgeColor','k')
    xlabel('Batch')
    xticks(theta_iter_vec)
    ylabel('Voltage RMSE (V)')
    title('Voltage RMSE vs. Batch')
    set(gca,'Fontsize',fs)
    box on
    grid on
     
    savefig(strcat(results_folder,date,'_',baseline,'_','voltage_rmse.fig'))
    print(strcat(results_folder,date,'_',baseline,'_','voltage_rmse'),'-dpng')
    
    
    %% Plot RMSE vs Computational Time
%     figure('Position', [100 100 900 700])
%     plot(wallclock, rmse_vec,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k');
%     title('RMSE vs. Wall Clock time')
%     xlabel('Wall Clock time (Hr)')
%     ylabel('RMSE (V)')
%     set(gca,'Fontsize',fs)
%     box on
%     grid on
% 
%     savefig(strcat(plots_folder,'rmse_evolution.fig'))
%     print(strcat(plots_folder,'rmse_evolution'),'-dpng')

    %% Plot Normalized Distance between True and Estimated Parameters & percentage error
    figure('Position', [100 100 900 700])
    plot(theta_iter_vec, norm_param_dist,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k');
    title('Normalized Parameter Guess Error')
    xlabel('Batch')
    xticks(theta_iter_vec)
    ylabel('$|| \theta^* - \hat\theta ||$','Interpreter','Latex')  
    set(gca,'Fontsize',fs)

    box on
    grid on

    savefig(strcat(results_folder,date,'_',baseline,'_','norm_param_dist.fig'))
    print(strcat(results_folder,date,'_',baseline,'_','norm_param_dist'),'-dpng')
    
    % grouped bar chart has 1 group per row in the y_data 
    param_err = [];
    for zz = 1:length(per_param_norm_error)
       param_err = horzcat(param_err,per_param_norm_error{zz}); 
    end
    param_err = param_err';
    
    % bpcombine2(:), bpcombine3(:)];
    figure('Position', [100 100 900 700])
    hb = bar(theta_iter_vec, param_err,'grouped');
    ylabel('Per Parameter Error (%)')
    xlabel('Batch')
    xticks(theta_iter_vec)
    legend('R_s^+','ElecFactorDA',horzcat(char(949),'_e^-'),'t_+','R_f^-','R_f^+','Location', 'northoutside','Orientation','horizontal');
    set(gca,'Fontsize',fs)

    box on
    grid on
    
    savefig(strcat(results_folder,date,'_',baseline,'_','per_param_err.fig'))
    print(strcat(results_folder,date,'_',baseline,'_','per_param_err'),'-dpng')
end