%% For visualizing model mismatch between spmet and dfn
% Voltage Rmse
figure('Position', [100 100 900 700])
plot(V_true,'LineWidth', 2.5);
hold on
plot(V_sim_initial,'LineWidth', 2.5);
hold off
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Voltage')
legend('DFN','SPMeT')
set(gca,'Fontsize',fs)

v_rmse = rmse(V_true,V_sim_initial);

% etas 
figure('Position', [100 100 900 700])
plot(etas_true,'LineWidth', 2.5);
hold on
plot(etas_initial+0.4,'LineWidth', 2.5);
hold off
xlabel('Time (s)')
ylabel('\eta_s (V)','Fontsize',fs);%,'Interpreter','latex')
title('Side Reaction Overpotential')
legend('DFN','SPMeT')
set(gca,'Fontsize',fs)

etas_rmse = rmse(etas_true,etas_initial);
etas_rmse_fixed = rmse(etas_true,etas_initial+0.4);

% css 
figure('Position', [100 100 900 700])
hold on
plot(cssn_true,'LineWidth',2.5)
plot(cssn_initial,'LineWidth',2.5)
hold off
xlabel('Time (s)')
ylabel('c_{ss}^- (mol/m^3)','Interpreter','tex')
legend('DFN','SPMeT')
title('Anodic & Cathodic Surface Concentration')
set(gca,'FontSize',fs)

savefig('cssn_fit.fig')
print('cssn_fit','-dpng')

figure('Position', [100 100 900 700])
hold on
plot(cssp_true,'LineWidth',2.5)
plot(cssp_initial,'LineWidth',2.5)
hold off
xlabel('Time (s)')
ylabel('c_{ss}^+ (mol/m^3)','Interpreter','tex')
legend('DFN','SPMeT')
title('Anodic & Cathodic Surface Concentration')
set(gca,'FontSize',fs)

savefig('cssp_fit.fig')
print('cssp_fit','-dpng')

cssn_rmse = rmse(cssn_true,cssn_initial);
cssp_rmse = rmse(cssp_true,cssp_initial);

% ce0
figure('Position', [100 100 900 700])
hold on
plot(ce0n_true,'LineWidth',2.5)
plot(ce0n_initial,'LineWidth',2.5)
hold off
legend('DFN','SPMeT')
xlabel('Time (s)')
ylabel('c_e(0^-) (mol/m^3)','Fontsize',fs);%,'Interpreter','latex')
title('Elyte Conc. at CC (Anode)')
set(gca,'FontSize',fs)

savefig('c_en0_fit.fig')
print('c_en0_fit','-dpng')

figure('Position', [100 100 900 700])
hold on
plot(ce0p_true,'LineWidth',2.5)
plot(ce0p_initial,'LineWidth',2.5)
hold off
legend('DFN','SPMeT')
xlabel('Time (s)')
ylabel('c_e(0^+) (mol/m^3)','Fontsize',fs);%,'Interpreter','latex')
title('Elyte Conc. at CC (Cathode)')
set(gca,'FontSize',fs)

savefig('c_ep0_fit.fig')
print('c_ep0_fit','-dpng')

ce0n_rmse = rmse(ce0n_true,ce0n_initial);
ce0p_rmse = rmse(ce0p_true,ce0p_initial);

%% Save data
spmet_vs_dfn.V_true = V_true;
spmet_vs_dfn.V_sim_initial = V_sim_initial;
spmet_vs_dfn.etas_true = etas_true;
spmet_vs_dfn.etas_initial = etas_initial;
spmet_vs_dfn.cssn_true = cssn_true;
spmet_vs_dfn.cssn_initial = cssn_initial;
spmet_vs_dfn.cssp_true = cssp_true;
spmet_vs_dfn.cssp_initial = cssp_initial;

spmet_vs_dfn.v_rmse = v_rmse;
spmet_vs_dfn.etas_rmse = etas_rmse;
spmet_vs_dfn.etas_rmse_fixed = etas_rmse_fixed;
spmet_vs_dfn.cssn_rmse = cssn_rmse;
spmet_vs_dfn.cssp_rmse = cssp_rmse;
spmet_vs_dfn.ce0n_rmse = ce0n_rmse;
spmet_vs_dfn.ce0p_rmse = ce0p_rmse;

save('comparison_data.mat','spmet_vs_dfn')
