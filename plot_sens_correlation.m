%% Need to Normalize the Sensitivity

clear all
close all
clc

%% File I/O and parameters
run params/params_NCA
run params/params_bounds
addpath('param_remove/') % functions needed for parameter elimination

num_params = 6;

sens_results_folder = 'sens-outputs/';
files = dir('/Users/ztakeo/Documents/GitHub/sens_ParamID/sens-outputs');
files=files(~ismember({files.name},{'.','..','.DS_Store'}));
% filenames = files.name;

%%%% load below for 6 drive cycles
% load('/Users/ztakeo/Desktop/ACC2020/sens_data.mat','sens_initial_spmet','sens_initial_dfn')
% num_exp = 6;

num_exp = length(files);

%% iterate over SPMeT sensitivity
% sens_norm_spmet = cell(num_exp,1);
sens_norm_spmet_diag = zeros(num_exp,num_params); 
for ii = 1:num_exp
    load(strcat(sens_results_folder,files(ii).name),'sens_initial_spmet');
    [~,sens_norm_spmet_diag(ii,:),~] = compute_sens_variables(p,p_bounds,sens_initial_spmet);
%     sens_current_exp = sens_initial_dfn{ii};
%     [~,sens_norm_spmet_diag(ii,:),~] = compute_sens_variables(p,p_bounds,sens_current_exp); 
end
x = median(sens_norm_spmet_diag);
xneg = x - min(sens_norm_spmet_diag);
xpos = max(sens_norm_spmet_diag) - x;

% each column of sens_norm_spmet_diag is the norm of the corresponding
% parameter for each experiment (rows)

%% iterate over DFN sensitivity
% sens_norm_dfn = cell(num_exp,1);
sens_norm_dfn_diag = zeros(num_exp,num_params); 
for ii = 1:num_exp
    load(strcat(sens_results_folder,files(ii).name),'sens_initial_dfn');
    [~,sens_norm_dfn_diag(ii,:),~] = compute_sens_variables(p,p_bounds,sens_initial_dfn);
%     sens_current_exp = sens_initial_dfn{ii};
%     [~,sens_norm_dfn_diag(ii,:),~] = compute_sens_variables(p,p_bounds,sens_current_exp);
end
y = median(sens_norm_dfn_diag);
yneg = y - min(sens_norm_dfn_diag);
ypos = max(sens_norm_dfn_diag) - y;


%% plot for parameters grouped
% fs = 25;
% figure('Position', [100 100 900 700])
% e = errorbar(x,y,yneg,ypos,xneg,xpos,'o');
% e.MarkerSize = 9;
% e.MarkerFaceColor = 'r';
% e.MarkerEdgeColor = 'r';
% e.LineWidth = 1.5;
% xlabel('SPMeT Sensitivity, $\frac{\partial V}{\partial \theta}$','Interpreter','latex')
% ylabel('DFN Sensitivity, $\frac{\partial V}{\partial \theta}$','Interpreter','latex')
% set(gca, 'XScale','log', 'YScale','log','FontSize',fs)
% % box on
% grid on

%% plot for no parameter grouping
fs = 25;
figure('Position', [100 100 900 700])
% e = errorbar(x,y,yneg,ypos,xneg,xpos,'o');
% e.MarkerSize = 9;
% e.MarkerFaceColor = 'r';
% e.MarkerEdgeColor = 'r';
% e.LineWidth = 1.5;

% plot sensitivities
sens_spmet_col_vec = [];
sens_dfn_col_vec = [];
hold on
for ii = 1:num_params
    scatter(sens_norm_spmet_diag(:,ii),sens_norm_dfn_diag(:,ii),250,'LineWidth',3)
    sens_spmet_col_vec = vertcat(sens_spmet_col_vec,sens_norm_spmet_diag(:,ii));
    sens_dfn_col_vec = vertcat(sens_dfn_col_vec,sens_norm_dfn_diag(:,ii));
end

% plot line of best fit
[linearCoefficients,S] = polyfit(log(sens_spmet_col_vec),log(sens_dfn_col_vec),1);

z = polyval(linearCoefficients,log(sens_spmet_col_vec));
plot(sens_spmet_col_vec,exp(z),'k','LineWidth',1.5)

hold off

% Display R^2 value
R = corrcoef(sens_spmet_col_vec,sens_dfn_col_vec); % 0.9966
R = round(R(1,2),3);
% R_2 = R^2; %0.993
% txt = text(1,10^-3,['R = ',num2str(R)]);
% txt(1).FontSize = fs;
% txt(1).FontWeight = 'bold';

% title('No Parameter Grouping')
xlabel('SPMeT Sensitivity, $\partial V / \partial \theta$','Interpreter','latex')
ylabel('DFN Sensitivity, $\partial V / \partial \theta$','Interpreter','latex')
set(gca, 'XScale','log', 'YScale','log','FontSize',fs)
% box on
% grid on

function exp_idx = find_outlier_exp(sens_norm_spmet_diag,sens_norm_dfn_diag)
    x = log(sens_norm_spmet_diag(:,1));
    y = log(sens_norm_dfn_diag(:,1));

    log_err_vec = y - linearCoefficients(1)*x - linearCoefficients(2);
    log_abs_err_vec_sorted = sort(abs(log_err_vec),'descend');

    exp_idx = find(abs(log_err_vec) > 2);
   
end
