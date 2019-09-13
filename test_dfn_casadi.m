%% Simulate DFN 
% By: Zach Gima
% Date: 2019-8-2

clearvars
close all
clc

%% Load Parameters & Inputs
run param/params_NCA % loads p struct
input_folder = 'SOC60_noV/';

ID_p = struct(); 
ID_p.num_events = 1;
ID_p.event_budget = 3; % batch budget
ID_p.num_batches = 1; % num batches; batches are made of events
ID_p.collinearity_thresh = 0.7;

% Load & Parse Data for DFN syntax
data = load_data(p,input_folder);

% % Time_exp = cell(ID_p.num_events,1);
% % Current_exp = cell(ID_p.num_events,1);
% % V0 = cell(ID_p.num_events,1);
% % T_amb = cell(ID_p.num_events,1);
% % for mm = 1:ID_p.num_events
% %     Time_exp{mm} = data(mm).time;
% %     Current_exp{mm} = -data(mm).cur;
% %     V0{mm} = data(mm).V0;
% %     T_amb{mm} = data(mm).T_amb;
% % end
% 
% % Load nominal parameters
run param/params_nominal
theta_0 = Nominal_param;
theta_str = {'R_s_p';'ElecFactorDA';'epsilon_e_n';'t_plus';'R_f_n';'R_f_p'}; % exact names of the p struct fields, used for updating parameter values
% theta_str = {'t_plus'};
% 
% %% Setup selection vector and params to be identified; 
% NOTE: INDICES must line up exactly w/ theta_str
sel_k = [4;15;10;14;18;19];
% sel_k = [4;10;14;15;18;19];
% selection_vector = zeros(25,1); %G1
% selection_vector(sel_k) = 1;
% SensFlag = 0; %Turn sensitivity calculation off or on
% 
% selected_params = theta_0(sel_k);
% perturb_factor_batch = 0.95; % each batch move parameters 5%

%% For new DFN implementation
%%% Set Model Discretization Parameters
p = set_discretization(p);

%% Call DFN function to simulate "truth" voltage
for batch_idx = 1:ID_p.num_batches
  
    % simulate model
    theta_guess_current = theta_0(sel_k);
    for event_idx = 1:ID_p.num_events
%        [V_true{event_idx}, States_true{event_idx}] = DFN_sim_casadi(p,Current_exp{event_idx}, Time_exp{event_idx}, V0{event_idx}, T_amb{event_idx}, selection_vector, selected_params,SensFlag,p.R_c);
        [v_sim_cell{event_idx},states_cell{event_idx},sens_cell{event_idx}] = dfn_casadi(p, data(event_idx),theta_guess_current,theta_str);
    end
    
    % store data
%     for jj = 1:ID_p.num_events
%         data(jj).V_exp = V_true{jj};
%         data(jj).states_true = States_true{jj};
%     end
%     
%     save(strcat('DFN_truth_data_',num2str(perturb_factor_batch),'_batch_',num2str(batch_idx),'.mat'),'data')
% 
%     % Update truth params
%     selected_params = selected_params*perturb_factor_batch;
    
end

