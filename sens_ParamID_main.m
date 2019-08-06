%% Sensitivity-based Online Parameter Estimation for SPMeT
% By: Zach Gima, 2019-4-12
clear;
clc;
close all;

% FontSize
fs = 25;

datetime_initial = datetime('now','TimeZone','America/Los_Angeles');
date_txt = strrep(datestr(datetime_initial), ':', '_');

% addpath('/Users/ztakeo/Documents/MATLAB/casadi') % Mac Laptop
% addpath('C:/Users/Zach/Documents/MATLAB/casadi_windows') % HPC-1
% addpath('C:/Users/zgima/Documents/MATLAB/casadi_windows') % HPC-2
% addpath('/global/home/users/ztakeo/modules/casadi-matlab');    % For Savio

%% instantiate global variables that will track iterations within fmincon & optimization options object
global history;
global searchdir;
history.x = [];
history.fval = [];
searchdir = [];

% opt sets the algorithm, "output function" this was an addition for debugging that
%helps us look inside fmincon
opt = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
'SpecifyObjectiveGradient',true,'Diagnostics','on','OutputFcn',@outfun);  %'CheckGradients',true,'FiniteDifferenceType','central',

%% User Inputs

%%% Load Parameters
run param/params_NCA
run param/params_bounds
% run param/params_LCO % loads p struct

%%% Set input/output paths & load data
baseline = 'A'; %'C'; %'B' 'C'
data_select_logic = 0; % 1 for selecting data according to sensitivity content
truth_model = 'DFN';
soc_0 = 'SOC60'; %SOC#
input_folder = strcat('input-data/',truth_model,'/Training Data/',soc_0,'/');
% input_filename = 'US06x3_batt_ObsData.mat';
% input_path = strcat(input_folder,input_filename);
% data = load_data(p,input_folder);
truth_filename = strcat(input_folder,'DFN_truth_data_1.05_batch_1.mat');
load(truth_filename,'data')

output_folder = strcat('output-data/','Baseline',baseline,'/',date_txt,'/');
mkdir(output_folder); %create new subfolder with current date in output_folder
output_path = strcat(output_folder,date_txt,'_results');

% For saving errors:
error_filename = strcat(output_folder,date_txt,'_sim_log.txt');
diary(error_filename)

%%% setup structure to hold all of param-ID related hyperparameters
ID_p = struct(); 
ID_p.num_events = length(data);
ID_p.event_budget = 3; % batch budget
ID_p.num_batches = 5; % num batches; batches are made of events
ID_p.collinearity_thresh = 0.7;

%%% Set Model Discretization Parameters
p = set_discretization(p);

%% Specify the parameters to identify:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_s_p -- "high" sens. param
% p.ElecFactorDA -- "medium" sens param
% t_plus -- "low" sens. param
% p.epsilon_e_n -- idk...
% R_f_n and R_f_p -- one should be eliminated via collinearity because
% they're directly related
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify parameters to be identified -- make sure each theta_ variable specifies the parameters
% in the same order
perturb_factor_initial = 1; %1.3;
perturb_factor_batch = 1.05; % each batch move parameters 5%
theta(ID_p.num_batches) = struct();
theta(1).truth = [p.R_s_p;p.ElecFactorDA;p.epsilon_e_n;p.t_plus;p.R_f_n;p.R_f_p]; % initial value
theta(1).guess = perturb_factor_initial*theta(1).truth;
theta(1).final_ID = theta(1).guess; % initialize the final ID as the guess value; parameters selected to be ID'ed will overwrite their corresponding index later

ID_p.np = length(theta(1).guess);
params_all_idx = (1:ID_p.np)'; % vector of indices for all parameters, used for initial model simulation


%% Run ParamID routine
% Initialize variables
sens_data(ID_p.num_batches) = struct();
paramID_idx = cell(ID_p.num_batches,1);
datetime_final = cell(ID_p.num_batches,1);

fprintf('\n')
disp('*********BEGINNING PARAMETER ID ROUTINE*********')
fprintf('Baseline %s \n',baseline)
fprintf('Data selection: %i \n', data_select_logic)
fprintf('Start Time: %s \n \n',datetime_initial);
tic

for batch_idx = 1:ID_p.num_batches
    fprintf('\n')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('Batch #%i \n',batch_idx)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('\n')

    %% Setup parameter structure for this batch
    theta(batch_idx).char = {'R_s^+';'ElecFactorDA';horzcat(char(949),'_e^-');'t_+';'R_f^-';'R_f^+'}; % char version of the params, used for plotting/fprintf purposes
    theta(batch_idx).str = {'R_s_p';'ElecFactorDA';'epsilon_e_n';'t_plus';'R_f_n';'R_f_p'}; % exact names of the p struct fields, used for updating parameter values
    
    % update p struct with initial guess
    p = update_p_struct(p,theta(batch_idx).guess,theta(batch_idx).str);
    
    %% Simulate SPMeT for truth parameters
%     % need to do this if changing truth parameters each batch
%     theta_truth_current = theta(batch_idx).truth;
%     p_truth = update_p_struct(p,theta_truth_current,theta(batch_idx).str);
% 
%     disp('%%%%%%%%%%%%%')
%     disp('Simulating truth data')
%     disp('%%%%%%%%%%%%%')
%     parfor mm = 1:ID_p.num_events
%         [V_true{mm},States_true{mm}] = spmet_casadi(p_truth,data(mm))
%     end
%     
%     % store data
%     for jj = 1:ID_p.num_events
%         data(jj).V_exp = V_true{jj};
%         data(jj).states_true = States_true{jj};
%     end
    
    %% If using DFN as truth model
    truth_filename = strcat(input_folder,'DFN_truth_data_1.05_batch_',num2str(batch_idx),'.mat');
    load(truth_filename,'data')

    %% Simulate SPMeT for initial parameter guess: voltage, sensitivity
    % Used for event_select
    V_sim_initial = cell(ID_p.num_events,1);
    sens_initial = cell(ID_p.num_events,1);
    states_initial = cell(ID_p.num_events,1);
    
    % Simulate model for Voltage, Sensitivity
    % parfor loops don't behave well with structures; assign temp variable
    theta_guess_initial = theta(batch_idx).guess; 
    theta_str_initial = theta(batch_idx).str; % initially pass in all variables because we haven't decided yet which params to ID
    
    fprintf('\n')
    disp('%%%%%%%%%%%%%')
    disp('Simulating for initial parameter guess')
    disp('%%%%%%%%%%%%%')
    fprintf('\n')
    parfor ii = 1:ID_p.num_events
        [V_sim_initial{ii},states_initial{ii},sens_initial{ii}] = spmet_casadi(p,data(ii),theta_guess_initial,theta_str_initial);
    end
    clear theta_guess_initial theta_str_initial
    
    % Compute STS_norm for selecting optimal data
    STS_norm_initial = cell(ID_p.num_events,1);
    for ii = 1:ID_p.num_events
        data(ii).V_sim_initial = V_sim_initial{ii};
        data(ii).sens_initial = sens_initial{ii};
        %%% ZTG Note 2019-7-11: compute_sens_variables recomputes the
        %%% normalize_sens_factor -- should this factor update with
        %%% parameter updates?
        [STS_norm_initial{ii},~,~] = compute_sens_variables(p,p_bounds,sens_initial{ii});
    end

    % Debugging -- save data so don't have to regenerate every time
%     save(strcat(output_folder,date_txt,'_initial_sim'),'data','STS_norm_initial','states_initial');
%     load('initial_sim.mat');

    %% Select events
    % opt_event_idx contains the indices corresponding to events that were
    % selected for having the highest sensitivity content w.r.t. the parameters
    % of interest
    if data_select_logic == 1
        opt_event_idx = event_select(ID_p,STS_norm_initial);
    else 
        opt_event_idx = (1:ID_p.num_events)'; % all events selected
    end
    ID_p.num_opt_events = length(opt_event_idx);
    opt_data = data(opt_event_idx); % create new struct of just the optimal data

    %% Compute additional sensitivity data needed for parameter collinearty + sensitivity analysis
    % Compute other sensitivity information for the concatenated sens data
    sens_data(batch_idx).sens = vertcat(opt_data.sens_initial); 

    [STS_norm,STS_norm_diag,corr_coeff_matrix] = compute_sens_variables(p,p_bounds,sens_data(batch_idx).sens);
    sens_data(batch_idx).STS_norm = STS_norm;
    sens_data(batch_idx).STS_norm_diag = STS_norm_diag;
    sens_data(batch_idx).corr_coeff_matrix = corr_coeff_matrix;
    
   %% Eliminate parameters from identification routine according to collinearity + sensitivity analysis
    % For each batch of optimal data, determine which parameters to identify 
   % Baseline A: All Params
    if strcmp('A',baseline) 
        paramID_idx{batch_idx} = (1:ID_p.np)';
    % Baseline B and C
    elseif strcmp('B',baseline) || strcmp('C',baseline)
        [paramID_idx{batch_idx}] = param_remove(sens_data(batch_idx),ID_p,theta(batch_idx).char);
    else
        error('Baseline improperly specified. Set baseline = ''A'',''B'', or ''C'' ')
    end
           
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DEBUG   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure('Position', [100 100 900 700])
%     plot(opt_data(1).V_exp,'LineWidth', 2.5);
%     hold on
%     plot(opt_data(1).V_sim_initial,'LineWidth', 2.5);
%     hold off
%     xlabel('Time (s)')
%     ylabel('Voltage (V)')
%     title('Initial Fit')
%     legend('Truth Data','Simulated Data')
%     set(gca,'Fontsize',fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DEBUG   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Identify Parameters
    % Print Initial Cost Function Eval
    cost_initial = 0;  
    for ii = 1:ID_p.num_opt_events
        cost_initial = cost_initial + sqrt(mean((opt_data(ii).V_exp - opt_data(ii).V_sim_initial).^2));
    end
    fprintf('\n')
    disp('%%%%%%%%%%%%%')
    fprintf('Initial Summed Voltage RMSE = %1.6f \n',cost_initial);
    fprintf('Initial Average Voltage RMSE = %1.6f \n',cost_initial/ID_p.num_opt_events);
    disp('%%%%%%%%%%%%%')
    fprintf('\n')
    
    % Set parameters to be identified based on the batch
    current_params_idx = paramID_idx{batch_idx};
    lb = p_bounds.min(current_params_idx);
    ub = p_bounds.max(current_params_idx);    
    theta_guess_current = theta(batch_idx).guess(current_params_idx);
    theta_str_current = theta(batch_idx).str(current_params_idx);
    
    % create anonymous function to pass in other parameters needed for V_obj
    fh = @(x)V_obj(x,opt_data,theta_str_current,p);
    
    fprintf('\n')
    disp('%%%%%%%%%%%%%')
    disp('Identifying parameters')
    disp('%%%%%%%%%%%%%')
    fprintf('\n')

    [theta_ID,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN]= ...
       fmincon(fh,theta_guess_current,[],[],[],[],lb,ub,[],opt);

    % Update computation time-related metrics
    wallclock = toc;
    datetime_final{batch_idx} = datetime('now','TimeZone','America/Los_Angeles');
    fprintf('Finished in %i seconds at %s \n',wallclock,datetime_final{batch_idx})
    
    % Update p struct with newly ID'ed values
    p = update_p_struct(p,theta_ID,theta_str_current);
    
    % Store identified values
    theta(batch_idx).final_ID(current_params_idx) = theta_ID;
    
    clear theta_guess_current theta_str_current
    
    %% re-simulate model for the final parameter set identified: need updated voltage profile and sens data for Conf. Int.
    V_sim_final = cell(ID_p.num_events,1);
    states_final = cell(ID_p.num_events,1);
    sens_final = cell(ID_p.num_events,1);
    
    % create temp params to pass to spmet_casadi
    %%%%%% ZTG NOTE: CHECK THIS MAY NOT BE BUG FREE
    theta_ID_final =  theta(batch_idx).final_ID(current_params_idx);
    theta_str_final = theta(batch_idx).str(current_params_idx);
    
    fprintf('\n')
    disp('%%%%%%%%%%%%%')
    disp('Re-simulating model for final parameter values')
    disp('%%%%%%%%%%%%%')
    fprintf('\n')
    parfor ii = 1:ID_p.num_events
        [V_sim_final{ii},states_final{ii},sens_final{ii}] = spmet_casadi(p,data(ii),theta_ID_final,theta_str_final);
    end
    
    clear theta_ID_final theta_str_final
        
    %% save current batch data
    ID_out.baseline = baseline;
    ID_out.data_select_logic = data_select_logic;
    ID_out.truth_model = truth_model;
    ID_out.soc_0 = soc_0;
    ID_out.V_sim_initial = V_sim_initial;
    ID_out.V_sim_final = V_sim_final;
    ID_out.states_initial = states_initial;
    ID_out.states_final = states_final;
    ID_out.sens_final = sens_final;
    
    ID_out.data = data;
    ID_out.params_final_idx = paramID_idx;
    ID_out.opt_event_idx = opt_event_idx;
    ID_out.datetime_initial = datetime_initial;
    ID_out.datetime_final = datetime_final;
    ID_out.wallclock = wallclock;
    ID_out.theta = theta(batch_idx);
    ID_out.perturb_factor_initial = perturb_factor_initial;
    ID_out.perturb_factor_batch = perturb_factor_batch;
    ID_out.fmincon_out = OUTPUT;
    ID_out.fmincon_exitflag = EXITFLAG;
    ID_out.fmincon_history = history;
    ID_out.fmincon_searchdir = searchdir;
    
    save(strcat(output_path,'_batch_',num2str(batch_idx)),'ID_out')
    clear ID_out

    %% Reset fmincon global variables
    history.x = [];
    history.fval = [];
    searchdir = [];
    
    %% Update truth params, guess, and final_ID vector for next batch
    if batch_idx < ID_p.num_batches
        theta(batch_idx+1).guess = theta(batch_idx).final_ID;
        theta(batch_idx+1).final_ID = theta(batch_idx+1).guess;
        theta(batch_idx+1).truth = perturb_factor_batch*theta(batch_idx).truth;
    end
end

diary off

%% The purpose of this function is to save the iteration values 
%found and slightly modified from online documentation
function stop = outfun(x,optimValues,state)
    stop = false;
    global history;
    global searchdir;
    
    switch state
        case 'init'
            figure('Position', [100 100 900 700])
            hold on
        case 'iter'
            % Concatenate current point and objective function
            % value with history. x must be a row vector.
            history.fval = [history.fval; optimValues.fval];
            history.x = [history.x; x'];
            % Concatenate current search direction with 
            % searchdir.
            searchdir = [searchdir;... 
                        optimValues.searchdirection'];
%             plot(x(1),x(2),'o');
            % Label points with iteration number and add title.
            % Add .15 to x(1) to separate label from plotted 'o'
            text(x(1)+.15,x(2),... 
                num2str(optimValues.iteration));
            title('Sequence of Points Computed by fmincon');
        case 'done'
            hold off
    otherwise
    end
end