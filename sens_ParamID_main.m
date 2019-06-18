%% Sensitivity-based Online Parameter Estimation for SPMeT
%%% FUNCTIONALIZE THIS
% By: Zach Gima, 2019-4-12
clear;
clc;
close all;

% FontSize
fs = 25;

datetime_initial = datetime('now','TimeZone','America/Los_Angeles');
date_txt = strrep(datestr(datetime_initial), ':', '_');

addpath('/Users/ztakeo/Documents/MATLAB/casadi') % Mac Laptop
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
baseline = 'C'; %'C'; %'B' 'C'
data_select_logic = 1; % 1 for selecting data according to sensitivity content
soc_0 = 'SOC60'; %SOC#
input_folder = strcat('input-data/SPMeT/',soc_0,'/');
% input_filename = 'US06x3_batt_ObsData.mat';
% input_path = strcat(input_folder,input_filename);
data = load_data(p,input_folder);

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
ID_p.num_batches = 1; % num batches; batches are made of events
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
perturb_factor = 1.3;
theta.truth = [p.R_s_p;p.ElecFactorDA;p.epsilon_e_n;p.t_plus;p.R_f_n;p.R_f_p]; % initial value
theta.guess = perturb_factor*theta.truth;

ID_p.np = length(theta.guess);
params_all_idx = (1:ID_p.np)'; % vector of indices for all parameters, used for initial model simulation

theta.history = zeros(ID_p.np,ID_p.num_batches+1); % array used to track identified values for each parameter over all of the batches
theta.history(:,1) = theta.guess;
theta.char = {'R_s^+';'ElecFactorDA';horzcat(char(949),'_e^-');'t_+';'R_f^-';'R_f^+'}; % char version of the params, used for plotting/fprintf purposes
theta.str = {'R_s_p';'ElecFactorDA';'epsilon_e_n';'t_plus';'R_f_n';'R_f_p'}; % exact names of the p struct fields, used for updating parameter values

% update p struct with initial guess
p = update_p_struct(p,theta.guess,theta.str);

%% Run ParamID routine
% Initialize variables
sens_data(ID_p.num_batches) = struct();
paramID_idx = cell(ID_p.num_batches,1);
datetime_final = cell(ID_p.num_batches,1);

fprintf('\n')
disp('*********BEGINNING PARAMETER ID ROUTINE*********')
fprintf('Baseline %s \n',baseline)
fprintf('Start Time: %s \n \n',datetime_initial);
tic

for batch_idx = 1:ID_p.num_batches
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('Batch #%i \n',batch_idx)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    %% Simulate SPMeT for initial parameter guess: voltage, sensitivity
    % Used for event_select
    V_sim_initial = cell(ID_p.num_events,1);
    sens_initial = cell(ID_p.num_events,1);
    state_info_initial = cell(ID_p.num_events,1);
    
    % Simulate model for Voltage, Sensitivity
    % parfor loops don't behave well with structures; assign temp variable
    theta_guess_initial = theta.guess; 
    theta_str_initial = theta.str;
    
    parfor ii = 1:ID_p.num_events
        [V_sim_initial{ii},state_info_initial{ii},sens_initial{ii}] = spmet_casadi(p,data(ii),theta_guess_initial,theta_str_initial);
        [V_check{ii},state_check{ii}] = spmet_casadi(p,data(ii));

    end
%     clear theta_guess_initial theta_str_initial
    
    % Compute STS_norm for selecting optimal data
    STS_norm_initial = cell(ID_p.num_events,1);
    for ii = 1:ID_p.num_events
        data(ii).V_sim_initial = V_sim_initial{ii};
        data(ii).sens_initial = sens_initial{ii};

        %%% ZTG Note 2019-7-11: compute_sens_variables recomputes the
        %%% normalize_sens_factor -- should this factor update with
        %%% parameter updates?
        [STS_norm_initial{ii},~,~] = compute_sens_variables(p,bounds,sens_initial{ii});
    end

    % Debugging -- save data so don't have to regenerate every time
    save(strcat(output_folder,date_txt,'_initial_sim'),'data','STS_norm_initial','state_info_initial');
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
    opt_data = data(opt_event_idx); % create new struct of just the optimal data

    %% Compute additional sensitivity data needed for parameter collinearty + sensitivity analysis
    % Compute other sensitivity information for the concatenated sens data

    sens_data(batch_idx).sens = vertcat(opt_data.sens_initial); 

    [STS_norm,STS_norm_diag,corr_coeff_matrix] = compute_sens_variables(p,bounds,sens_data(batch_idx).sens);
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
        [paramID_idx{batch_idx}] = param_remove(sens_data(batch_idx),ID_p,theta.char);
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
    for ii = 1:ID_p.event_budget
        cost_initial = cost_initial + sqrt(mean((opt_data(ii).V_exp - opt_data(ii).V_sim_initial).^2));
    end
    fprintf('Initial Voltage RMSE = %1.6f \n \n',cost_initial);
    
    % Set parameters to be identified based on the batch
    current_params_idx = paramID_idx{batch_idx};
    lb = bounds.min(current_params_idx);
    ub = bounds.max(current_params_idx);    
    theta_guess_current = theta.guess(current_params_idx);
    theta_str_current = theta.str(current_params_idx);
    
    % create anonymous function to pass in other parameters needed for V_obj
    fh = @(x)V_obj(x,opt_data,theta_str_current,p);
    
    [theta_ID,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN]= ...
       fmincon(fh,theta_guess_current,[],[],[],[],lb,ub,[],opt);

    % Update p struct with newly ID'ed values
    p = update_p_struct(p,theta_ID,theta_str_current);
    
    % Store identified values & update guess for next event
    theta.guess(current_params_idx) = theta_ID;
    theta.history(:,batch_idx+1) = theta.history(:,batch_idx); % for parameters not identified this batch, this stores the previous value (i.e. it stays constant)
    theta.history(current_params_idx,batch_idx+1) = theta_ID;
    
    % Update computation time-related metrics
    wallclock = toc;
    datetime_final{batch_idx} = datetime('now','TimeZone','America/Los_Angeles');
    fprintf('Finished in %i seconds at %s \n',wallclock,datetime_final{batch_idx})
    
    % re-simulate model for this batch and final parameter set identified
    V_sim = cell(ID_p.num_events,1);
    state_info = cell(ID_p.num_events,1);
    sens = cell(ID_p.num_events,1);
    
    parfor ii = 1:ID_p.num_events
        %%%% ZTG Note 2019-6-17: should I be passing in just the parameters
        %%%% we sought to identify? or all of the params?
        % difference between passing theta_guess_current, theta_str_current
        % vs passing in all the params i.e.:
        %     theta_guess_initial = theta.guess; 
        %     theta_str_initial = theta.str;
        [V_sim{ii},state_info{ii},sens{ii}] = spmet_casadi(p,data(ii),theta_guess_initial,theta_str_initial);
%         [V_check{ii},state_check{ii}] = spmet_casadi(p,data(ii));
    end
    
    % save current batch data
    ID_out.V_sim = V_sim;
    ID_out.state_info = state_info;
    ID_out.sens = sens;
    ID_out.data = data;
    ID_out.params_final_idx = paramID_idx;
    ID_out.opt_event_idx = opt_event_idx;
    ID_out.datetime_initial = datetime_initial;
    ID_out.datetime_final = datetime_final;
    ID_out.wallclock = wallclock;
    ID_out.theta = theta;
    ID_out.perturb_factor = perturb_factor;
    ID_out.fmincon_out = OUTPUT;
    ID_out.fmincon_fval = FVAL;
    ID_out.fmincon_exitflag = EXITFLAG;
    ID_out.fmincon_history = history;
    ID_out.fmincon_searchdir = searchdir;
    
    save(strcat(output_path,'_batch_',num2str(batch_idx)),'ID_out')
    clear ID_out

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
            plot(x(1),x(2),'o');
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