% SPMeT Voltage Obj. Fcn for parameter identification using fmincon

% By: ZTG 2019-4-23
% Purpose: Take in current parameter estimate and map to an RMSE value
% Inputs: (1) Initial parameter estimate
%         (2) Input profile

function [ cost,fgrad ] = V_obj(theta_0, data, theta_str, p)
    global function_evals
    
    num_events = length(data);
    
    %% Simulate SPMeT
    v_sim_cell = cell(num_events,1);
    sens_cell = cell(num_events,1);
    
    parfor ii = 1:num_events
        [v_sim_cell{ii},~,sens_cell{ii}] = spmet_casadi(p,data(ii),theta_0,theta_str);        
    end
    
    cost = 0;
    fgrad = 0;
    
    % Cost function is the sum of non-contiguous events -- add up cost and gradient for
    % each event!
    for ii = 1:num_events
        
        % Parse current event data
        v_dat = data(ii).V_exp; % truth/measured voltage data we're fitting to   
        v_sim = v_sim_cell{ii};
        sens = sens_cell{ii};
        N = length(v_dat);
        
        % Compute MSE & RMSE along with their gradients
        MSE = mean((v_dat - v_sim).^2);    
        RMSE = sqrt(MSE);

        fgrad_MSE = -2/N*sens'*(v_dat - v_sim); % MSE gradient
        fgrad_RMSE = 1/(2*sqrt(MSE))*fgrad_MSE; %RMSE gradient

        % Sum these quantities for non-contiguous event data
        fgrad = fgrad_RMSE + fgrad;
        cost = RMSE + cost; %MSE
    end
    
    function_evals = function_evals + 1;
    
    fprintf('Total voltage RMSE = %1.6f \n \n',cost);
end