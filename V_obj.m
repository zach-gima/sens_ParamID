% SPMeT Voltage Obj. Fcn for parameter identification using fmincon

% By: ZTG 2019-4-23
% Purpose: Take in current parameter estimate and map to an RMSE value
% Inputs: (1) Initial parameter estimate
%         (2) Input profile

function [ cost,fgrad ] = V_obj(theta_0, data, theta_sx, p,ID_p)
    num_events = ID_p.event_budget;
    
    %% Simulate SPMeT
    v_sim_cell = cell(num_events,1);
    sens_cell = cell(num_events,1);
    
    for ii = 1:num_events
        [v_sim_cell{ii},sens_cell{ii},~] = spmet_casadi(data(ii),theta_0,theta_sx,p,1);
    end
    
    cost = 0;
    fgrad = 0;
    
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
    
    fprintf('Total voltage RMSE = %1.6f \n \n',cost);
end