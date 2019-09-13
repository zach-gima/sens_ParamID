%% param_remove
% Function for removing parameters from initial ID set according to
% collinearity and sensitivity analysis
% Inputs:
% Outputs: 

function [paramID_idx] = param_remove(data,ID_p,theta_char)
    %% Parse Inputs
%     STS_norm = data.STS_norm; % STS_batch is already normalized
    corr_coeff_matrix = data.corr_coeff_matrix;
    STS_norm_diag = data.STS_norm_diag;
    np = ID_p.np;
    NT = length(data.sens); % NT should correspond to the batch size
        
    %% Run collinearity analysis + elimination
    [coll_removed_idx] = findclusterparams(STS_norm_diag,corr_coeff_matrix,ID_p);
    remain_p = remaining(np, coll_removed_idx);

    disp('Removed the following parameters during collinearity clustering:');
    for pp = 1:length(coll_removed_idx)
        fprintf('%s \n', theta_char{coll_removed_idx(pp)});
    end
    fprintf('\n');
    
    %% Eliminate based on sensitivity threshold
    mu = 1E-4;
    sens_thresh = mu*sqrt(NT);
    STS_norm_diag_no_collinear = STS_norm_diag(remaining(np,coll_removed_idx));
    
    % find indices in Sens_mag_new that are below sens threshold
    sens_remove_idx = find(STS_norm_diag_no_collinear < sens_thresh); 
    
    disp('Removed the following parameters according to noise sensitivity threshold:');
    for zz = 1:length(sens_remove_idx)
        fprintf('%s \n', theta_char{remain_p(sens_remove_idx(zz))});
    end
    fprintf('\n');

    %% Determine remaining parameters
    tot_removed = [coll_removed_idx, remain_p(sens_remove_idx)];
    
    np_final = np - length(tot_removed);
    paramID_idx = remaining(np,tot_removed)';

    disp('Remaining parameters:')
    fprintf('%s \n', theta_char{paramID_idx});
    fprintf('\n');

end