%% Compute Sensitivity-related variables for each event in a batch of data
% By: Zach Gima, 2019-6-6
function [STS_norm,STS_norm_diag,corr_coeff_matrix] = compute_sens_variables(p,bounds,sens)
    
    % Used for evaluating sensitivity content of 
    normalize_sens_factor = normalizesens(p,bounds);
    sens_norm = normalize_sens_factor.*sens; % normalized sensitivity 
    STS_norm = sens_norm'*sens_norm;
    
    % Used for collinearity analysis
    STS_norm_diag = diag(STS_norm);
    ttt = inv(diag(sqrt(diag(STS_norm)))); % matrix of the diagonal elements of STS_norm -- these elements are the normalized (w.r.t. each other parameter) sensitivity magnitude
    corr_coeff_matrix = abs(ttt*STS_norm*ttt);  % Correlation coefficient matrix
end
