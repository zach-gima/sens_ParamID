%% Function to overwrite parameters to be identified in the p struct with the CasADi SX type
% This is necessary for sensitivity computations

% By: Zach Gima 2019-6-10
function [p,theta_sx] = update_casadi_vars(p,theta_str)
    import casadi.*
    np = length(theta_str);
    theta_sx = SX.zeros(np,1);
    
    % iterate through all of the parameters and generate casadi variables
    for mm = 1:np
        theta_sx(mm) = SX.sym(theta_str{mm});
        p.(theta_str{mm}) = theta_sx(mm);
    end
end