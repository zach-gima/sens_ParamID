%% Function to overwrite parameters to be identified in the p struct with the CasADi SX type
% This is necessary for sensitivity computations

% By: Zach Gima 2019-6-10
function [p] = update_casadi_vars(p,theta,ID_p,current_params_idx,jj)
    % iterate through all of the parameters intially specified to potentially be ID'ed
    for mm = 1:ID_p.np 
        if ismember(mm,current_params_idx) % parameter being identified
            % then overwrite with SX type
            p.(theta.str{mm}) = theta.sx(mm);
        else % (parameter not being identified)
            % want to use most up-to-date estimated parameter value -- could
            % have re-identified a batch or two ago, so want to use that value
            p.(theta.str{mm}) = theta.history(mm,jj);
        end
    end
    
    % Update parameter dependencies: some parameters are functions of the above
    % specified parameters
    p = update_dependencies(p);
end