%% By: Zach Gima 2019-6-12
% Function to update parameters in the p struct with newly identified
% values
function p = update_p_struct(p,theta_ID,theta_str)

    np = length(theta_ID);
    % iterate through all of the parameters intially specified to potentially be ID'ed
    for mm = 1:np
        p.(theta_str{mm}) = theta_ID(mm);
    end

    %% Also Update Parameter Dependencies
    p = update_dependencies(p);
end