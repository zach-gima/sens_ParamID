%% Function to update parameters that are functions of other parameter values
% By: Zach Gima 2019-5-17
function p = update_dependencies(p)

    % make element to caclulate phi_{s} by Saehong Park 
    p.epsilon_f_n = 1 - p.epsilon_s_n - p.epsilon_e_n;  % Volume fraction of filler in neg. electrode
    p.epsilon_f_p = 1 - p.epsilon_s_p - p.epsilon_e_p;  % Volume fraction of filler in pos. electrode

    % Specific interfacial surface area
    p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
    p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]
    
    % Electrolyte concentration matrices: These are parameter dependent so need to re-generate everytime we
    % change the parameters
    p = update_c_e_mats(p);
end