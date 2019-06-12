%% Initial Conditions for SPMeT
% By: Zach Gima, 2019-5-1
function [p,x0] = init_spmet_nocasadi(p,V0,T_amb,print_flag)
 
    % Solid concentration
    [csn0,csp0] = init_cs(p,V0);
    c_n0 = csn0 * ones(p.Nr-1,1);
    c_p0 = csp0 * ones(p.Nr-1,1);

    % Electrolyte concentration
    Nx = p.Nx - 3;
    ce0 = p.c_e*ones(Nx,1);

    % Temperature -- MAKE SURE IN KELVIN
%     T10 = p.T_amb;
%     T20 = p.T_amb;
    T10 = T_amb + 273.15;
    T20 = T_amb + 273.15;

    % SEI layer
    delta_sei0 = 0;

    x0 = [c_n0; c_p0; ce0; T10; T20; delta_sei0];
    
    %% Generate Constant System Matrices
    % Electrolyte concentration matrices
    % These are parameter dependent so need to re-generate everytime we
    % update the parameters
    p = update_c_e_mats_nocasadi(p);
    
    if print_flag
        % Output Initial Conditions
        disp('Initial Conditions:');
        fprintf(1,'Voltage : %1.3f V\n',V0);
        fprintf(1,'Normalized Solid Concentration in Anode | Cathode : %1.2f | %1.2f\n',csn0/p.c_s_n_max,csp0/p.c_s_p_max);
        fprintf(1,'Electrolyte Concentration : %2.3f kmol/m^3\n',ce0(1)/1e3);
        fprintf(1,'Temperature in Roll | Can : %3.2f K | %3.2f K \n',T10,T20);
        fprintf(1,'SEI Layer in Anode : %2f um \n',delta_sei0*1e6);
        disp(' ');
    end
end