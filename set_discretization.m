%% Function to set discretization of echem model
% By: Zach Gima, 2019-5-1

function [p] = set_discretization(p, truth_model)
    if strcmp('SPMeT',truth_model)
        %%% SPMeT Specific Finite Difference Scheme
        p.Nr = 3; % 100 Make this very large so it closely approximates the true model
        p.delta_r_n = 1/p.Nr;
        p.delta_r_p = 1/p.Nr;
        p.r_vec = (0:p.delta_r_n:1)';
        %     r_vecx = r_vec(2:end-1);

        % Finite difference points along x-coordinate
        p.Nxn = 5; %70; 
        p.Nxs = 3; %35; 
        p.Nxp = 5; %70;
        p.Nx = p.Nxn+p.Nxs+p.Nxp;
    %     Nx = p.Nx - 3;
    %     x_vec_spme = linspace(0,1,Nx+4);

        p.delta_x_n = 1 / p.Nxn;
        p.delta_x_s = 1 / p.Nxs;
        p.delta_x_p = 1 / p.Nxp;

%         % Output Discretization params
%         disp('Discretization Params:');
%         fprintf(1,'No. of FDM nodes in Anode | Separator | Cathode : %1.0f | %1.0f | %1.0f\n',p.Nxn,p.Nxs,p.Nxp);
%         fprintf(1,'No. of FDM nodes in Single Particles : %1.0f\n',p.Nr);
%     %     fprintf(1,'Time Step : %2.2f sec\n',p.delta_t);
%         disp(' ');
    
    elseif strcmp('DFN',truth_model)
        %% Discretization parameters
        %%% DFN Specific Finite Difference Scheme
        % Pade Order
        p.PadeOrder = 3;

        % Finite difference points along x-coordinate
        p.Nxn = 70;
        p.Nxs = 35;
        p.Nxp = 70;
        p.Nx = p.Nxn+p.Nxs+p.Nxp;

        p.delta_x_n = 1 / p.Nxn;
        p.delta_x_s = 1 / p.Nxs;
        p.delta_x_p = 1 / p.Nxp;
    else
        error('Please specify either "SPMeT" or "DFN"')
    end
end