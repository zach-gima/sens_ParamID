%% Initial Conditions for DFN
% By: ZTG, 2019-9-11

function[x0,z0,SOC,c_avg_n] = init_dfn(p,NT,V0,T_amb,SensFlag)
    import casadi.*
    % Given V0, obtain 
    [csn0,csp0] = init_cs(p,V0);
%     V0 = refPotentialCathode(csp0/p.c_s_p_max)-refPotentialAnode(csn0/p.c_s_n_max);

    % Electrolyte concentration
    ce0 = p.c_e;

    % Vector lengths
    Ncsn = p.PadeOrder * (p.Nxn-1);
    Ncsp = p.PadeOrder * (p.Nxp-1);
    Nn = p.Nxn - 1;
    Ns = p.Nxs - 1;
    Np = p.Nxp - 1;
    Nx = p.Nx - 3;

    c_s_n0 = zeros(p.PadeOrder,1);
    c_s_p0 = zeros(p.PadeOrder,1);

    %%%%% Initial condition based on Jordan form
    c_s_n0(3) = csn0;
    c_s_p0(3) = csp0;

    c_s_n = zeros(Ncsn,NT);
    c_s_p = zeros(Ncsp,NT);

    c_s_n(:,1) = repmat(c_s_n0, [Nn 1]);
    c_s_p(:,1) = repmat(c_s_p0, [Nn 1]);

    % Electrolyte concentration
    if SensFlag == 1
        c_e = SX.zeros(Nx,NT);    
    else
        c_e = zeros(Nx,NT); %CASADI CHANGE
    end
    c_e(:,1) = ce0 * ones(Nx,1);

    % Temperature
    % T = zeros(NT,1);
    % T(1) = T0;
%     T10 = p.T_amb;
%     T20 = p.T_amb;
    T10 = T_amb(1) + 273.15; % T_amb expected in Celsius
    T20 = T_amb(1) + 273.15;

    p.T_amb = T_amb(1) + 273.15;
    
    % Solid Potential
    Uref_n0 = refPotentialAnode(csn0(1)*ones(Nn,1) / p.c_s_n_max);
    Uref_p0 = refPotentialCathode(csp0(1)*ones(Np,1) / p.c_s_p_max);

    phi_s_n = zeros(Nn,NT);
    phi_s_p = zeros(Np,NT);
    phi_s_n(:,1) = Uref_n0;
    phi_s_p(:,1) = Uref_p0;

    % Electrolyte Current
    i_en = zeros(Nn,NT);
    i_ep = zeros(Np,NT);

    % Electrolyte Potential
    phi_e = zeros(Nx+2,NT);

    % Molar Ionic Flux
    jn = zeros(Nn,NT);
    jp = zeros(Np,NT);

    %% Return values below
    % Volume average concentration
    c_avg_n = zeros(Nn,NT);
    c_avg_n(:,1) = repmat(csn0, [Nn 1]);
    
    % SOC (Bulk Anode SOC)
    SOC = zeros(NT,1);
    SOC(1) = (mean(c_avg_n(:,1)) - p.cn_low) / (p.cn_high - p.cn_low);%soc00;%(mean(c_avg_n(:,1)) - cn_low) / (cn_high - cn_low);

    % Initial Conditions
    x0 = [c_s_n(:,1); c_s_p(:,1); c_e(:,1); T10; T20];

    z0 = [phi_s_n(:,1); phi_s_p(:,1); i_en(:,1); i_ep(:,1);...
          phi_e(:,1); jn(:,1); jp(:,1)];
end