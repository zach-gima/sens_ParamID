%% Nonlinear output for voltage in Single Particle Model
%   Created July 21, 2011 by Scott Moura

function [V] = nonlinearSPMOutputVoltage(p,c_se_n,c_se_p,I,VCE)

% Stochiometric Concentration Ratio
theta_n = c_se_n / p.c_s_n_max;
theta_p = c_se_p / p.c_s_p_max;

% Equilibrium Potential
Unref = refPotentialAnode(p,theta_n);
Upref = refPotentialCathode(p,theta_p);

% Exchange Current Density
c_e0 = p.c_e * ones(p.Nx,1);     % Fixed electrolyte concentration [mol/m^3]
[i_0n,i_0p] = exch_cur_dens(p,c_se_n,c_se_p,c_e0);

V = (p.R*p.T_amp)/(p.alph*p.Faraday) * asinh(I / (2*p.a_s_p*p.Area*p.L_p*i_0p(1))) ...
    -(p.R*p.T_amp)/(p.alph*p.Faraday) * asinh(I / (2*p.a_s_n*p.Area*p.L_n*i_0n(1))) ...
    + Upref - Unref - (p.R_f_n/(p.a_s_n*p.L_n*p.Area) + p.R_f_p/(p.a_s_p*p.L_p*p.Area))*I ...
    - (p.L_p + 2*p.L_s + p.L_n)*I + VCE;