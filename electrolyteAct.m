function [Act,varargout] = electrolyteAct(c_e,T,p)

% (d ln f+-)/(d ln c_e) 
% DFN thermal factor as function of temperature and electrolyte
% concentration
% v(c,T) = (1 - t_+)*(1+ (d ln f+-)/(d ln c_e)) - Eq. [3] in paper
% Source: Transport Properties of LiPF6-Based Li-Ion Battery Electrolytes -
% Valoen and Reimers
% Created by Dong Zhang
% 5/11/2018

% T is in Kelvin
% c is in mol/L
% v is unit-less

c_e = c_e/1000;

v00 = 0.601; v01 = 0; 
v10 = -0.24; v11 = 0;
v20 = 0; v21 = 0;
v30 = 0.982; v31 = -0.0052; % Table I

v0 = v00*(1+v01*(T-p.T_ref));
v1 = v10*(1+v11*(T-p.T_ref));
v2 = v20*(1+v21*(T-p.T_ref));
v3 = v30*(1+v31*(T-p.T_ref));

v = v0 + v1*c_e.^0.5 + v2*c_e + v3*c_e.^(1.5); % Eq. [4]
Act = v./(1-p.t_plus) - 1; % Eq. [3]

if(nargout == 2)
    dAct = (0.5*v10*(c_e).^(-0.5)+(1.5)*v30*(1+v31*(T-p.T_ref))*(c_e).^(0.5))/(1-p.t_plus);
    varargout{1} = dAct;
end