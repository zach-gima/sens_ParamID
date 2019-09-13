function [D_e,varargout] = electrolyteDe(T,c_e)

% DFN Electrolyte Diffusion Coefficient
% Created by Dong Zhang
% Source: Transport Properties of LiPF6-Based Li-Ion Battery Electrolytes -
% Valoen and Reimers
% 4/24/2018

% T is in Kelvin
% c is in mol/L - not mol/m^3
% De is in cm^2/s

% Transfer unit from mol/m^3 to mol/L
c_e = c_e/1000;

% Parameters from Table II
D00 = -4.43;
D01 = -54;
D10 = -0.22;
D11 = 0;
Tg0 = 229;
Tg1 = 5.0;

Tg = Tg0 + c_e*Tg1;
D0 = D00 + D01./(T-Tg);
D1 = D10 + D11./(T-Tg);

D_e = 10^(-4)*10.^(D0+D1.*c_e);

if(nargout == 2)
    dD_e = 10^(-7)*(10.^(D00 - D01./(Tg0 - T + Tg1.*c_e) + c_e.*(D10 - D11./(Tg0 - T + Tg1.*c_e))).*log(10).*(D10 - D11./(Tg0 - T + Tg1.*c_e) + (D01.*Tg1)./(Tg0 - T + Tg1.*c_e).^2 + (D11.*Tg1.*c_e)./(Tg0 - T + Tg1.*c_e).^2));
    varargout{1} = dD_e;
end