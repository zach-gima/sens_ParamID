function [kappa,varargout] = ElectrolyteCond(T,c)

% DFN Electrolyte Conductivity
% Created by Dong Zhang
% Source: Transport Properties of LiPF6-Based Li-Ion Battery Electrolytes -
% Valoen and Reimers
% 5/10/2018

% T is in Kelvin
% c is in mol/L

% Transfer unit from mol/m^3 to mol/L
c = c/1000;

% Parameters from Table III
k = zeros(3,3);
k(1,1) = -10.5;
k(1,2) = 0.0740;
k(1,3) = -6.96*10^(-5);
k(2,1) = 0.668;
k(2,2) = -0.0178;
k(2,3) = 2.80*10^(-5);
k(3,1) = 0.494;
k(3,2) = -8.86*10^(-4);
k(3,3) = 0;

% Kappa from Eq. [17]
K = 0;
for i = 0:2
    for j = 0:2
        K = K + k(i+1,j+1).*c.^i.*T.^j;
    end
end

kappa = c.*K.^2;

% Compute derivative of K wrt c
dK = k(2,1) + k(2,2).*T + k(2,3).*T.^2 + 2*(k(3,1).*c + k(3,2).*c.*T + k(3,3).*c.*T.^2);

if(nargout == 2)
    dkappa = 2.*K.*dK.*c + K.^2;
    varargout{1} = dkappa;
end
