%% Exchange Current Density function i_0
%   Created July 12, 2011 by Scott Moura

function [i_0n,i_0p,varargout] = exch_cur_dens(p,c_ss_n,c_ss_p,c_e)
import casadi.*

% Parse out concentrations in anode and cathode
%%% CASADI CHANGE -- original spmet code uses average electrolyte
%%% concetration values in the anode, cathode, and separator then
%%% interpolates the values across a finer spatial discretization
% c_e_value = full(evalf(c_e)); %can't interpolate across SX variable so need to extract numerical value

% ce_interp=interp1(linspace(0,1,length(c_e)),c_e,linspace(0,1,p.Nx),'linear');
% c_e_n = ce_interp(1:p.Nxn-1);
% c_e_p = ce_interp(p.Nxn+p.Nxs-1:end);

c_e_n = c_e(1:p.Nxn-1); %%%% ZTG NOTE 2019-5-23
c_e_p = c_e(p.Nxn+p.Nxs-1:end);

% Compute exchange current density
i_0n = p.k_n * ((p.c_s_n_max - c_ss_n) .* c_ss_n .* c_e_n).^p.alph;
di_0n = p.k_n * (p.alph*((p.c_s_n_max - c_ss_n) .* c_ss_n .* c_e_n).^(p.alph-1).*( c_e_n.*(p.c_s_n_max - c_ss_n) -  c_ss_n .* c_e_n) );
i_0p = p.k_p * ((p.c_s_p_max - c_ss_p) .* c_ss_p .* c_e_p).^p.alph;
di_0p = p.k_p * (p.alph*(max((p.c_s_p_max - c_ss_p),0) .* c_ss_p .* c_e_p).^(p.alph-1).*( c_e_p.*(p.c_s_p_max - c_ss_p) -  c_ss_p .* c_e_p) );

if(nargout > 2)
    varargout{1}=di_0n;
    varargout{2}=di_0p;
end