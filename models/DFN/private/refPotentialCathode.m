%% Reference Potential for Pos. Electrode: Unref(theta_p)
%   Created May 5, 2017 by Saehong Park & Dylan Kato
%   NCA cell, (NCR-18650pf) 
%   Curve-fitting w/ Bosch OCP data

function [Uref,varargout] = refPotentialCathode(theta)

c_p=[ -40.045585568588542
 -62.042811084183654
  52.447046217508564
 -11.281882678497276
  63.276043910291172
  21.159687366489823
  37.390508845897301
  22.261671639629835
   8.574181451931103
  10.133001156239731
  -3.313604725236584
   1.977856101799057
  -3.046181118828750
  -0.087883198680019
  -0.836818408057937
  -0.072435003409969
  -0.069320106210594
   4.456159792325947];

    w=c_p(end);
    
    Uref=c_p(1) + c_p(2).*cos(theta.*w) + c_p(3).*sin(theta.*w) + ...
      c_p(4).*cos(2.*theta.*w) + c_p(5).*sin(2.*theta.*w) + c_p(6).*cos(3.*theta.*w) + c_p(7).*sin(3.*theta.*w) + ...
      c_p(8).*cos(4.*theta.*w) + c_p(9).*sin(4.*theta.*w) + c_p(10).*cos(5.*theta.*w) + c_p(11).*sin(5.*theta.*w) + ...
      c_p(12).*cos(6.*theta.*w) + c_p(13).*sin(6.*theta.*w) + c_p(14).*cos(7.*theta.*w) + c_p(15).*sin(7.*theta.*w) + ...
      c_p(16).*cos(8.*theta.*w) + c_p(17).*sin(8.*theta.*w);

  % Gradient of OCP wrt theta
if(nargout == 2)
    
%     % Polynomail Fit
%     dUref = ppvalFast(p.dUppp,theta);
%     varargout{1} = dUref / p.c_s_p_max;

 dUref = -c_p(2).*w.*sin(w .*theta) - 2.*c_p(4).*w .*sin(2 .*w .*theta) - 3.* c_p(6).* w.* sin(3 .*w .*theta) - ...
        4 .*c_p(8) .*w .*sin(4 .*w .*theta) - 5 .*c_p(10) .*w .*sin(5 .*w .*theta) - 6 .*c_p(12) .*w .*sin(6 .*w .*theta) - ...
        7.* c_p(14) .*w .*sin(7.* w.* theta) - 8 .*c_p(16) .*w .*sin(8 .*w.* theta)+c_p(1) .*w.* cos(w.* theta) + ...
        2.* c_p(3).* w .*cos(2.* w.* theta) + 3 .*c_p(5) .*w .*cos(3 .*w .*theta) + ...
        4.* c_p(7) .*w .*cos(4.* w .*theta) + 5.* c_p(9) .*w .*cos(5 .*w .*theta) + 6.* c_p(11).* w.* cos(6.* w .*theta) + ...
        7.* c_p(13).* w .*cos(7.* w.* theta) + 8 .*c_p(15).* w.* cos(8 .*w .*theta);

varargout{1} = dUref;
    
end

% Gradient of OCP wrt temperature
if(nargout >= 3)
    
    dUdT = 0;
    varargout{2} = dUdT;
    
end