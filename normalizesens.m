%%% Created by Zach Gima; Adapted from SHP origin_to_norm code

function [normalize_sens_factor] = normalizesens(p,bounds)
    % % % % Log normalization
    % % % % Multiply (1/log10(exp(1))) * log(\theta_max / \theta_min) * theta_0

    % % % % minmax normalization
    % % % % Multiply (\theta_max - \theta_min)

    % % % % nominal range normalization
    % % % % Multiply norminal param (\theta_0)

    % % % % Parameter Min/Max values

    
     %emptyblk = NaN;
%      normalize_sens_factor = [ (1/log10(exp(1)))*log10(bounds.max(1) / bounds.min(1))*p.D_s_n0; ... %(1) Dsn (log)
%                 (1/log10(exp(1)))*log10(bounds.max(2) / bounds.min(2))*p.D_s_p0; ... %(2) Dsp (log)
%                 bounds.max(3) - bounds.min(3);%(3) Rsn (Min/Max)
%                 bounds.max(4) - bounds.min(4);%(4) Rsp (Min/Max)
% %                 emptyblk;%(5)
% %                 emptyblk;%(6)
%                 (1/log10(exp(1)))*log10(bounds.max(7)/bounds.min(7)) * p.sig_n; ... %(7) sig_n (log)
%                 (1/log10(exp(1)))*log10(bounds.max(8)/bounds.min(8)) * p.sig_p; ... %(8) sig_p (log)
%                 (bounds.max(9) - bounds.min(9)); ... %(9) D_e (Min/Max) 
%                 bounds.max(10) - bounds.min(10); ... %(10) eps_e_n (Min/Max)
%                 bounds.max(11) - bounds.min(11); ... %(11) eps_e_s (Min/Max)
%                 bounds.max(12) - bounds.min(12); ... %(12) eps_e_p (Min/Max)
%                 (bounds.max(13) - bounds.min(13)) ; ... %(13) kappa (Min/Max)
%                 (bounds.max(14) - bounds.min(14)) ; ... %(14) t_plus (Min/Max) 
%                 (bounds.max(15) - bounds.min(15)); ... %(15) dactivity (Min/Max)
%                 (1/log10(exp(1)))*log10(bounds.max(16) / bounds.min(16)) * p.k_n0; ... %(16) k_n (log)
%                 (1/log10(exp(1)))*log10(bounds.max(17) / bounds.min(17)) * p.k_p0; ... %(17) k_p (log)
%                 bounds.max(18) - bounds.min(18); ... %(18) Rfn (Min/Max)
%                 bounds.max(19) - bounds.min(19); ... %(19) Rfp (Min/Max)
% %                 emptyblk; %(20)
%                 bounds.max(21) - bounds.min(21); %(21) ce_0 (Min/Max)                     
%                 bounds.max(22) - bounds.min(22); %(22) E.Dsn (Min/Max)               
%                 bounds.max(23) - bounds.min(23); %(23) E.Dsp (Min/Max)                
%                 bounds.max(24) - bounds.min(24); %(24) E.kn (Min/Max)              
%                 bounds.max(25) - bounds.min(25); %(25) E.kp (Min/Max)               
%                 ];

    normalize_sens_factor = [bounds.max(1) - bounds.min(1);%(4) Rsp (Min/Max)
                bounds.max(2) - bounds.min(2); ... %(15) dactivity (Min/Max)
                bounds.max(3) - bounds.min(3); ... %(10) eps_e_n (Min/Max)
                bounds.max(4) - bounds.min(4); ... %(14) t_plus (Min/Max) 
                bounds.max(5) - bounds.min(5); ... %(18) Rfn (Min/Max)
                bounds.max(6) - bounds.min(6); ... %(19) Rfp (Min/Max)            
                ];
            
    normalize_sens_factor = normalize_sens_factor'; % row maxtrix

end