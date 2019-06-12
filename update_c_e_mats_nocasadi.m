% By: Zach Gima 2019-4-30
% Function that updates the electrolyte concentration matrices (c_e_mats)
% according to the parameters. Function was written in the context of
% parameter estimation, where re-identifying some electrolyte-related
% parameters would consequently require these matrices (stored in the p struct) to be updated
function p = update_c_e_mats_nocasadi(p)
%     [M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p);
    [M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats_nocasadi(p);
    p.ce.M1n = M1n;
    p.ce.M2n = M2n;
    p.ce.M3n = M3n;
    p.ce.M4n = M4n;
    p.ce.M5n = M5n;

    p.ce.M1s = M1s;
    p.ce.M2s = M2s;
    p.ce.M3s = M3s;
    p.ce.M4s = M4s;

    p.ce.M1p = M1p;
    p.ce.M2p = M2p;
    p.ce.M3p = M3p;
    p.ce.M4p = M4p;
    p.ce.M5p = M5p;

    p.ce.C = C_ce;

    clear M1n M2n M3n M4n M5n M1s M2s M3s M4s M1p M2p M3p M4p M5p C_ce;
end