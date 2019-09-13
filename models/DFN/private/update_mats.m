% By: Zach Gima 2019-9-11
% Function that updates the electrolyte phase concentration
% matrices (c_e_mats), solid & electrolyte phase electric
% potentials (phi_s_mats, phi_e_mats), and elyte current (i_e_mats)
% according to the parameters. Function was written in the context of
% parameter estimation, where re-identifying some electrolyte-related
% parameters would consequently require these matrices (stored in the p struct) to be updated
function p = update_mats(p,SensFlag)
    
    % Electrolyte concentration matrices
    [M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p,SensFlag);

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
    
    % ZTG Change
    Nn = p.Nxn - 1;
    Ns = p.Nxs - 1;
    Np = p.Nxp - 1;
    
    rM3 = [Nn; Ns; Np];
    cM3 = rM3';
    p.ce.M3 = sparse(blkdiagFast(rM3, cM3, p.ce.M3n, p.ce.M3s, p.ce.M3p));

    clear M1n M2n M3n M4n M5n M1s M2s M3s M4s M1p M2p M3p M4p M5p C_ce rM1 cM1;

    % Solid Electric Potential
    [F1_psn,F1_psp,F2_psn,F2_psp,G_psn,G_psp,...
        C_psn,C_psp,D_psn,D_psp] = phi_s_mats(p);
    p.F1_psn = F1_psn;
    p.F1_psp = F1_psp;
    p.F2_psn = F2_psn;
    p.F2_psp = F2_psp;
    p.G_psn = G_psn;
    p.G_psp = G_psp;
    p.C_psn = C_psn;
    p.C_psp = C_psp;
    p.D_psn = D_psn;
    p.D_psp = D_psp;

    clear F1_psn F1_psp F2_psn F2_psp G_psn G_psp C_psn C_psp D_psn D_psp;

    % Electrolyte Current
    [F1_ien,F1_iep,F2_ien,F2_iep,F3_ien,F3_iep] = i_e_mats(p,SensFlag);
    p.F1_ien = F1_ien;
    p.F1_iep = F1_iep;
    p.F2_ien = F2_ien;
    p.F2_iep = F2_iep;
    p.F3_ien = F3_ien;
    p.F3_iep = F3_iep;

    clear F1_ien F1_iep F2_ien F2_iep F3_ien F3_iep;

    % Electrolyte Electric Potential
    p.M1_pen_skel = sparse(diag(ones(p.Nxn-2,1),1) + diag(-ones(p.Nxn-2,1),-1));
    p.M1_pes_skel = sparse(diag(ones(p.Nxs-2,1),1) + diag(-ones(p.Nxs-2,1),-1));
    p.M1_pep_skel = sparse(diag(ones(p.Nxp-2,1),1) + diag(-ones(p.Nxp-2,1),-1));

    [M1_pe,M2_pe,M3_pe,M4_pe,C_pe] = phi_e_mats(p);
    p.M1_pe = M1_pe;
    p.M2_pe = M2_pe;
    p.M3_pe = M3_pe;
    p.M4_pe = M4_pe;
    p.C_pe = C_pe;

    clear M1_pe M2_pe M3_pe M4_pe C_pe
end