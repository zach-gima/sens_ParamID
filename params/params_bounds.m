%% Bounds for the "Nominal" Parameters from Literature Used to Initialize ParamID by Saehong Park in Bosch-UC Berkeley Yr. 1
% Corresponds to NCA cell (param/params_NCA_final.m)
% By: Zach Gima 2018-8-9
% Version: 3 (Updated transference number bounds and slimmed down for online paramID 2019-4-19)
%% May/may not need to add in the 2-state thermal dynamics parameters into this

p_bounds.min = [  1.0e-6;   % R_s_p
                0.5;          %15.d_activity
                0.18;        %10.eps_e_n
                0.1;       %14.t_plus
                1e-5;      %18.R_f_n
                1e-5;      %19.R_f_p         ZTG NOTE originally set at 1e-4; expanding
                ];

             
p_bounds.max = [  100.0e-6;   %R_s_p
                1.5;          %d_activity
                0.45;        %eps_e_n                
                0.8;       %t_plus
                1e-3;      %R_f_n
                1e-3;      %R_f_p        
                ];
% 
% bounds.min = [2.25E-16;    %1.D_s_n
%                  2.00E-16;    %2.D_s_p
%                  1.0e-6;   %3.R_s_n
%                  1.0e-6;   %4.R_s_p
%                  NaN;%nan;        %5.eps_s_n (x), Equil. Struct.
%                  NaN;%nan;        %6.eps_s_p (x), Equil. Struct.
%                  50;        %7.sig_n
%                  50;        %8.sig_p
%                  0.5;          %9.D_e
%                  0.18;        %10.eps_e_n
%                  0.45;        %11.eps_e_s
%                  0.18;        %12.eps_e_p
%                  0.5;          %13.Kappa
%                  0.1;       %14.t_plus
%                  0.5;          %15.d_activity
%                  0.1*7.5e-04;    %16.k_n0
%                  0.1*2.3e-03;    %17.k_p0
%                  1e-5;      %18.R_f_n
%                  1e-4;      %19.R_f_p
%                  NaN;%nan;        %20.n_Li_s (x), Equil. Struct.
%                  500;         %21.ce0
%                  22e3;       %22.E.Dsn
%                  20e3;       %23.E.Dsp
%                  37.48e3;       %24.E.kn
%                  30e3        %25.E.kp                 
%                  ];  % assign nominal vlaues that selects sensitivity
% 
%              
% bounds.max = [1.05E-12;    %1.D_s_n
%                  1.00E-12;    %2.D_s_p
%                  100.0e-6;   %3.R_s_n
%                  100.0e-6;   %4.R_s_p
%                  NaN;%nan;        %5.eps_s_n (x), Equil. Struct.
%                  NaN;%nan;        %6.eps_s_p (x), Equil. Struct.
%                  500;        %7.sig_n
%                  500;        %8.sig_p
%                  1.5;          %9.D_e
%                  0.45;        %10.eps_e_n
%                  0.5;        %11.eps_e_s
%                  0.33;        %12.eps_e_p
%                  1.5;          %13.Kappa
%                  0.8;       %14.t_plus
%                  1.5;          %15.d_activity
%                  10*7.5e-04;    %16.k_n0
%                  10*2.3e-03;    %17.k_p0
%                  1e-3;      %18.R_f_n
%                  1e-3;      %19.R_f_p
%                  NaN;%nan;        %20.n_Li_s (x), Equil. Struct.
%                  1500;         %21.ce0
%                  48.9e3;       %22.E.Dsn
%                  80.6e3;       %23.E.Dsp
%                  67.5e3;       %24.E.kn
%                  43.6e3        %25.E.kp                     
%                  ];  % assign nominal vlaues that selects sensitivity
             