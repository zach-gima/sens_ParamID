%% call_spmet_casadi
% ZTG: meant for testing spmet simulation against simulink CCCV controller implementation 
clear all
close all
clc
% load necessary parameters and data
run params/params_NCA;

% Add path to model of choice: SPMeT or DFN
addpath(genpath('models/'))

%%% Set Model Discretization Parameters
p = set_discretization(p,'SPMeT');

t_length = 3500;
dt = 1;

Current = -2.9; %1.45 = 0.5C

data.cur = Current*ones(t_length,1);
data.cur(1) = 0;
data.V0 = 2.8;
data.time = (0:dt:(dt*t_length-dt))';
data.T_amb = 25+273.15;

% Call function
[V,state_info] = spmet_casadi(p,data);

% plot stuff
plot(V);