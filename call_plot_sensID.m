%% call_plot_sensID
% Purpose: script for calling plot_sensID to generate figures related to
% parameter identification
% By: ZTG
% Date: 2019-7-25
clear all
close all
clc

%% Initialize I/O and other variables
num_batches = 5;
baseline = 'Baseline2b';
date = '14-Aug-2019 10_12_13';
results_folder = strcat('/Users/ztakeo/Box Sync/HPC/HPC2/',baseline,'/',date,'/');
% results_folder = strcat('/Users/ztakeo/Documents/GitHub/sens_ParamID/output-data/',baseline,'/',date,'/');
partial_path = strcat(results_folder,date,'_results_batch_');

%% get metrics and plot
% extract and compute metrics for plotting
metrics = compute_metrics(num_batches,partial_path);
save(strcat(results_folder,'metrics.mat'),'metrics');

% run plot_sensID
plot_sensID(metrics,results_folder,baseline,date)

