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
baseline = 'BaselineC';
date = '07-Aug-2019 12_44_36';
results_folder = strcat('/Users/ztakeo/Box Sync/HPC/HPC2/',baseline,'/',date,'/');
partial_path = strcat(results_folder,date,'_results_batch_');

%% get metrics and plot
% extract and compute metrics for plotting
metrics = compute_metrics(num_batches,partial_path);
save(strcat(results_folder,'metrics.mat'),'metrics');

% run plot_sensID
plot_sensID(metrics,results_folder,baseline,date)

