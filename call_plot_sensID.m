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
results_folder = '/Users/ztakeo/Box Sync/HPC/HPC2/BaselineA/01-Aug-2019 10_49_13/';
partial_path = strcat(results_folder,'01-Aug-2019 10_49_13_results_batch_');

%% get metrics and plot
% extract and compute metrics for plotting
metrics = compute_metrics(num_batches,partial_path);
save(strcat(results_folder,'metrics.mat'),'metrics');

% run plot_sensID
plot_sensID(metrics,partial_path)

