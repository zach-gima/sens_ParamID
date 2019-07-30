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
results_folder = '/Users/ztakeo/Box Sync/HPC/HPC2/BaselineA/24-Jul-2019 22_42_43/';
partial_path = strcat(results_folder,'24-Jul-2019 22_42_43_results_batch_');

%% get metrics and plot
% extract and compute metrics for plotting
metrics = compute_metrics(num_batches,partial_path);

% run plot_sensID
plot_sensID(metrics,partial_path)

