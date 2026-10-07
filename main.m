%% OH - Satellite Image Optimization
% Optimization and Heuristics
% Sentinel-2 Satellite Image Processing
%
% Main script of the project.

clear;
clc;
close all;

fprintf('============================================\n');
fprintf('  SATELLITE IMAGE OPTIMIZATION PROJECT\n');
fprintf('============================================\n\n');

%% Project paths

projectRoot = fileparts(mfilename('fullpath'));

addpath(fullfile(projectRoot, 'src'));

rawDataDir       = fullfile(projectRoot, 'data', 'raw');
processedDataDir = fullfile(projectRoot, 'data', 'processed');
figuresDir       = fullfile(projectRoot, 'results', 'figures');
metricsDir       = fullfile(projectRoot, 'results', 'metrics');

fprintf('Project initialized successfully.\n');
fprintf('Raw data: %s\n', rawDataDir);
fprintf('Processed data: %s\n', processedDataDir);
fprintf('Figures: %s\n', figuresDir);
fprintf('Metrics: %s\n\n', metricsDir);
