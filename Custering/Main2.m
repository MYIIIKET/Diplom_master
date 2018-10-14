clear all;
close all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));

%% Load train data
D = [cello'; viola'; violin']';
%% Projection
tic
dim = 2;
% proj = getProjection('sammon', D, dim);
proj = D';
toc
toReplace = 15;

%% Prepare data
sD = som_data_struct(proj);
% sD = som_normalize(sD, 'var');
sD.labels = [trainDataStruct1.som_data.labels; trainDataStruct2.som_data.labels];

%% Train map
% sM = som_make(sD,...,
%     'init', 'randinit',...
%     'algorithm', 'batch',...
%     'lattice', 'hexa',...
%     'shape', 'sheet',...
%     'neigh', 'gaussian',...
%     'training', 'long',...
%     'mapsize', 'big',...
%     'tracking', 3);
sM = som_supervised(sD, 'tracking', 3, 'big', 'msize', [50 50]);

%% Add hits 
sD = loadHits(sM, sD);

%% Visualization
visualize(sM, sD, sD.H)