clear all;
close all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));
d = 'C:\Users\MYlll\Desktop\Features\celloCustom13';

%% Load train data
% files = dir(d);
% names = {files.name};
% for i=3:length(names)
%     load(names{i});
% end
D = [F1';F2']';

%% Projection
tic
dim = 2;
proj = getProjection('sammon', D, dim);
% proj = D';
toc

%% Prepare data
sD = som_data_struct(D');
sD = som_normalize(sD, 'var');

sD = som_label(sD,'add',[1:length(F1)]','F1');
sD = som_label(sD,'add',[length(F1)+1:length(F1)+length(F2)]','F2');

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
sD = loadHits(sM, sD);

%% Visualization
figure;
som_show(sM,'umat','all');
showHits(sM, sD.H);

[q,t] = som_quality(sM,sD);
adm = som_distortion(sM,sD);
fprintf('QE = %f; TE = %f; ADM = %f;\n', q, t, adm);

figure;
scatter(sD.data(:,1), sD.data(:,2))