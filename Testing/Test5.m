clear all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));

creaturesNumber = 10;
range = 100;
speed = 0.1;

positions = (range*rand(creaturesNumber, 2));
directions = zeros(creaturesNumber, 2);
targets = (range*rand(creaturesNumber, 2));
classes = zeros(creaturesNumber, 1);
classes(1,1) = 1;
classes = classes';

creature =  [0 0 0 0 0 0 0 0];
predator =  [0 0 0 0 0 0 0 1];
herbivore = [0 0 0 0 0 0 1 0];

D = [   predator herbivore;
        predator predator;
        herbivore herbivore;
        herbivore predator;];

labels = {'run', 'chase', 'search', 'eat', 'default'};

sD = som_data_struct(D);
sD = som_label(sD, 'add', 1, 'chase');
sD = som_label(sD, 'add', 2, 'search');
sD = som_label(sD, 'add', 3, 'search');
sD = som_label(sD, 'add', 4, 'run');

sM = som_make(sD,...,
    'init', 'lininit',...
    'algorithm', 'batch',...
    'lattice', 'hexa',...
    'shape', 'sheet',...
    'neigh', 'gaussian',...
    'training', 'short',...
    'mapsize', 'big',...
    'tracking', 3);
sM = som_autolabel(sM,sD, 'vote');

sM = som_supervised(sD, 'big');
[bmu, err] = som_bmus(sM, [predator herbivore], 'all');
bmu = bmu(find(min(err(:))));
co = sM.codebook(bmu,:);
sM.labels(bmu, :)