clear all;
close all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));

D1 = rand(500, 2)+0.2;
D2 = rand(500, 2)-0.2;

% figure;
% hold on;
% plot(D1(:,1), D1(:,2), 'ro');
% plot(D2(:,1), D2(:,2), 'bo');
% hold off;

D = [D1; D2];

sD = som_data_struct(D);
sD = som_normalize(sD, 'var');

sM = som_make(sD,...,
    'init', 'lininit',...
    'algorithm', 'seq',...
    'lattice', 'hexa',...
    'shape', 'sheet',...
    'neigh', 'gaussian',...
    'training', 'long',...
    'mapsize', 'big',...
    'tracking', 3);

%%
Dt = 2*rand(500, 2)-1;

sDt = som_data_struct(Dt);
% sDt = som_normalize(sDt, 'var');

[bmu err] = som_bmus(sM, sDt, 'all');
bmu = bmu(find(min(err(:))));
co = sM.codebook(bmu,:);


figure;
hold on;
som_grid(sM,'Coord',sM.codebook,...
    'Markersize',2,'Linecolor','k','Surf',sM.codebook(:,2));
plot(co(1),co(2), 'or');
plot(Dt(:,1), Dt(:,2), 'b+');
hold off;

%%
