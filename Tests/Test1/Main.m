clear all;
close all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));
%% 
size = 1000;
D1 = rand(size, 3);
D2 = rand(size, 3)/2.5+1;
D = [D1;D2];
% D = rand(size, 3);
% D = table2array(Chainlink);

% plot3(D(:,1),D(:,2),D(:,3), '.');

sM1 = som_randinit(D);
sM2 = sM1;
sD = som_data_struct(D);
% sM = som_prototrain(sM,D, 20);

t = 100;
l = 100;
o = ones(l,1);
r = (1-(1:t)/t);
figure;
for i=1:t,
sM1 = som_seqtrain(sM1,sD,'tracking',0,...
    'trainlen',l,'samples',...
    'alpha',0.1*o,'radius',(4*r(i)+1)*o,...
    'con');
som_grid(sM1,'Coord',sM1.codebook)
%   hold on, plot3(D(:,1),D(:,2),D(:,3),'.'), hold off
hold on, plot(D(:,1),D(:,2),'.'), hold off
drawnow
end
figure; som_show(sM1, 'umat','all');
[qe, te]=som_quality(sM1, sD)
bmus_indxs_indata = som_bmus(sM1, sD, 'best');
bmus_indxs = length(unique(bmus_indxs_indata))

figure;
for i=1:t,
sM2 = som_seqtrain(sM2,sD,'tracking',0,...
    'trainlen',l,'samples',...
    'alpha',0.1*o,'radius',(4*r(i)+1)*o);
som_grid(sM2,'Coord',sM2.codebook)
%   hold on, plot3(D(:,1),D(:,2),D(:,3),'.'), hold off
hold on, plot(D(:,1),D(:,2),'.'), hold off
drawnow
end
figure; som_show(sM2, 'umat','all');
[qe, te]=som_quality(sM2, sD)
bmus_indxs_indata = som_bmus(sM2, sD, 'best');
bmus_indxs = length(unique(bmus_indxs_indata))