clear all;
close all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));
%% 
sz = 500;
D1 = rand(sz, 3);
m1 = find(abs(D1-mean(D1))<0.01);
[x1,y1]=ind2sub(size(D1),m1);

D2 = rand(sz, 3)/2.5+1;
m2 = find(abs(D2-mean(D2))<0.01);
[x2,y2]=ind2sub(size(D2),m2);

D = [D1;D2];

%% 

sM = som_randinit(D);
sD = som_data_struct(D);
%% 
sD = som_label(sD,'add',[1:500]','c1');
sD = som_label(sD,'add',[501:1000]','c2');
% 
% sM = som_make(sD,...
%     'tracking', 0,...
%     'algorithm', 'seq',...
%     'mapsize', 'big',...
%     'lattice', 'hexa',...
%     'shape','sheet',...
%     'neigh','gaussian');
% 
% sM.labels = cell(length(sM.codebook(:,1)),1);
% 
% bmus_indxs_indata = som_bmus(sM, D1, 'best');
% bmus_indxs1 = unique(bmus_indxs_indata);
% for i=1:length(bmus_indxs1)
%     sM.labels{bmus_indxs1(i,1),1}='c1';
% end
% 
% 
% bmus_indxs_indata = som_bmus(sM, D2, 'best');
% bmus_indxs2 = unique(bmus_indxs_indata);
% for i=1:length(bmus_indxs2)
%     sM.labels{bmus_indxs2(i,1),1}='c2';
% end
% 
% C = intersect(bmus_indxs1,bmus_indxs2)

%%
sM2 = som_make(sD,...
    'tracking', 1,...
    'algorithm', 'seq',...
    'mapsize', 'big',...
    'lattice', 'hexa',...
    'shape','sheet',...
    'neigh','gaussian');

%%
% sM3 = som_make(sD,...
%     'tracking', 0,...
%     'algorithm', 'seq',...
%     'mapsize', 'big',...
%     'lattice', 'hexa',...
%     'shape','sheet',...
%     'neigh','gaussian');
% 
% sM3 = som_label(sM3,'clear','all');
% sM3 = som_autolabel(sM3,sD);

%% 
lab = unique(sD.labels(:,1));
mu = length(lab)*5;
% sM2 = lvq1(sM2, sD, 50*mu, 0.01);
% figure; som_show(sM2, 'umat','all');

%%
sM = lvq3(sM,sD,50*mu,0.05,0.2,0.3);
figure; som_show(sM, 'umat','all');

%%
D3 = rand(sz, 3)+1;
bmus_indxs_indata = som_bmus(sM2, D3, 'best');
bmus_indxs = unique(bmus_indxs_indata);

pos = 0;
neg = 0;
for i=1:length(bmus_indxs)
    if strcmp(sM2.labels{bmus_indxs(i,1),1},'c1')
        pos = pos+1;
    else
        neg = neg+1;
    end
end
pos
neg