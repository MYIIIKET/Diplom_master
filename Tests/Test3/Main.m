clear all;
close all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));

%% Init 
sz = 1000;

%% Red
trainR = getR(sz);

%% Green
trainG = getG(sz);

%% Blue
trainB = getB(sz);

%%
trainY = getY(sz);

%%
trainM = getM(sz);

%%
trainC = getC(sz);

%% Train data 
trainData = [trainR;trainG;trainB;trainY;trainM;trainC];

%% Plot
% figure;
% hold on;
% for i=1:sz
%     plot3(testR(i,1),testR(i,2),testR(i,3),'.','color',testR(i,:));
%     plot3(testG(i,1),testG(i,2),testG(i,3),'.','color',testG(i,:));
%     plot3(testB(i,1),testB(i,2),testB(i,3),'.','color',testB(i,:));
% end

%% Train
sD = som_data_struct(trainData);
sD = som_label(sD,'add',[1:sz]','R');
sD = som_label(sD,'add',[sz+1:2*sz]','G');
sD = som_label(sD,'add',[2*sz+1:3*sz]','B');
sD = som_label(sD,'add',[3*sz+1:4*sz]','Y');
sD = som_label(sD,'add',[4*sz+1:5*sz]','M');
sD = som_label(sD,'add',[5*sz+1:6*sz]','C');

sM = som_supervised(sD,...
    'tracking', 3,...
    'algorithm', 'batch',...
    'mapsize', 'big',...
    'lattice', 'hexa',...
    'shape','sheet',...
    'neigh','gaussian');

% sM = som_label(sM,'clear','all');
% sM = som_autolabel(sM,sD);

som_show(sM,'umat','all');
% som_show_add('label',sM.labels,'TextSize',8,'TextColor','r')

%% 
figure; som_show(sM,'umat','all');

%% LVQ 
% testR = getR(sz);
% sD_test = som_data_struct(testR);
lab = unique(sD.labels(:,1));
mu = length(lab)*5;
sM = lvq1(sM, sD, 50*mu, 0.01);
figure; som_show(sM, 'umat','all');

%% 
% sM = lvq3(sM,sD,50*mu,0.01,0.25,0.1);
% figure; som_show(sM, 'umat','all');

%% Test
% testR = getC(sz);
% bmus_indxs_indata = som_bmus(sM, testR, 'best');
% bmus_indxs = unique(bmus_indxs_indata);
% 
% pos = 0;
% neg = 0;
% for i=1:length(bmus_indxs)
%     if strcmp(sM.labels{bmus_indxs(i,1),1},'C')
%         pos = pos+1;
%     else
%         neg = neg+1;
%     end
% end
% pos
% neg

