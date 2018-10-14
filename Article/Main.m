clear all;
close all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));
addpath(genpath('C:\Users\MYlll\Desktop\Features\flatData'));
addpath(genpath('C:\Users\MYlll\Desktop\Features\out'));

%% 
load m1_train m1_train;
load f1_train f1_train;
load m1_test m1_test;
load f1_test f1_test;

%% Load train data 1
trainDataStruct1.flat_data = f1_train;

%% Load train data 2
trainDataStruct2.flat_data = m1_train;

%% Load test data 1
testDataStruct1.flat_data = f1_test;

%% Load test data 2
testDataStruct2.flat_data = m1_test;

%% Unite train data
D = [m1_train'; f1_train';]';
D = som_normalize(D', 'range');
sD = som_data_struct(D);

%% 
a = length(m1_train);
b = a+length(f1_train);
% c = b+length(trainDataStruct3.flat_data);
hold on;
plot3(D(1:a,1),D(1:a,2),D(1:a,3),'.r');
plot3(D(a+1:b,1),D(a+1:b,2),D(a+1:b,3),'.g');
% plot3(D(b+1:c,1),D(b+1:c,2),D(b+1:c,3),'.b');
hold off;

%% Add labels
a = length(m1_train);
b = a+length(f1_train);
c = b+length(f1_test);
sD = som_label(sD,'add',[1:a]','m');
sD = som_label(sD,'add',[a+1:b]','f');
% sD = som_label(sD,'add',[b+1:c]','f');
% sD = som_label(sD,'add',[b+1:c]','f');

lab = unique(sD.labels(:,1));
mu = length(lab)*5;

%% Init map
sM = som_randinit(sD);
%% 
for i=3:length(sM.codebook(:,1))
    if mod(i,2)==0
        lab{end+1,1}='f';
    else
        lab{end+1,1}='m';
    end
    
end
sM.labels = lab;
sM = lvq1(sM,sD,50*mu,0.05);
sM = lvq3(sM,sD,50*mu,0.05,0.2,0.3);
%%
sM = som_supervised(sD,...
    'tracking', 3,...
    'algorithm', 'batch',...
    'mapsize', 'small',...
    'lattice', 'hexa',...
    'shape','sheet',...
    'neigh','gaussian');

%% labeling
sM = sM;
sM = som_label(sM,'clear','all');
sM = som_autolabel(sM,sD);

%% LVQ 
lab = unique(sD.labels(:,1));
mu = length(lab)*5;
sD2 = som_data_struct(f1_test');
sM = lvq1(sM, sD2, 50*mu, 0.01);
figure; som_show(sM, 'umat','all');
%% 
sM = lvq3(sM,sD,50*mu,0.05,0.2,0.3);
figure; som_show(sM, 'umat','all');

%%
s1 = sD.data(1:a,:);
s2 = sD.data(a+1:b,:);
h1 = som_hits(sM, s1);
h2 = som_hits(sM, s2);
figure;
colormap(1-gray);
som_show(sM,'umat','all');
% som_show_add('hit',[h1, h2],'MarkerColor',[1 0 0; 0 1 0],'Subplot',1);

%%
figure;
h = som_hits(sM, sD);
colormap(1-gray);
som_show(sM,'umat','all');
som_show_add('label',sM.labels,'TextSize',8,'TextColor','r')

%%
[Pd,V,me] = pcaproj(D,2);
Pm = pcaproj(sM.codebook,V,me);
C = som_colorcode(Pm);  
figure;
som_show(sM,'umat','all');
som_cplane(sM,C);

%% 
sD2 = som_data_struct(m1_test');
% sD2 = som_label(sD,'clear','all');
sD2 = som_autolabel(sD2,sM);

labels_train = sD.labels(b+1:end,1);
labels_test = sD2.labels;

sz = length(sD2.labels);
(sum(strcmp(labels_test,'f'))/sz)*100
(sum(strcmp(labels_test,'m'))/sz)*100

% ok = strcmp(labels_train,labels_test);
% 100*(1-sum(ok)/length(ok))

%%
D_test = m1;
D_test = som_normalize(D_test','range');
sD_test = som_data_struct(D_test);
% sD_test = som_label(sD_test,'add',[1:length(D_test(:,1))]','f');

%% 
hf0 = som_hits(sM, sD, 'fuzzy');

som_show(sM,'umat','all');
som_show_add('hit',hf0,'Subplot',1,'MarkerColor','r');
%% 

[bmus qerr] = som_bmus(sM, sD_test.data);

sum(strcmp(sM.labels(bmus),'f'))
sum(strcmp(sM.labels(bmus),'m'))

%%
sD_test = som_label(sD_test,'clear','all');

%% 
sD_test = som_autolabel(sD_test,sM);
sum(strcmp(sD_test.labels,'f'))
sum(strcmp(sD_test.labels,'m'))
% ok = strcmp(sD_test.labels,sM.labels);
% 100*(1-sum(ok)/length(ok))

%%
len = length([b+1:c]);
(sum(strcmp(sD2.labels(b+1:c),sD.labels(b+1:c)))/len)

%%
D3 = m1_test';
bmus_indxs_indata = som_bmus(sM, D3, 'all');
bmus_indxs = unique(bmus_indxs_indata)

pos = 0;
neg = 0;
for i=1:length(bmus_indxs)
    if strcmp(sM.labels{bmus_indxs(i,1),1},'m')
        pos = pos+1;
    else
        neg = neg+1;
    end
end
pos
neg