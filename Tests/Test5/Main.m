close all;
clear all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));

%%
load chainlink.mat chainlink
data_struct = struct();
firs_chain = chainlink(1:500,:);
second_chain = chainlink(501:1000,:);
data_struct.data{1}=firs_chain;
data_struct.data{2}=second_chain;

%% 
sD = som_data_struct(chainlink);

sM = som_randinit(sD);
sM = som_seqtrain(sM,sD,...
    'tracking', 3,...
    'trainlen',100,...
    'con_thresh',1,...
    'con_rest',2,...
    'delta',1,...
    'con');
figure;
som_grid(sM,'Coord',sM.codebook);