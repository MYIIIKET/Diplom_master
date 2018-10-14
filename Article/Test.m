clear all;
clc;
addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));

n=1000;
dim = 3;
data = rand(n, dim);

% hold on;
% for i=1:length(data)
%     plot(data(i,1),data(i,2),'.', 'Color', [data(i,1) data(i,2) data(i,3)]);
% end
data(1,:)=[1 0 0];
data(2,:)=[0 1 0];
data(3,:)=[0 0 1];
sD = som_data_struct(data);
sD = som_label(sD,'add',[1]','r');
sD = som_label(sD,'add',[2]','g');
sD = som_label(sD,'add',[3]','b');
sM = som_supervised(sD,...
    'tracking', 3,...
    'algorithm', 'batch',...
    'mapsize', 'big',...
    'lattice', 'hexa',...
    'shape','sheet',...
    'neigh','bubble',...
    'training','long');

%% 
C = som_colorcode(sM,'rgb4');
som_cplane(sM,C);

%%
figure;
colormap(1-gray);
som_show(sM,'umat','all');
som_show_add('label',sM.labels,'TextSize',8,'TextColor','r')