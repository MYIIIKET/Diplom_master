clear all; close all; %echo off; %clc;
disp('_______________________________________________________________');
data_path_train='C:\Users\123\Documents\MATLAB'; cd(data_path_train);
% �������� � ���������� ������� ������ %%
load('Chainlink.mat');
[Dlen, Ddim]=size(data);
D_init_all=data(:,[1:Ddim-1,Ddim+1:end]); % ������� ������� � �������
Ddim=Ddim-1; % ����������� ������������ ���������
sD_init=som_data_struct([D_init_all]);
%% ��������� �����:
mShape = 'sheet'; % sheet/toroid
mWidth=64; mHeight=64; 
mUnits=mWidth*mHeight; %mUnits=Dlen; %mUnits = round((5*Dlen^0.54321)); F=factor(mUnits);
% ������:
R00 = sqrt((mWidth/2)^2+(mHeight/2)^2)/3;%mWidth/6;
R01 = 1.5;
R11 = R01;
% �������� ��������:
LR0 = 0.8;
LR1 = LR0;
% ����������������� ��������:
trainlen1 = Dlen+10000;  %length(D_init_all)%(msizeX*msizeY)
trainlen2 = Dlen+100000; %500*length(D_init_all)%(msizeX*msizeY)
fprintf('Shape: %s\nNeurons: %d\nRadius: %.2f %.2f %.2f\nLR: %.2f %.2f\n',...
    mShape,mUnits,R00,R01,R11,LR0,LR1);
% ��������� ���������:
global conThreshold_init; conThreshold_init = 2; %1
global actFreq_init;      actFreq_init = 3;      %3
global actDelta_init;     actDelta_init = 0.1;   %0
global conRest_init;      conRest_init = 4; %2
global delta_init;        delta_init = 0.1;      %0.01
%% ���������� ����������� ������� ��� ������� ������ %%
%{
figure(1); clf;
%
tic; disp('2D Sammon projection of input data');
projD2d = sammon(sD_init.data, 2, 10^(-5), '3'); toc; % 2D
sD2d = som_data_struct([projD2d]);
%
figure(2); clf;
tic; disp('3D Sammon projection of input data');
projD3d = sammon(sD_init.data, 3, 10^(-5), '3'); toc; % 3D
sD3d = som_data_struct([projD3d]);
%}
%% ��������� SOM %
% 1) ������ ����������:
%sD = sD2d;
%sD = sD3d;
tic; figure(2); clf; %title('Fig.2');
sTopol = som_set('som_topol', 'msize',[mWidth mHeight], ...
    'lattice','hexa', 'shape',mShape);
sM_init = som_lininit(sD_init, sTopol); %randinit/lininit
sTrain1 = som_train_struct(sM_init, sTopol, 'rough');
sTrain1 = som_set(sTrain1, 'algorithm','seq', 'neigh','gaussian', ...
    'radius_ini',R00, 'radius_fin',R01, ...
    'alpha_ini',LR0, 'alpha_type','power',...
    'trainlen',trainlen1);
disp('> Rough training phase...'); clf;
fprintf('Training length (rough): %d\n',sTrain1.trainlen);
%
sM_train_rough = som_seqtrain(sM_init, sD_init, sTopol, sTrain1, ...
    'trainlen_type','samples', ... % samples/epochs
    'sample_order','random', ... % random/ordered
    'tracking',3); toc;
%[qe, te] = som_quality(sM_train, sD); fprintf('QE = %f, TE = %f\n',qe,te);

% 2) ������ ����������:
tic; figure(2); clf;
sTrain2 = som_train_struct(sM_train_rough, sTopol, 'finetune');
sTrain2 = som_set(sTrain2, 'algorithm','seq', 'neigh','gaussian', ...
    'radius_ini',R01, 'radius_fin',R11, ...
    'alpha_ini',LR1, 'alpha_type','power',...
    'trainlen',trainlen2);
disp('> Finetuning...');
fprintf('Training length (finetune): %d\n',sTrain2.trainlen);
fprintf('Training length (total): %d\n',(sTrain1.trainlen+sTrain2.trainlen));
%
sM_train = som_seqtrain(sM_train_rough, sD_init, sTopol, sTrain2, ...
    'trainlen_type','samples', ... % samples/epochs
    'sample_order','random', ... % random/ordered
    'tracking',3); toc;
% ������ ����������:
[fqe, fte] = som_quality(sM_train, sD_init); 
fprintf('FQE = %f\nFTE = %f\n',fqe,fte);
%adm = som_distortion(sM_train,sD); fprintf('ADM = %f\n',adm);
%tp = somq_topographic_product(sM_train); fprintf('TP = %f\n',tp);

% ��������� ���������������� ��������� %
figure; clf; %hold on
plot(sD_init.data(:,1),sD_init.data(:,2),'k.'); %pause
% �������-����������:
bmus_indxs_indata = som_bmus(sM_train, sD_init, 'best'); % ������� BMU � ������������ ������
%bmus_data = D_init_all(bmus_indxs,:); % [Dlen,Ddim] - ???
%data_bmus_indxs=horzcat(sD_init.data,bmus_indxs_indata); % ������ + ������ �����. BMU
bmus_indxs = unique(bmus_indxs_indata); % ������� BMUs
bmu_num_best = length(bmus_indxs);
bmus_vecs = sM_train.codebook(bmus_indxs,:); % ������� BMU
bmus_vecs_indata = sM_train.codebook(bmus_indxs_indata,:); % ������� BMU � ������������ ������
plot(bmus_vecs(:,1),bmus_vecs(:,2),'bo'); %pause
bmusAll_indxs_indata = som_bmus(sM_train, sD_init, 'all');
bmusAll_indxs = unique(bmusAll_indxs_indata);
bmu_num_all = length(bmusAll_indxs);
bmusWorst_indxs_indata = som_bmus(sM_train, sD_init, 'worst');
bmusWorst_indxs = unique(bmusWorst_indxs_indata);
bmu_num_worst = length(bmusWorst_indxs);
fprintf('BMUall: %d\nBMUbest: %d\nBMUworst: %d\n',bmu_num_all,bmu_num_best,bmu_num_worst);
%% �������� ������� ����������:
[dmatmin_indxs_insom,dmatmin_indxs] = som_dmatclusters(sM_train,'neighf','N1');
% dmatmin_indxs - ������� ��������-������� ���������
% dmatmin_indxs_insom - ������� ��������-������� ��������� � ������������ SOM
dmatmin_vecs = sM_train.codebook(dmatmin_indxs,:); % ������� ��������-������� ���������
norms=[]; norm_min = []; norm_min_indxs = [];
for i=1:Dlen
    for j=1:numel(dmatmin_vecs(:,1))
        norms(j) = norm(sD_init.data(i,:)-dmatmin_vecs(j,:),1);
    end
    norm_min = min(norms);
    norm_min_indxs(i)=find(norms==norm_min); % ������� ��������-������� ��������� � ������������ ������
end
norms = norms'; norm_min = norm_min'; norm_min_indxs = norm_min_indxs';
dmatmin_vecs_indata = dmatmin_vecs(norm_min_indxs,:);
plot(dmatmin_vecs_indata(:,1),dmatmin_vecs_indata(:,2),'rx','MarkerSize',10);
hold off
%% Sammon of SOM
sM = sM_train; figure(5); clf;
%
tic; projM2d = sammon(sM, 2, 10^(-5), '3'); toc; % in 2D
%{
tic; projM3d = sammon(sM, 3, 10^(-9), '2'); toc; % in 3D
%}
%% ����������� SOM ������ ����������� ������� %%
%sM = sM_init;
sM = sM_train;
%sM = projM2d;
%sD = sD2d;
%sD = sD3d;

figure(6); clf; %title('Visualization of map projection'); 
disp('Visualization of map projection');

som_grid(sM,'Coord', sM.codebook, 'Linecolor','r'); hold on; 
scatter(sD.data(:,1),sD.data(:,2),'b.'); %for 2D
%scatter3(sD.data(:,1),sD.data(:,2),sD.data(:,3),'b.'); %for 3D
%% ������������ 2D %%
sM = sM_train; % init/train
sD = sD_init;
%sD = sD2d;
%sD = sD3d;
figure(7); clf; %title('Fig.4');
%colormap(flipud(gray));
colormap(gray);
som_show(sM, 'comp','all', 'edge','off');
disp('Pic.4.1 Component planes'); pause;
subplot(2,1,1)
som_show(sM, 'umat','all'); 
disp('Pic.4.2.1 U-matrix'); pause;
%{
som_show_add('traj', Bmus);
disp('Pic.4.2.2 U-matrix + BMUs'); pause;
%}
h1 = som_hits(sM, sD.data(1:Dlen/2,:));      % 100 red
h2 = som_hits(sM, sD.data(Dlen/2+1:Dlen,:)); % 010 green
%h3 = som_hits(sM, sD.data(301:Dlen,:)); % 001 blue
%h4 = som_hits(sM, sD.data(762:764,:)); % 011 cyan
%h5 = som_hits(sM, sD.data(765:767,:)); % 101 magenta
%h6 = som_hits(sM, sD.data(768:Dlen,:)); % 110 yellow

som_show_add('hit', [h1,h2], ...
    'MarkerColor', [1 0 0;0 1 0], 'Subplot', 1);
disp('Pic.4.2.3 U-matrix + Sample hits histogram');

%% ���������� ���������� ������� SOM %%
%clf;
%som_grid(sM_init, 'Linecolor', 'k'); %view(0, -90);
%title('Map grid'); disp('SOM grid');

%% ������������ 3D %%
sM = sM_train; 
sD = sD_init;
%sD = sD2d;
%sD = sD3d;
figure(8); clf; grid on; colormap(gray);
%subplot(2,1,1)
Co = som_unit_coords(sM);
U = som_umat(sM);
U = U(1:2:size(U,1), 1:2:size(U,2));
som_grid(sM, 'Coord', [Co, U(:)], 'Surf', U(:),...
    'Marker','.', 'MarkerSize',4, 'MarkerColor','r',...
    'LineWidth',0.3, 'LineColor','k');
view(-80,45), axis tight, title('Distance matrix')
disp('3D surface of the Distance Matrix');
%{
subplot(2,1,2)
%plot(sD3d.data(:,1), 'r+'), hold on
som_grid(sM, 'Coord', [Co, U(:)], 'Linecolor', 'k');
view(-80,45), axis tight, title('Map grid');
disp('Pic.4.2 SOM grid surface'); 
%}
%% Clustering SOM
[c, p, err, ind] = kmeans_clusters(sM_train, 2, 5); % find clusterings: round(sqrt(mUnits))
[dummy,i] = min(ind); % select the one with smallest index
som_show(sM_train,'color',{p{i},sprintf('%d clusters',i)}); % visualize
colormap(jet(i)), som_recolorbar % change colormap
h1 = som_hits(sM_train, sD_init.data(1:4160,:));
h2 = som_hits(sM_train, sD_init.data(4161:Dlen,:));
%h3 = som_hits(sM_train, sD.data(301:Dlen,:));
som_show_add('hit', [h1,h2], ...
    'MarkerColor', [0.7 0 0; 0 0.7 0], 'Subplot', 1);
%% ������������� SOM
%X = D_init_all; % ������������� �������� ������
%X = sD_init.data; % ������������� �������� ������
%X = dmatmin_vecs_indata; % ������������� �� ��������� ������� ���������� SOM
X = bmus_vecs_indata; % ������������� �� ����� ���������� ������������ SOM
%X = sM_train.codebook; % ������������� �� ���� �������� SOM
numClasses = 2; % ����� �������������� �������
ClSmpls = Dlen/2; % ����� ��������� � ������
%% ���������� ������� ������� k-�������
opts = statset('Display','final');
[cidx,C] = kmeans(X,numClasses,'Distance','cityblock',...
    'Start','plus','Replicates',10,'Options',opts);
Xlabled = horzcat(X,cidx);
%{
Xlabled1 = sortrows((Xlabled(1:4160,:)),22);
Xlabled2 = sortrows((Xlabled(4161:Dlen,:)),22);
%}
figure; clf; grid on; hold on
plot(X(cidx==1,1),X(cidx==1,2),'r.','MarkerSize',12)
plot(X(cidx==2,1),X(cidx==2,2),'g.','MarkerSize',12)
plot(X(cidx==3,1),X(cidx==3,2),'b.','MarkerSize',12)
plot(X(cidx==4,1),X(cidx==4,2),'c.','MarkerSize',12)
plot(X(cidx==5,1),X(cidx==5,2),'m.','MarkerSize',12)
plot(X(cidx==6,1),X(cidx==6,2),'y.','MarkerSize',12)

plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3) 
title 'Cluster Assignments and Centroids'
%plot(D_init_all(:,1),D_init_all(:,2),'k.');
hold off

%sil = silhouette(X,cidx,'cityblock');
%figure; silhouette(X,cidx);
%% 
% ���������� ������� ������� SINGLE-LINKAGE
Z = linkage(X,'single'); %,'cityblock'
cutoff = median([Z(end,3)]);
cidx = cluster(Z,'maxclust',numClasses); % 'cutoff',cutoff
dendrogram(Z,'ColorThreshold',cutoff);
%dendrogram(Z);
Xlabled = horzcat(X,cidx);

%% �������� ������������ ������������� ����� �������� ���������
fprintf('---\n');
f = @(a) struct('value', num2cell(unique(a)), 'times', arrayfun(@(b) {sum(a==b)}, unique(a)));
arrayfun(@(a)fprintf('����� %g ����������� %d ���\n', a.value, a.times),...
    f(cidx(1:ClSmpls)));
f = @(a) struct('value', num2cell(unique(a)), 'times', arrayfun(@(b) {sum(a==b)}, unique(a)));
arrayfun(@(a)fprintf('����� %g ����������� %d ���\n', a.value, a.times),...
    f(cidx(ClSmpls+1:Dlen)));

%% ������ �������� �������������
targetsVec = data(:,3)'; % ������� ����� �������
myTargets = zeros(numClasses,Dlen);
targetsIdx = sub2ind(size(myTargets), targetsVec, 1:Dlen);
myTargets(targetsIdx) = 1;
outputsVec = cidx'; % ����������������� ����� �������
myOutputs = zeros(numClasses,Dlen);
outputsIdx = sub2ind(size(myOutputs), outputsVec, 1:Dlen);
plotconfusion(myTargets,myOutputs); % ������� ��������������
%% ������������ ������ ��� BMU
DT = delaunayTriangulation(bmus_vecs);
clf; triplot(DT); hold on; grid on;
plot(bmus_vecs(:,1),bmus_vecs(:,2),'r.','MarkerSize',4);
