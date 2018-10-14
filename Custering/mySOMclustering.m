clear all; 
%clc; 
echo off
disp('_______________________________________________________________');
%% Загрузка и подготовка обучающей выборки, настройка SOM:
data_path_train = ...
    'D:\myworkspace\data\speech\workingproc';
cd(data_path_train)
%
% voice1:
%load('fsp1train1_mfcc.mat');
load('fsp1train2_mfcc.mat');
load('fsp1train3_mfcc.mat');
load('fsp1train4_mfcc.mat');
load('fsp1train5_mfcc.mat');
load('fsp1train6_mfcc.mat');
load('fsp1train7_mfcc.mat');
%
D_init_all = horzcat(...
    fsp1train2_mfcc_orig, fsp1train3_mfcc_orig, ...
    fsp1train4_mfcc_orig, fsp1train5_mfcc_orig, ...
    fsp1train6_mfcc_orig, fsp1train7_mfcc_orig);
%{
% voice2:
%load('msp1train1_mfcc.mat');
load('msp1train2_mfcc.mat');
load('msp1train3_mfcc.mat');
load('msp1train4_mfcc.mat');
load('msp1train5_mfcc.mat');
load('msp1train6_mfcc.mat');
load('msp1train7_mfcc.mat');
%{
D_init_all = horzcat(...
    msp1train2_mfcc_orig, msp1train3_mfcc_orig, ...
    msp1train4_mfcc_orig, msp1train5_mfcc_orig, ...
    msp1train6_mfcc_orig, msp1train7_mfcc_orig);
%}
D_init_all = horzcat(...
    fsp1train2_mfcc_orig, fsp1train3_mfcc_orig, ...
    msp1train6_mfcc_orig, msp1train7_mfcc_orig);
%}
%% Подготовка данных:
D_init_all = D_init_all';
% Создание структуры данных и нормализация:
sD_init = som_data_struct([D_init_all]);
sD_init = som_normalize(sD_init, 'var');
% Генерация и настройка SOM:
%{
sM_init = som_make(sD_init, ...
    'init', 'randinit', 'neigh', 'gaussian', 'algorithm', 'seq', ...
    'msize', [39 45], 'training', 'long', 'tracking', 2); %mapsize big
%}
%% Отображение Сэммона:
D_init_sammon = [sD_init.data zeros(size(sD_init.data,1),1)];
[x,idx] = unique(D_init_sammon(:,1:13), 'rows');
t = D_init_sammon(idx,13+1); %idx,N-элементов вектора + 1 метка
n = size(x, 1);

y = sammon(D_init_sammon,3);
%
figure(1); %clf;
set(axes, 'FontSize', 9);
plot(y(t == 0,1), y(t == 0,2), 'r.');
title(['Sammon Mapping']);

%% Настройка карты 1):
msizeY = 10; msizeX = 20; munits = msizeX*msizeY;
sTopol = som_set('som_topol', 'msize',[msizeX msizeY], ...
    'lattice','hexa', 'shape','sheet');
sM_init = som_randinit(sD_init, sTopol); %randinit/lininit

disp('> Rough training phase...');
sTrain1 = som_train_struct(sM_init, sTopol, 'rough');
sTrain1 = som_set(sTrain1, 'algorithm','seq', ...
    'neigh','gaussian', ...
    'radius_ini',20, 'radius_fin',1, ...
    'alpha_ini',0.5, ...
    'trainlen',2);
fprintf('training length: %d\n',sTrain1.trainlen);
sM_train = som_seqtrain(sM_init, sD_init, sTopol, sTrain1, ...
    'trainlen_type','epoches', ... % samples/epochs
    'sample_order','random', ...
    'tracking',3);
%[qe, te] = som_quality(sM_train, sD_init);
%fprintf('QE = %f\nTE = %f\n',qe,te);
%% Настройка карты 2):
disp('> Finetuning...');
sTrain2 = som_train_struct(sM_train, sTopol, 'finetune');
sTrain2 = som_set(sTrain2, 'algorithm','seq', ...
    'neigh','gaussian', ...
    'radius_ini',1, 'radius_fin',0.1, ...
    'alpha_ini',0.07, ...
    'trainlen',munits*length(D_init_all)); %500
fprintf('training length: %d\n',sTrain2.trainlen);
sM_train = som_seqtrain(sM_train, sD_init, sTopol, sTrain2, ...
    'trainlen_type','samples', ... % samples/epochs
    'sample_order','random', ...
    'tracking',3);
%% Оценка построения:
[fqe, fte] = som_quality(sM_train, sD_init);
adm = som_distortion(sM_train,sD_init);
fprintf('FQE = %f\nFTE = %f\nADM = %f\n',fqe,fte,adm)
%}
%% Загрузка и подготовка тестовой выборки, настройка SOM:
%{
data_path_test = ...
    'D:\myworkspace\data\speech\testset\workingproc'; 
cd(data_path_test)
% voice1:
D11_test = 'ANprob_0500Hz.mat'; load(D11_test);
D12_test = 'ANprob_4000Hz.mat'; load(D12_test);
% voice2:
D21_test = 'ANprob_0500Hz.mat'; load(D21_test);
D22_test = 'ANprob_4000Hz.mat'; load(D22_test);
%
D_test_all = horzcat(...
    ANprob_1000Hz, ANprob_2000Hz, ...
    ANprob_6000Hz, ANprob_8000Hz);
D_test_all = D_test_all';
%
sD_test = som_data_struct([D_test_all]);
sD_test = som_normalize(sD_test, 'var');
%
sM_test = som_seqtrain(sM_init, sD_test, 'tracking', 2);
%
[qe, te] = som_quality(sM_test, sD_test);
fprintf('QE = %f\nTE = %f\n',qe,te)
%}
%% Визуализация 2D после обучения:
sM2d = sM_train;
sD2d = sD_init;

colormap(gray);
som_show(sM2d, 'comp', ...
    [1 2 3 4 5 6 7 8 9 10 11 12 13], ...
    'edge', 'off'); % отображение компонент вектора, кол-во = размерности
disp('> Component planes'); fprintf('Press any key to continue\n'); pause
subplot(2,1,1)
som_show(sM2d, 'umat', 'all');
disp('> U-matrix'); fprintf('Press any key to continue\n'); pause
%{
Bmus = som_bmus(sM2d, sD2d, 'best');
som_show_add('traj', Bmus);
disp('> SOM BMUs'); fprintf('Press any key to continue\n'); pause
%}

h11 = som_hits(sM2d, sD2d.data(1:1509,:));
h12 = som_hits(sM2d, sD2d.data(1509:2647,:));
h13 = som_hits(sM2d, sD2d.data(2647:4519,:));
h14 = som_hits(sM2d, sD2d.data(4519:5347,:));
h15 = som_hits(sM2d, sD2d.data(5347:6245,:));
h16 = som_hits(sM2d, sD2d.data(6245:7063,:));
%{
h21 = som_hits(sM2d, sD2d.data(8061:9039,:));
h22 = som_hits(sM2d, sD2d.data(9039:10497,:));
h23 = som_hits(sM2d, sD2d.data(10497:11535,:));
h24 = som_hits(sM2d, sD2d.data(11535:12973,:));
h25 = som_hits(sM2d, sD2d.data(12973:13681,:));
h26 = som_hits(sM2d, sD2d.data(13681:14429,:));
%}
som_show_add('hit', [h11, h12, h13, h14, h15, h16], ...
    'MarkerColor', [0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0; 0 1 0], 'Subplot', 1);
disp('> Sample hits histogram'); fprintf('Press any key to continue\n'); pause

%% Визуализация 2D после тестирования:
%{
sM = sM_test;
sD = sD_test;

colormap(gray);
som_show(sM, 'comp', ...
    [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21], ...
    'edge', 'off'); 
disp('Component planes'); fprintf('\nPress any key to continue\n'); pause
subplot(2,1,1)
som_show(sM, 'umat', 'all');
disp('U-matrix'); fprintf('\nPress any key to continue\n'); pause
%{
Bmus = som_bmus(sM, sD, 'best');
som_show_add('traj', Bmus);
disp('SOM BMUs'); fprintf('\nPress any key to continue\n'); pause
%}

h11 = som_hits(sM, sD.data(1:22100,:));
h12 = som_hits(sM, sD.data(22101:44200,:));

h21 = som_hits(sM, sD.data(132601:154700,:));
h22 = som_hits(sM, sD.data(154701:176800,:));

som_show_add('hit', [h11, h12, h21, h22], ...
    'MarkerColor', [0 1 0; 0 1 0; 1 0 0; 1 0 0], 'Subplot', 1);
disp('Sample hits histogram'); fprintf('\nPress any key to continue\n'); pause
%}
%% Отобразить полученную решетку SOM:
clf
som_grid(sM_init, 'Linecolor', 'k');
%view(0, -90)
title('Map grid');
disp('SOM grid');
%% Визуализация 3D:
clf 
sM3d = sM_train;
sD3d = sD_init;

colormap(gray);
subplot(2,1,1)
Co = som_unit_coords(sM3d);
U = som_umat(sM3d);
U = U(1:2:size(U,1), 1:2:size(U,2));
som_grid(sM3d, 'Coord', [Co, U(:)], 'Surf', U(:), 'Marker', 'none');
view(-80,45), axis tight, title('Distance matrix')
disp('Surface of distance matrix'); fprintf('\nPress any key to continue'); pause

subplot(2,1,2)
%plot(sD3d.data(:,1), 'r+'), hold on
som_grid(sM3d, 'Coord', [Co, U(:)], 'Linecolor', 'k');
view(-80,45), axis tight, title('Map grid')
disp('SOM grid surface'); fprintf('\nPress any key to continue');

