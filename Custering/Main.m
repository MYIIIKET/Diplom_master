clear all;
close all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));
addpath(genpath('C:\Users\MYlll\Desktop\Features\flatData'));
addpath(genpath('C:\Users\MYlll\Desktop\Features\out'));

DATA = '_data'; MFCC = '_mfcc'; LPC = '_lpc'; PLP = '_plp';
ORIG = '_orig'; PINK = '_pink'; WHITE = '_white';
TRAIN = 'train'; TEST = 'test';
FSP = 'fsp'; MSP = 'msp';
PLUS_5_DB = '_add_5dB'; MINUS_5_DB = '_sub_5dB'; ZERO_DB = '_0dB';

%% 
load m1 m1;
load f1 f1;
load m2 m2;
load f2 f2;

%% Load train data 1
trainDataStruct1 = loadDataStruct([1], [1], TRAIN, MFCC, ORIG,...
    MSP, MINUS_5_DB, 'dictor');
% f1(:,1:round(end/2))
trainDataStruct1.flat_data = m1(:,1:end);

%% Load train data 2
trainDataStruct2 = loadDataStruct([5], [1], TRAIN, MFCC, ORIG,...
    MSP, MINUS_5_DB, 'dictor');
trainDataStruct2.flat_data = m2(:,1:end);

%% Load train data 3
trainDataStruct3 = loadDataStruct([5], [1], TRAIN, MFCC, ORIG,...
    MSP, MINUS_5_DB, 'dictor');
trainDataStruct3.flat_data = f3(:,1:end);
%% Unite train data
D = [trainDataStruct1.flat_data'; trainDataStruct2.flat_data']';
%% Projection
% tic
% dim = 2;
% proj = getProjection('sammon', D, dim);
% proj = D;
% toc

%% Prepare data
% sD = som_data_struct(D');
% sD = som_normalize(sD, 'var');
% sD.labels = [trainDataStruct1.som_data.labels; trainDataStruct2.som_data.labels];

%% Train map
% sM = som_randinit(D');
% while true
%     sM = som_seqtrain(sM,D','tracking',3,...
%         'trainlen',100,'samples');
% sM = som_make(sD,...,
%     'init', 'randinit',...
%     'algorithm', 'seq',...
%     'lattice', 'hexa',...
%     'shape', 'sheet',...
%     'neigh', 'gaussian',...
%     'training', 'long',...
%     'mapsize', 'big',...
%     'tracking', 3);
% end
%% 
D = D';
%% 
D = som_normalize(D, 'range');
%% 
sM = som_randinit(D);

%%
a = length(trainDataStruct1.flat_data);
b = a+length(trainDataStruct2.flat_data);

k=100;
trainlen = 1000;
alpha = 0.1;
o = ones(k,1);
r = (1-(1:k)/k);
% r = 10*max(D(:));
% r = 10;
dr = 0.05;
da = 0.00;
i=1;
[q_prev, t] = som_quality(sM,D);
q_cur = 1000;
%% 
% while true
    for i=1:k
        sM = som_seqtrain(sM,D,'tracking',0,...
            'trainlen',trainlen,'samples',...
            'alpha',alpha*o,'radius',(4*r(i)+1)*o,...
            'neigh','gaussian', 'mapsize', 'big');
        
        if mod(i,1)==0
            som_grid(sM,'Coord',sM.codebook(:,[1 2]));
            hold on;
            plot(D(1:a,1),D(1:a,2),'.r');
            plot(D(a+1:b,1),D(a+1:b,2),'.b');
            hold off;
            drawnow
            tmp = q_cur;
%             [q_cur, t] = som_quality(sM,D);
            fprintf('QE = %f; TE = %f; iter = %f; dQE = %f\n', q_cur, t, i, abs(q_cur-q_prev));
            q_prev = tmp;
%             if alpha>0.1
%                 alpha = alpha - alpha*da;
%             end
        end
%         i=i+1;
%         r = r-r*dr;
%         alpha = alpha - alpha*da;
    end
    [q_cur, t] = som_quality(sM,D)
% end
%% Add hits 
[q_cur, t] = som_quality(sM,D)
% sD = loadHits(sM, sD);

%% Visualization
sD = som_data_struct(D);
visualize(sM, sD)

%%
figure;
hold on;
a = length(trainDataStruct1.flat_data);
b = a+length(trainDataStruct2.flat_data);
c = b+length(trainDataStruct3.flat_data);
% som_grid(sM,'Coord',sM.codebook(:,[1 2]));
plot(D(1:a,1),D(1:a,2),'.r');
plot(D(a+1:b,1),D(a+1:b,2),'.b');
plot(D(b+1:c,1),D(b+1:c,2),'.g');
hold off;

%% Cluster structure
figure;
Z  = linkage(pdist(sM.codebook));
sC = som_clstruct(Z);
som_clplot(sC,sM);

%% kmean
figure;
[color,b]=som_kmeanscolor(sM,10);
figure;
som_show(sM,'color',{color(:,:,b),strcat('"Best clustering - "', num2str(b), ' Clusters')});
figure;
som_show(sM,'color',{color(:,:,2),'"3 Clusters"'});

% [c, p, err, ind] = kmeans_clusters(sM);
% [dummy,i] = min(ind);
% figure;
% som_show(sM,'color',{p{i},sprintf('%d clusters',i)});
% colormap(jet(i)), som_recolorbar

%% cPlane
figure;
for i=1:length(sM.codebook(1,:))
    figure;
    som_cplane(sM,'k',sM.codebook(:,i));
end


%% Dend
figure;
Z = som_linkage(sM); 
som_dendrogram(Z,sM);

%%
[cPCAarg, Pdata, Pproto] = som_projections(D,sM);
figure;
som_projections_plot('scatter',Pproto(:,1:3),Pproto(:,5:7),5,sM)



%% distances
% f1n = som_normalize(f1','range')';
pos = [D(:,1)'; D(:,2)']';
% A = sum(pos.^2, 2);
% B = sum(pos.^2, 2)';
% AB = pos*pos';
% distances = sqrt(bsxfun(@plus, B, bsxfun(@minus, A, 2*AB)));

%% 
f1n = som_normalize(f1','range')'; pos1 = [f1n(1,:); f1n(2,:)]';
f2n = som_normalize(f2','range')'; pos2 = [f2n(1,:); f2n(2,:)]';
f3n = som_normalize(f3','range')'; pos3 = [f3n(1,:); f3n(2,:)]';

grid = 128;
% minvals = min(pos1);
% maxvals = max(pos1);
% rangevals = maxvals - minvals;
% xidx = 1 + round((pos1(:,1) - minvals(1)) ./ rangevals(1) * (grid-1));
% yidx = 1 + round((pos1(:,2) - minvals(2)) ./ rangevals(2) * (grid-1));
% density = accumarray([yidx, xidx], 1, [grid,grid]);  %note y is rows, x is cols
% figure;
% imagesc(density, 'xdata', [minvals(1), maxvals(1)], 'ydata', [minvals(2), maxvals(2)]);

minvals = min(pos2);
maxvals = max(pos2);
rangevals = maxvals - minvals;
xidx = 1 + round((pos2(:,1) - minvals(1)) ./ rangevals(1) * (grid-1));
yidx = 1 + round((pos2(:,2) - minvals(2)) ./ rangevals(2) * (grid-1));
density = accumarray([yidx, xidx], 1, [grid,grid]);  %note y is rows, x is cols
figure;
imagesc(density, 'xdata', [minvals(1), maxvals(1)], 'ydata', [minvals(2), maxvals(2)]);


% minvals = min(pos3);
% maxvals = max(pos3);
% rangevals = maxvals - minvals;
% xidx = 1 + round((pos3(:,1) - minvals(1)) ./ rangevals(1) * (grid-1));
% yidx = 1 + round((pos3(:,2) - minvals(2)) ./ rangevals(2) * (grid-1));
% density = accumarray([yidx, xidx], 1, [grid,grid]);  %note y is rows, x is cols
% figure;
% imagesc(density, 'xdata', [minvals(1), maxvals(1)], 'ydata', [minvals(2), maxvals(2)]);


%%
pos= [f1(:,1)'; D(:,2)']';

grid = 256;
minvals = min(D);
maxvals = max(D);
rangevals = maxvals - minvals;
xidx = 1 + round((D(:,1) - minvals(1)) ./ rangevals(1) * (grid-1));
yidx = 1 + round((D(:,2) - minvals(2)) ./ rangevals(2) * (grid-1));
density = accumarray([yidx, xidx], 1, [grid,grid]);  %note y is rows, x is cols
figure;
imagesc(density, 'xdata', [minvals(1), maxvals(1)], 'ydata', [minvals(2), maxvals(2)]);
colorbar

%%
hist(pos3)
