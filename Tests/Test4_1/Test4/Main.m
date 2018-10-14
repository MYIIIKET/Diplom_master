clear all;
close all;
clc;
addpath(genpath('SOM-Toolbox-master'));

%% Chainlink data
load chainlink.mat chainlink
data1 = load('Fs44kHzT250ms_250Hz.mat');
data1 = data1.data;
data2 = load('Fs44kHzT250ms_500Hz.mat');
data2 = data2.data;
data_struct = struct();
% data_struct.data{1}=data1;
% data_struct.data{2}=data2;
firs_chain = chainlink(1:500,:);
second_chain = chainlink(501:1000,:);
data_struct.data{1}=firs_chain;
data_struct.data{2}=second_chain;

%% Init data
sz = 100;
D1 = rand(sz, 3);
D2 = rand(sz, 3)/2.5+1;
data = [data1';data2'];
len = length(data(:,1));

%% Map init and train
sD = som_data_struct(chainlink);
sD = som_label(sD,'add',[1:1]','c1');
sD = som_label(sD,'add',[len/2+1:len/2+1]','c2');
type = 2;
if type == 1
    sM = som_make(sD,...
        'tracking', 3,...
        'algorithm', 'seq',...
        'mapsize', 'big',...
        'lattice', 'hexa',...
        'shape','sheet',...
        'neigh','gaussian',...
        'training', 'long',...
        'con_thresh',5,...
        'con_rest',20,...
        'delta',1,...
        'con');
elseif type == 2
    %%
    t = 1000;
    l = 50;
    o = ones(l,l);
    r = (1-(1:t)/t)/5;
    sM = som_randinit(sD);
    sM = som_seqtrain(sM,sD,...
        'tracking', 3,...
        'trainlen',l,...
        'con_thresh',1,...
        'con_rest',2,...
        'delta',1,...
        'con');
    figure;
    som_grid(sM,'Coord',sM.codebook)
%     for i=1:t,
%         sM = som_seqtrain(sM,sD,'tracking',0,...
%             'trainlen',l,'samples',...
%             'alpha',0.4*o,'radius',(4*r(i)+1)*o,...
%             'con_thresh',1,...
%             'con_rest',5,...
%             'delta',1,...
%             'con');
%         som_grid(sM,'Coord',sM.codebook)
%         hold on, plot(sD.data(:,1),sD.data(:,2),'.'), hold off
%         if mod(i,100)==0
%           [q t] = som_quality(sM,sD)  
%         end
%         drawnow
%     end
end
%% 
sM = som_autolabel(sM,sD);
%% Optimal sigma (?)
% delta = 10e-8;
% U = som_umat(sM);
% sigma = getSigma(U,delta,100000);

%%
clc;
% clear all;
close all;
% load sM sM
% load sD sD

delta = 10e-9;
U = som_umat(sM);
sigma = getSigma(U,delta,100000);

d=som_eucdist2(sM,sM);

sM = som_label(sM,'clear','all'); 

sM = labelByBMU(sM,data_struct);

cluster_number = length(data_struct.data);
labels = convert2numeric(sM,cluster_number);

labels = labelprop(sM,labels,sigma);
sM = labelByNumeric(sM, labels);

sD = som_label(sD,'clear','all');
sD_test = som_label(sD,'clear','all'); 
sD_test.data = awgn(sD_test.data,0);
figure;
% plot3(sD_test.data(:,1),sD_test.data(:,2),sD_test.data(:,3),'.');
sD_test = som_autolabel(sD_test,sM);


sD_train = labelByOriginalData(sD,data_struct,cluster_number);

test = som_label2num(sD_test.labels);
train = som_label2num(sD_train.labels);
confusion.getMatrix(test, train, true)


Co = som_unit_coords(sM);

U = U(1:2:size(U,1), 1:2:size(U,2));
figure;
som_grid(sM, 'Coord', [Co, U(:)], 'Surf', U(:),...
    'Marker','.', 'MarkerSize',4, 'MarkerColor','r',...
    'LineWidth',0.3, 'LineColor','k');
% U = som_umat(sM);
% W=-exp(-(U.^2)/sigma^2);
% W = W(1:2:size(W,1), 1:2:size(W,2));
% figure;
% som_grid(sM, 'Coord', [Co, W(:)], 'Surf', W(:),...
%     'Marker','.', 'MarkerSize',4, 'MarkerColor','r',...
%     'LineWidth',0.3, 'LineColor','k');

figure;
som_show(sM,'umat','all');
labels1=labels;
labels2=labels;
labels1(labels(:)==2)=0;
labels2(labels(:)==1)=0;
figure;
som_show(sM,'umat','all');
som_show_add('hit',labels1,'MarkerColor','r','Subplot',1);
som_show_add('hit',labels2,'MarkerColor','g','Subplot',1);