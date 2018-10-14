clear;clc;

load('CBCL3000');     gY = label(:);   gX = NormalizeFea(test,0);
TtrainList = [200 440 680 920 1040 1280 1400];
optLPGMM = []; optLPGMM.n_size = 6;  optLPGMM.ReducedDim = 40; %%% CBCL3000 data


%% data preprocessing
[Dim,N] = size(gX);
Classlabels=unique(gY);
Classnum = length(Classlabels);


%% NtrainList = Range of Varying l/C
NtrainList = [1 2 4 6 8 10 12]; 
%% Parameters for comparing algorithms
n_size = 6;

Nround = 20;                    %%% Repeat the experiment Nround times

AccLPGMM = [];    AccGGMC = [];   AccLGC = [];    AccGRF = [];

for T = 1:Nround
    
    %% Random permutation of the data points
    rnperm = randperm(N);
    dataX = gX(:,rnperm);
    labelY = gY(rnperm);
    %% index of each class
    Dind=cell(Classnum,1);
    for iter1=1:Classnum,
        Dind{iter1}=find(labelY==Classlabels(iter1));
    end
    %%
    
    [Wlap] = DoubleWforMR(dataX, 'k', n_size);   %% Weight matrix
    
    for iter1 = 1:length(NtrainList)
        %% index of labeled training: ind1
        ind1 =[]; 
        %% indenx of unlabeled training, which are also the test point(for SSL): ind2
        ind2=[];
        for c=1:Classnum  
            Ntrain = NtrainList(iter1);
            ind1 = [ind1; Dind{c}(1:Ntrain)];              %%%  training index
            ind2 = [ind2; Dind{c}((1+Ntrain):end)];        %%%  test index
        end
        dataTrain = dataX(:,ind1);    labelTrain = labelY(ind1);      %%% labeled training data
        dataTest  = dataX(:,ind2);    labelTest  = labelY(ind2);      %%% unlabeled traing Also the test data

        %% ---------LPGMM Classification--------

        AccLPGMM(T,iter1) = LPGMM(dataTrain, labelTrain, dataTest, labelTest,optLPGMM)
        
        %%% -------Learning with Local and Global Consistency
        graph=TransductionModel(Wlap);
        graph.prior=ones(1,Classnum);
        [predict, error]=lgc(dataX',labelY,ind1,graph);
        AccLGC(T,iter1) = 100*(1-error)
        %%
        
        %%% -------Learning with Local and Global Consistency graph transductive learning
        [predict, error]=grf(dataX',labelY,ind1,Wlap);
        AccGRF(T,iter1) = 100*(1-error)
        %%%
        
        %%% -------Greedy Gradient based Max-Cut algorithm GGMC
        [predict, error]=ggmc(dataX',labelY,ind1,graph);
        AccGGMC(T,iter1) = 100*(1-error)
        %%% -------
    end
end

pause;

% iind =alphaList;
iind =NtrainList;

figure; hold on;

errorbar(iind,mean(AccLPGMM,1), 0.1*std(AccLPGMM),...
    'ro-','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','y','MarkerSize',8);%%%

errorbar(iind,mean(AccLGC,1), 0.1*std(AccLGC),...
    'm<-','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',8) %%% Traditional MR

errorbar(iind,mean(AccGRF,1), 0.1*std(AccGRF),...
    'bs-','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',8) %%% Traditional MR

errorbar(iind,mean(AccGGMC,1), 0.1*std(AccGGMC),...
    'b^-','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',8) %%% Traditional MR

xlabel('Number of Labeled points (Varying l/C)','fontsize',10); ylabel('Accuracy (%)','fontsize',10);
grid on;
legend('LPGMM', 'LGC','GRF','GGMC');
hold off








%% Parameters for comparing algorithms
 Ntrain = 2;
 
 n_size = 6;   

Nround = 10;                    %%% Repeat the experiment Nround times

AccLPGMM = [];    AccGGMC = [];   AccLGC = [];    AccGRF = [];


for T = 1:Nround
    
    %% data preprocessing
    [Dim,N2] = size(gX);
    Classlabels=unique(gY);   
    Classnum = length(Classlabels);            
    rnperm = randperm(N2);     dataX1 = gX(:,rnperm);       labelY1 = gY(rnperm);    
    
    for iter2 = 1:length(TtrainList)
        
        Tnum = TtrainList(iter2);      ind3 =[];
        for iter1=1:Classnum
            indC = find(labelY1==Classlabels(iter1));
            ind3 = [ind3; indC(1:Tnum)];                     %%% Total training index, where only 2 points are labeled !! Different from previous experiments
        end
        dataX = dataX1(:,ind3);  labelY = labelY1(ind3);%%% OverWrites these variables
        
        [Dim,N] = size(dataX);
        %%
        dataY = zeros(Classnum,N);   Dind=cell(Classnum,1);
        for iter1=1:Classnum,
            Dind{iter1}=find(labelY==Classlabels(iter1));                           
        end
        %%
        
        
        [Wlap] = DoubleWforMR(dataX, 'k', n_size);   %% Weight matrix
 
        
        ind1 =[]; ind2=[];
        for iter1=1:Classnum
            ind1 = [ind1; Dind{iter1}(1:Ntrain)];              %%% Only 2 points are labeled!!
            ind2 = [ind2; Dind{iter1}((1+Ntrain):end)];        %%% unlabeled data index
        end
        dataTrain = dataX(:,ind1);  labelTrain = labelY(ind1);    %%% 
        dataTest  = dataX(:,ind2);  labelTest  = labelY(ind2);    %%%
        
        %% ---------LPGMM Classification--------

        AccLPGMM(T,iter2) = LPGMM(dataTrain, labelTrain, dataTest, labelTest,optLPGMM)
        
        
        %% -------Learning with Local and Global Consistency
        graph=TransductionModel(Wlap);
        graph.prior=ones(1,Classnum);
        [predict, error]=lgc(dataX',labelY,ind1,graph);
        AccLGC(T,iter2) = 100*(1-error)
        %%
        
        %% -------Learning with Local and Global Consistency graph transductive learning
        [predict, error]=grf(dataX',labelY,ind1,Wlap);
        AccGRF(T,iter2) = 100*(1-error)
        %%
        
        %% -------ggmc
        [predict, error]=ggmc(dataX',labelY,ind1,graph);
        AccGGMC(T,iter2) = 100*(1-error)
        %% -------

    end
end


% iind =alphaList;
iind =TtrainList;

figure; hold on;

errorbar(iind,mean(AccLPGMM,1), 0.1*std(AccLPGMM),...
    'ro-','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','y','MarkerSize',8);%%%

errorbar(iind,mean(AccLGC,1), 0.1*std(AccLGC),...
    'm<-','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',8) %%%  

errorbar(iind,mean(AccGRF,1), 0.1*std(AccGRF),...
    'bs-','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',8) %%% 

errorbar(iind,mean(AccGGMC,1), 0.1*std(AccGGMC),...
    'b^-','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','g','MarkerSize',8) %%% 

xlabel('Number of Training points (Varying N/C)','fontsize',10); ylabel('Accuracy (%)','fontsize',10);
grid on;
legend('LPGMM', 'LGC','GRF','GGMC');
hold off








