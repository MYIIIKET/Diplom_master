addpath('utils');

dataname = 'covtype';
data = load(['data/',dataname],'W','y'); % load data
W0 = data.W;
y = data.y;

T = 20;
times = zeros(T,1);
accs = zeros(T,1); 
aucs = zeros(T,1);

tic;    
idx = load(['data_ssl/',dataname],'L','U'); % load labeled indexes
Ls = idx.L;
Us = idx.U;

for t = 1:T
    L = Ls{t};
    U = Us{t};
    time = toc;
    fs = SLP(W0,y,L,0.1,6);
    times(t) = toc - time;
    accs(t) = Accuracy(fs,y,U);
    [~,~,~,aucs(t)] = perfcurve(y,fs,'1');
    fprintf('%d: acc: %.3f, auc: %.3f, time: %.3f\n',t,accs(t),aucs(t),times(t));    
end

fprintf('MEAN\nacc: %.3f, auc: %.3f, time: %.3f\n',mean(accs),mean(aucs),mean(times)); 

