%% process input data
data_path = ...
    'C:\Users\user_scc_106\Documents\MATLAB';
cd(data_path)

D_init_sammon = D_init_all;
D_init_sammon = [D_init_sammon zeros(size(D_init_sammon,1),1)];
[x,idx] = unique(D_init_sammon(:,1:13), 'rows');
t = D_init_sammon(idx,14); %idx,N-элементов вектора + 1 метка
n = size(x, 1);
%% modify options and perform Sammon mapping 
opts                = sammon;
opts.Display        = 'iter';
opts.TolFun         = 1e-9;
opts.Initialisation = 'pca';

tic
for i=1:1 % 1:N-эпохи
[y, E] = sammon(x, 2, opts);
end
toc
%% plot results
figure(1); %clf;
set(axes, 'FontSize', 9);
plot(y(t == 0,1), y(t == 0,2), 'r.');
title(['Sammon Mapping (stress = ' num2str(E) ')']);
%legend('dtmf 0', 'dtmf 1', 'dtmf 2', 'dtmf 3', 'dtmf 4', 'dtmf 5');