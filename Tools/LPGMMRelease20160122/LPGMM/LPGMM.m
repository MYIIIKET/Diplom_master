function [acc, ydata] = LPGMM(Xtag, labeltag, Xtest, labeltest,opt)

%%%  test :  optimize the objective function
Xtag = Xtag';   Xtest = Xtest';    X = [Xtag; Xtest];

labeltag = labeltag(:);  labeltest = labeltest(:);  

Classlabels = unique(labeltag);  Classnum = length(Classlabels);
 
n = size(X, 1);   % number of instances
Ntrain = length(labeltag);

% Initialize embedding coordinates randomly (close to origin)
R = 0.1;% R = 0.01;
ytag = -R*ones(length(labeltag),Classnum);    
for iter1=1:Classnum,
    ind1 = labeltag==Classlabels(iter1);
    ytag(ind1,iter1)=R;
end
ytest = 0.0001 * randn(length(labeltest), Classnum);  
ydata = [ytag; ytest];
opt = optsetfun(opt);
if ~isfield(opt,'ReducedDim'),           
    opt.ReducedDim = 40;       
end  % iteration where momentum changes
 
 
if isfield(opt,'ReducedDim')
    disp('Preprocessing data...');
    if size(X, 2) < size(X, 1)
        C = X' * X;
    else
        C = (1 / size(X, 1)) * (X * X');
    end
    [M, lambda] = eig(full(C));
    [lambda, ind] = sort(diag(lambda), 'descend');
    M = M(:,ind(1:opt.ReducedDim));
    lambda = lambda(1:opt.ReducedDim);
    if ~(size(X, 2) < size(X, 1))
        M = bsxfun(@times, X' * M, (1 ./ sqrt(size(X, 1) .* lambda))');
    end
    X = bsxfun(@minus, X, mean(X, 1)) * M;
    clear M lambda ind
end 


model = consGraphVectors(X', opt.n_size);
 
P = d2p(model.G, opt.perplexity, 1e-5);                                          
clear model;
P(1:n + 1:end) = 0;   
P = (P+P')/2;   
P = max(P ./ sum(P(:)), realmin); 

const = sum(P(:) .* log(P(:)));                    
P = P * 4;                                          
y_incs  = zeros(size(ydata));   gains = ones(size(ydata));

% Run the iterations
momentum = opt.momentum;

for iter=1:opt.max_iter
    sum_ydata =  sum(ydata .^ 2, 2);
    num = 1 + (bsxfun(@plus, sum_ydata, bsxfun(@plus, sum_ydata', -2 * ydata * ydata')) ./ opt.v);
    Q = num .^ -((opt.v + 1) / 2);                                                  
    Q(1:n+1:end) = 0;                                                           
    Q = Q ./ sum(Q(:));                                                         
    Q = max(Q, eps);                                                            
    num = 1 ./ num;                                                           
    num(1:n+1:end) = 0;                                  
    y_grads = zeros(size(ydata));
    stiffnesses = (4 .* ((opt.v + 1) ./ (2 .*opt.v))) .* (P - Q) .* num;
    for theta= Ntrain + 1 :n
        y_grads(theta,:) = sum(bsxfun(@times, bsxfun(@minus, ydata(theta,:), ydata), stiffnesses(:,theta)), 1);
    end 
 
    gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...        
        + (gains * .8) .* (sign(y_grads) == sign(y_incs));
    
    y_incs = momentum * y_incs - opt.epsilon * (gains .* y_grads);
    ydata = ydata + y_incs;
    
    ydata = bsxfun(@minus, ydata, mean(ydata, 1));
    ydata(1:Ntrain,:) = ydata(1:Ntrain,:) + ytag; 
    
       
    if iter == opt.mom_switch_iter
       momentum = opt.final_momentum;
    end
    if iter == opt.stop_lying_iter
        P = P ./ 4;
    end
    
 
    if ~rem(iter, 50)
        disp(['Iteration ' num2str(iter)]);
        Sy = ydata(1:Ntrain,:);
        if Classnum == 1
            scatter(ydata, ydata, 9, [labeltag; labeltest], 'filled');
        elseif Classnum == 2
            plot(Sy(:,1), Sy(:,2), 'rd'); hold on;
            scatter(ydata(:,1), ydata(:,2), 9, [labeltag; labeltest], 'filled');
        else
            plot3(Sy(:,1), Sy(:,2), Sy(:,3),'rd'); hold on;
            scatter3(ydata(:,1), ydata(:,2), ydata(:,3), 9, [labeltag; labeltest], 'filled');
        end
        hold off;
        drawnow
    end
end

[~,inds] = max(ydata(length(labeltag)+1:end,:),[],2);
acc = 100*mean(inds(:) == labeltest(:));

    