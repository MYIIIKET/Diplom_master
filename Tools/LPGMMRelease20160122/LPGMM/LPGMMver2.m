function [acc, ydata] = LPGMMver2(Xtag, labeltag, Xtest, labeltest,opt)

Xtag = Xtag';   Xtest = Xtest';    X = [Xtag; Xtest];

labeltag = labeltag(:);  labeltest = labeltest(:);  


Classlabels = unique(labeltag);  Classnum = length(Classlabels);

% Initialize some variables
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

% ---------- Argument defaults Ä¬ÈÏ²ÎÊý ----------
if ~exist('opt','var'),                              opt = [];      end;
if ~isfield(opt,'perplexity'),            opt.perplexity = 30;      end;
if ~isfield(opt,'epsilon'),                  opt.epsilon = 300;     end  % learning rate
if ~isfield(opt,'min_gain'),                opt.min_gain = 0.01;    end  % initial jitter
if ~isfield(opt,'momentum'),                opt.momentum = 0.5;     end  % initial momentum
if ~isfield(opt,'final_momentum'),    opt.final_momentum = 0.8;     end  % final momentum
if ~isfield(opt,'max_iter'),                opt.max_iter = 300;     end  % maximum number of iterations
if ~isfield(opt,'stop_lying_iter'),  opt.stop_lying_iter = 100;     end  % jitter decay
if ~isfield(opt,'mom_switch_iter'),  opt.mom_switch_iter = 250;     end  % iteration where momentum changes
if ~isfield(opt,'v'),                             opt.v = 1;        end  % iteration where momentum changes
% ---------- End of "argument defaults" ----------

% Perform preprocessing using PCA
if isfield(opt,'ReducedDim')
    disp('Preprocessing data using PCA...');
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
% Perform preprocessing using PCA


model = consGraphVectors(X', opt.n_size);
% Compute joint probabilities
P = d2p(model.G, opt.perplexity, 1e-5);                                           % compute affinities using fixed perplexity
clear model;
P(1:n + 1:end) = 0;   P = (P+P')/2;   P = max(P ./ sum(P(:)), realmin); 

const = sum(P(:) .* log(P(:)));                     % constant in KL divergence
P = P * 4;                                          % lie about the P-vals to find better local minima
y_incs  = zeros(size(ydata));   gains = ones(size(ydata));

% Run the iterations
momentum = opt.momentum;

for iter=1:opt.max_iter
    sum_ydata =  sum(ydata .^ 2, 2);
    num = 1 + (bsxfun(@plus, sum_ydata, bsxfun(@plus, sum_ydata', -2 * ydata * ydata')) ./ opt.v);
    Q = num .^ -((opt.v + 1) / 2);                                                  % apply power we 'forgot'
    Q(1:n+1:end) = 0;                                                           % set diagonal to zero
    Q = Q ./ sum(Q(:));                                                         % normalize to get probabilities
    Q = max(Q, eps);                                                            % make sure we don't have zero probs
    num = 1 ./ num;                                                             % this term is needed in the gradient
    num(1:n+1:end) = 0;                                                         % set diagonal to zero
    
    % Compute the gradients (faster implementation) 
    y_grads = zeros(size(ydata));
    stiffnesses = (4 .* ((opt.v + 1) ./ (2 .*opt.v))) .* (P - Q) .* num;
    for theta= Ntrain + 1 :n
        y_grads(theta,:) = sum(bsxfun(@times, bsxfun(@minus, ydata(theta,:), ydata), stiffnesses(:,theta)), 1);
    end 

    % Update the solution
    gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...         % note that the y_grads are actually -y_grads
        + (gains * .8) .* (sign(y_grads) == sign(y_incs));
    
    y_incs = momentum * y_incs - opt.epsilon * (gains .* y_grads);
    ydata = ydata + y_incs;
    
    ydata = bsxfun(@minus, ydata, mean(ydata, 1));
    ydata(1:Ntrain,:) = ydata(1:Ntrain,:) + ytag; 
    
      
    % Update the momentum if necessary
    if iter == opt.mom_switch_iter
       momentum = opt.final_momentum;
    end
    if iter == opt.stop_lying_iter
        P = P ./ 4;
    end
    
    % Print out progress
    if ~rem(iter, 50)
        % Compute the value of the cost function
        cost = sum(sum(P .* log((P + eps) ./ (Q + eps))));
        disp(['Iteration ' num2str(iter) ': cost  = ' num2str(cost)]);
    end
    
    % Display scatter plot (maximally first three dimensions)
    if ~rem(iter, 50)  
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

    