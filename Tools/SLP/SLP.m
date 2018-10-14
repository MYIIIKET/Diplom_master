function [pred] = SLP(W,y,l_ind,stepSize,T)
% Description
%
% SLP takes,
%   W           - affinity matrix, labels should be at the left-top corner,
%               - should be in sparse form
%   y           - label vector
%   l_ind       - labeled indexes
%   stepSize    - coefficient: step size
%   T           - coefficient: running epoches
%
% and returns,
%   pred        - prediction of f, in the original sort
% 
% [1] De-Ming Liang and Yu-Feng Li. Lightweight Label Propagation for 
%     Large-Scale Network Data. In: Proceedings of the 27th International 
%     Joint Conference on Artificial Intelligence (IJCAI'18).


%% default configuration
if nargin == 3
    stepSize = 0.01;
    T = 5;
end

N = size(W,1);
u_ind = 1:N;
u_ind(l_ind) = [];
ind = [l_ind,u_ind];
W1 = W(ind,ind); 
Y = y(l_ind);

fu = Sgd(W1,Y,stepSize,T);

pred = zeros(N,1);
pred(l_ind) = y(l_ind);
pred(u_ind) = fu; 

end


function fu = Sgd(W, fl, stepSize,T)
% minimize energy function using stochastic gradient descent
% W      - affinity matrix, labels should be on upper left corner
% fl     - labels, n x k label identify matrix
% stepSize 
% T      - iteration epoches

n = size(W, 1);
l = size(fl, 1);

f = zeros(n, 1);
f(1:l) = fl(1:l);

%% find2cells
[row, col] = find(W(l+1:n,:));
E = length(col);
lnonzeros = find2cells(row, col, E, n, l);

%% iteration
f = iteration(W, f, lnonzeros, stepSize, T, l);

%% Result
fu = f(l+1:n);
end

