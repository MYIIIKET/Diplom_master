function [Wlap Wlle] = DoubleWforMR(X, n_fcn, n_size,E)
% X is the original data set, whose columns are data points.
% n_fcn= 'k' or 'epsilon', n_size is a number.
%algo=3 for 'isomap' or 2 for 'lle' or 1 for 'lap'
% E is the 0-1 incidence matrix. If it exists, the neighborhood will
% be identified using this matrix.

%%%%%----Compute pair wise Euclidean distances
[Dim,N] = size(X);

if (N<3000)&(Dim<2000)
    D=L2_distance(X(:,1:N),X(:,1:N));
else
    D=[];
    for i=1:N
        for j=i:N
            D(i,j)=sqrt((X(:,i)-X(:,j))'*(X(:,i)-X(:,j)));
            D(j,i)=D(i,j);
        end
    end
end

dist=D;
%%%%%----Compute pair wise Euclidean distances Over!



%%%%% Laplace %%%%%%----Build the neighborhood graph G={X,W}
if ~exist('E','var')
    if n_fcn=='k'
        [tmp, ind] = sort(D);
        for i=1:N
            temp=ind((2+n_size):end,i);
            D(temp,i) = 0;
        end
        clear temp;
        D = sparse(D);
        D = max(D,D');    %% Make sure distance matrix is symmetric
    elseif n_fcn=='epsilon'
        warning off    %% Next line causes an unnecessary warning, so turn it off
        D =  D.*(D<=n_size);
        warning on
        D = sparse(D);
        D = max(D,D');    %% Make sure distance matrix is symmetric
    end
else
    D=E.*D;
end

De=D;
%%%---稀疏图，对于laplace方法，该稀疏图矩阵应该是对称的


A = D;
Wlap=sparse(N,N);
[A_i, A_j, A_v] = find(A);
for i = 1: length(A_i)
    if exist('sigma','var')
        Wlap(A_i(i), A_j(i)) = exp( -A(A_i(i), A_j(i))/(2*sigma^2) );
    else
        Wlap(A_i(i), A_j(i)) = 1;
    end
end;

%%%%%lle  %%%%%%----Build the neighborhood graph G={X,W}
D=dist;
if n_fcn=='k'
    [tmp, ind] = sort(D);
    for i=1:N
        temp=ind((2+n_size):end,i);
        D(temp,i) = 0;
    end
    clear temp;
    D = sparse(D);
    D = max(D,D');    %% Make sure distance matrix is symmetric
elseif n_fcn=='epsilon'
    warning off    %% Next line causes an unnecessary warning, so turn it off
    D =  D.*(D<=n_size);
    warning on
    D = sparse(D);
    D = max(D,D');    %% Make sure distance matrix is symmetric
end
De=D;  %%%---稀疏图

% NEIGHBORS
neighbors = cell(N);                      % extract neighbbours
for i=1:N
    neighbors{i} = find(D(:,i)~=0);
end
% RECONSTRUCTION WEIGHTS
tol=1e-5;                      % if tol isn't inputed, set tol to the default 1e-4
M = sparse(zeros(N,N));                                    % creat a zero matrix of size K by N i.e. initialize W
for i=1:N
    L = length(neighbors{i});
    z = X(:,neighbors{i})-repmat(X(:,i),1,L);      % compute the inner product of xi and its nbs
    C = z'*z;                                        % compute the Gram matrix G=C
    C = C + tol*trace(C)*eye(L)/L;                   % regularization
    invC = inv(C);                                   % compute the inverse of C
    M(neighbors{i},i) = sum(invC)'/sum(sum(invC));              % first sum columns then sum rows
    % sum(invC) is a row vector
end
Wlle = sparse(eye(N));                                        % initialize W
I=Wlle;
for i=1:N
    m = M(neighbors{i},i);
    Wlle(i,neighbors{i}) = Wlle(i,neighbors{i}) - m';
    Wlle(neighbors{i},i) = Wlle(neighbors{i},i) - m;
    Wlle(neighbors{i},neighbors{i}) = Wlle(neighbors{i},neighbors{i}) + m*m';
end;
Wlle=I-Wlle;

