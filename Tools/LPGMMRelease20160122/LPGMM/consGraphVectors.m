function model = consGraphVectors(X, n_size, Ci)
[~,N] = size(X); 
if (N<3000)
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
Dist = D; model.D = D;  model.n_size = n_size;

[~, ind] = sort(D); 
for i=1:N,     D(i,ind((2+n_size):end,i)) = 0;   end
D = sparse(D);   D = max(D,D');     model.De = D;    

landmarks = 1:N;    D = dijkstra(D, landmarks);     tcompsind = cell(0); 
[~, firsts] = min(D==inf);       %% first point each point connects to
[comps, I, J] = unique(firsts);    %% represent each connected component once
n_comps = length(comps)               %% number of connected components
size_comps = sum((repmat(firsts,n_comps,1)==(comps(1:n_comps)'*ones(1,N)))');
[~, comp_order] = sort(size_comps);  %% sort connected components by size
comps = comps(comp_order(end:-1:1));
for i = 1:n_comps
    tcompsind{i} = find(firsts==comps(i));    % indices of each components
end
G = D;
if n_comps>1
    if ~exist('Ci','var'),  Ci=4;     end
    for i = 1:n_comps-1
        for j = i+1:n_comps
            compsDist = Dist(tcompsind{i},tcompsind{j});
            indC=minmatrix(compsDist,Ci,1);
            for k = 1:size(indC,1)
                model.De(tcompsind{i}(indC(k,1)),tcompsind{j}(indC(k,2))) = compsDist(indC(k,1),indC(k,2));
                model.De(tcompsind{j}(indC(k,2)),tcompsind{i}(indC(k,1))) = compsDist(indC(k,1),indC(k,2));
            end    % for k = 1:Ci
        end    %for j = i+1:n_comps
    end    %for i = 1:n_comps-1
    D = model.De;  D = sparse(D);
    D = max(D,D');    %% Make sure distance matrix is symmetric
    G = dijkstra(D, 1:N);
    disp('reupdate distance matrix complete');
end   % if n_comps>1
clear dist;
model.G = G; model.X = X;
