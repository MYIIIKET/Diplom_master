function owner = zhou(X,slabel,sigma,disttype,alpha,nclass,iter)
    
    if (nargin < 7) || isempty(iter),
        iter = 10000; 
    end
    
    if (nargin < 6) || isempty(nclass),
        nclass = max(slabel); 
    end 
    
    if (nargin < 5) || isempty(alpha),
        alpha = 0.99;
    end
    
    if (nargin < 4) || isempty(disttype),
        disttype = 'euclidean';
    end
    
    qtnode = size(X,1);  
    
    W = exp(-squareform(pdist(X,disttype).^2)/2*sigma^2); 
    
    W = W - eye(qtnode);  
   
    D = diag(sum(W,2));
    
    DInv = sparse(D^(-1/2));
    
    clear D;
    
    S = DInv * W * DInv;
    
    clear W DInv;
    
    Y = zeros(qtnode,nclass); 
    
    labelednodes = find(slabel);
    
    Y(sub2ind(size(Y),labelednodes,slabel(labelednodes))) = 1;
    
    noch = 0; 
    
    F = Y;
    
    [~,owner] = max(double(F),[],2);
    
    alphaS = alpha * S;
    
    clear S;
    
    alphaY = (1 - alpha) * Y;
    
    clear Y;
   
    for i=1:iter 
        
        F = alphaS * F + alphaY;
        ownerbak = owner;
        [~,owner] = max(F,[],2);
        
        if sum(ownerbak~=owner)==0  
            noch = noch + 1;
            if noch>=100 
                break; 
            end
        else
            noch = 0;
        end
    
    end
   
end