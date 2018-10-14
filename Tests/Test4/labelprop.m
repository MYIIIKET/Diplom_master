function owner = labelprop(sM,slabel,sigma,disttype,nclass,iter)
    addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));
    if (nargin < 6) || isempty(iter),
        iter = 10000; 
    end
    
    if (nargin < 5) || isempty(nclass),
        nclass = max(slabel); 
    end
    
    if (nargin < 4) || isempty(disttype),
        disttype = 'euclidean';
    end
    
    qtnode = size(sM.codebook,1);  
    
%     
    d=som_eucdist2(sM,sM);
    W=-exp(-(d.^2)/sigma^2);
%     W = exp(-squareform(pdist(sM.codebook,disttype).^2)/2*sigma^2); 
    

    D = sparse(diag(sum(W,2)));
    
    Y = zeros(qtnode,nclass); 
    
    labelednodes = find(slabel);
    
    Y(sub2ind(size(Y),labelednodes,slabel(labelednodes))) = 1;
    
    noch = 0; 
    
    [~,owner] = max(Y,[],2);
    
    DW = D\W;
    
    clear D W;
    
    YI = Y;
    
    YTP = repmat(slabel~=0,1,nclass);
    
    YTN = 1-YTP;
    
    for j=1:iter  
        
        Y = DW * Y;
        Y = YTP .* YI + YTN .* Y;
        ownerbak = owner;
        [~,owner] = max(Y,[],2);
        
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