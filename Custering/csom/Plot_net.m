function  Plot_net( m,w )
% Строим сеть на графике
% Чтобы было быстрее, строим графики каждые 100 итераций
[w_rows, w_cols,dim] = size(w);
wNum = w_rows*w_cols;
    
 plot(m(1,:),m(2,:),'+g');
 
pts = reshape(w,wNum,length(m(:,1)));

wt = zeros(size(w,2),size(w,1),size(w,3));
for  ii=1:w_rows
    for  jj=1:w_cols
        wt(ii,jj,:) = w(jj,ii,:);
    end;
end;

ptst = reshape(wt,wNum,length(m(:,1)));

    for k=0:w_cols-1
        hold on;
        plot(pts(k*w_rows+1:(k+1)*w_rows,1), ...
            pts(k*w_rows+1:(k+1)*w_rows,2),'-r');
    end;    
 
    for k=0:w_rows-1
        plot(ptst(k*w_cols+1:(k+1)*w_cols,1), ...
            ptst(k*w_cols+1:(k+1)*w_cols,2),'.-r','MarkerSize',16);
    end;  
    
    hold off;
end