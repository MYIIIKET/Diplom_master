function [index]=minmatrix(D,Ci,top)
% great !!!!!!已经经过了完善的测试,效果完全令人满意.
% D是输入的两块流形上点之间的距离矩阵. Ci是不同块之间所需要连接的边数 
% index=(x_1,y_1;x_2,y_2,......;x_Ci,y_Ci)^T  分别是所需要矩阵的元素位置
% top=1;   %每个节点所连接其它流形上点的最高个数
[row,col]=size(D);
L=reshape(D,row*col,1);
[sortL,indexL]=sort(L);  %reshape总是一列用完再用下一列.
index=[];iter=0; i=1; timesI=0;timesJ=0;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

while(iter<Ci)&(i<=row*col) %没有加入对i上界的限制,一般来说不会出问题的.首先求出当前边indeL(i)的位置(indeI,indeJ)
    temp=indexL(i);
    remind=rem(temp,row);                                    %行数
    exactdivide=floor(temp/row);                             %列数
    if(remind==0)
        indexI=row; 
        indexJ=exactdivide;
    else
        indexI=remind;
        indexJ=exactdivide+1;
    end
%判断每个节点连接其它节点的个数,给出上界top必须满足 
if ~isempty(index)
 block1=index(:,1);
 block2=index(:,2);
 timesI=sum(block1==indexI);
 timesJ=sum(block2==indexJ);
end  %if ~isempty(index)
%-------------------------------------------------------------- 
if (timesI<top)&&(timesJ<top)
index((iter+1),:)=[indexI indexJ];
iter=iter+1;                %更新发生一次.
end  % if(timesI<top)&&(timesJ<top) 
%假如满足条件,就更新index.否则就继续试下一个距离.
i=i+1;timesI=0;timesJ=0;
end % while(iter<Ci)
 



