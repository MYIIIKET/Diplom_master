function [index]=minmatrix(D,Ci,top)
% great !!!!!!�Ѿ����������ƵĲ���,Ч����ȫ��������.
% D����������������ϵ�֮��ľ������. Ci�ǲ�ͬ��֮������Ҫ���ӵı��� 
% index=(x_1,y_1;x_2,y_2,......;x_Ci,y_Ci)^T  �ֱ�������Ҫ�����Ԫ��λ��
% top=1;   %ÿ���ڵ����������������ϵ����߸���
[row,col]=size(D);
L=reshape(D,row*col,1);
[sortL,indexL]=sort(L);  %reshape����һ������������һ��.
index=[];iter=0; i=1; timesI=0;timesJ=0;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

while(iter<Ci)&(i<=row*col) %û�м����i�Ͻ������,һ����˵����������.���������ǰ��indeL(i)��λ��(indeI,indeJ)
    temp=indexL(i);
    remind=rem(temp,row);                                    %����
    exactdivide=floor(temp/row);                             %����
    if(remind==0)
        indexI=row; 
        indexJ=exactdivide;
    else
        indexI=remind;
        indexJ=exactdivide+1;
    end
%�ж�ÿ���ڵ����������ڵ�ĸ���,�����Ͻ�top�������� 
if ~isempty(index)
 block1=index(:,1);
 block2=index(:,2);
 timesI=sum(block1==indexI);
 timesJ=sum(block2==indexJ);
end  %if ~isempty(index)
%-------------------------------------------------------------- 
if (timesI<top)&&(timesJ<top)
index((iter+1),:)=[indexI indexJ];
iter=iter+1;                %���·���һ��.
end  % if(timesI<top)&&(timesJ<top) 
%������������,�͸���index.����ͼ�������һ������.
i=i+1;timesI=0;timesJ=0;
end % while(iter<Ci)
 



