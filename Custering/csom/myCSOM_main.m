clear all;
clc;
addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));
%% ��������� ������ ������

% ���������� ��������

trainDataStruct = loadDataStruct([1:2], [1], 'train', '_mfcc', '_orig',...
    'fsp', '_sub_5dB', 'dictor');
dim = 2;
% proj = getProjection('sammon', trainDataStruct.flat_data, dim);
m = trainDataStruct.flat_data';
samples_num = length(m(:,1));
%% 
% som_grid(sMt, 'coord', sMt.codebook(:,[1 2]))
% som_show(sMt,'umat', 'all')
%% ��������� ��������� ����
m=m';
w_num = 15; % ���������� �������� �� ������� �������

% ��� ��������� ������� ������ ������������� ������ � "�������"
sigma0 = 2; % �������� ������ �������������� �����������
tau1 = 8000/log(sigma0); % �������� ��������
nu0 = 0.15; % �������� ��������
% �������� �������� �������� �������� (��� ������ ��� �������� �������):
tau2 = 4000;

S = 0.05; % ����������� "�������"
% ������� ��� ����� "�������":
p = zeros(w_num, w_num);
b = zeros(w_num, w_num);

%% ������������� ������� �������������

% ������������� ���������� ����������
w = rand(w_num, w_num,length(m(:,1)))*0.1;

%% ������� ���� ��������

n = 0; % ���������� ����� ��������
while(1)
    
    % ����� �������� �������
    % (��������� �������� �� ��������� 1-samples_num)
    
    cur_sample_num = ceil(rand()*samples_num);
    
    % �������� ������ �� �������
    cur_sample = m(:,cur_sample_num);
    
    % ��������� ���������� �� �������� ������� �� ������� �������
    d = zeros(w_num,w_num,length(m(:,1)));
%     dist = zeros()
    for i=1:length(m(:,1))
        d(:,:,i) = bsxfun(@minus,w(:,:,i),cur_sample(i,1));
    end
    dist=(sqrt(d(:,:,1).^2+d(:,:,2).^2)-b) ; % ��������� ����������
    
    % ������� ��������� � �������� ������� ������ (������-����������)
    [A,B] = min(dist(:));
    % ������� ������� � ������� �� ��������� �������
    [I,J] = ind2sub([w_num, w_num],B);
    figure(1);
    
    % ����� ������������� ����� ��������� ���������:
    delta_sum = zeros(length(m(:,1)),1);
    
    sigma = sigma0*exp(-n/tau1); % �����������
    nu = nu0*exp(-n/tau2); % �������� ��������
    
    % ��������� �� ���� ��������
    for j=1: w_num
        for i=1: w_num
            % ������� ���������� ����� ������� � ���������� ��������� � �������
            % (� ������������ ��������)
            d_sq = (I-i)^2+(J-j)^2;
            % ������ �������������� ����������� �������
            h = exp(-d_sq/(2*(sigma^2)));
            
            % ������� ���� �������� �������
            w_cur = w(i,j,:);
            % ����������� ��� � ������ �������
            w_cur = reshape(w_cur,length(m(:,1)),1);
            % ������� �������� ��������� �����
            delta = nu*h*(cur_sample-w_cur);
            % ����� ����� ��������� ��������� ����� (��� �������� �������� ��������)
            delta_sum = delta_sum+delta.^2;
            % ���������� ���������
            w_cur = w_cur+delta;
            % ���������� ������������������ ��� �� �����
            w(i,j,:) = w_cur;
            % "�������" ����
            y=((I==i) && (J==j));
            p(i,j) = p(I,J)+0.5*(y-p(I,J));
            b(i,j) = S*(1.0/(w_num^2)-p(I,J));
        end;
    end;
    
    % ������ ������ ������ 100 ��������
    if(mod(n,100)==0)
        Plot_net(m,w);
    end;
      
    % �������� �������� ��������
    if (sqrt(delta_sum(1)+delta_sum(2))<1e-3)
        break;
    end;
    
    n=n+1; % ���������� ������ �������� (���������� ����� ��������)
    
end;
%% 
[w_rows, w_cols,dim] = size(w);
wNum = w_rows*w_cols;
 
pts = reshape(w,wNum,length(m(:,1)));

wt = zeros(size(w,2),size(w,1),size(w,3));
for  ii=1:w_rows
    for  jj=1:w_cols
        wt(ii,jj,:) = w(jj,ii,:);
    end;
end;

ptst = reshape(wt,wNum,length(m(:,1)));

sMap = som_map_struct(13,'msize',[15 14],'hexa');
sMap.codebook = ptst;
sMap.labels = cell(size(ptst));
sMap.topol.msize = [15 15];

som_grid(sMap, 'coord', sMap.codebook(:,[1 2 3]))
% som_show(sMap)

%% 
sprintf('������������� %d ��������',n);
