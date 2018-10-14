clear all;
clc;
addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));
%% Параметры набора данных

% Количество примеров

trainDataStruct = loadDataStruct([1:2], [1], 'train', '_mfcc', '_orig',...
    'fsp', '_sub_5dB', 'dictor');
dim = 2;
% proj = getProjection('sammon', trainDataStruct.flat_data, dim);
m = trainDataStruct.flat_data';
samples_num = length(m(:,1));
%% 
% som_grid(sMt, 'coord', sMt.codebook(:,[1 2]))
% som_show(sMt,'umat', 'all')
%% Параметры нейронной сети
m=m';
w_num = 15; % количество нейронов на сторону решетки

% для двумерной решетки обычно устанавливают равным её "радиусу"
sigma0 = 2; % параметр ширины топологической окрестности
tau1 = 8000/log(sigma0); % параметр убывания
nu0 = 0.15; % скорость обучения
% Скорость убывания скорости обучения (чем больше тем медленее убывает):
tau2 = 4000;

S = 0.05; % Коэффициент "совести"
% массивы для учета "совести":
p = zeros(w_num, w_num);
b = zeros(w_num, w_num);

%% Инициализация весовых коэффициентов

% Инициализация случайными значениями
w = rand(w_num, w_num,length(m(:,1)))*0.1;

%% Главный цикл обучения

n = 0; % Дискретное время процесса
while(1)
    
    % Номер текущего примера
    % (случайное значение из диапазона 1-samples_num)
    
    cur_sample_num = ceil(rand()*samples_num);
    
    % Извлекли пример из выборки
    cur_sample = m(:,cur_sample_num);
    
    % Вычисляем расстояния от текущего примера до каждого нейрона
    d = zeros(w_num,w_num,length(m(:,1)));
%     dist = zeros()
    for i=1:length(m(:,1))
        d(:,:,i) = bsxfun(@minus,w(:,:,i),cur_sample(i,1));
    end
    dist=(sqrt(d(:,:,1).^2+d(:,:,2).^2)-b) ; % евклидово расстояние
    
    % Находим ближайший к входному вектору нейрон (нейрон-победитель)
    [A,B] = min(dist(:));
    % Находим индексы в матрице по сквозному индексу
    [I,J] = ind2sub([w_num, w_num],B);
    figure(1);
    
    % здесь накапливается сумма квадратов коррекции:
    delta_sum = zeros(length(m(:,1)),1);
    
    sigma = sigma0*exp(-n/tau1); % окрестность
    nu = nu0*exp(-n/tau2); % скорость обучения
    
    % Пробегаем по всем нейронам
    for j=1: w_num
        for i=1: w_num
            % Квадрат расстояния между текущим и победившим нейронами в решетке
            % (в пространстве индексов)
            d_sq = (I-i)^2+(J-j)^2;
            % ширина топологической окрестности нейрона
            h = exp(-d_sq/(2*(sigma^2)));
            
            % Достаем веса текущего нейрона
            w_cur = w(i,j,:);
            % Преобразуем его в вектор столбец
            w_cur = reshape(w_cur,length(m(:,1)),1);
            % Находим величину коррекции весов
            delta = nu*h*(cur_sample-w_cur);
            % Копим сумму квадратов изменений весов (для проверки критерия останова)
            delta_sum = delta_sum+delta.^2;
            % Производим коррекцию
            w_cur = w_cur+delta;
            % Возвращаем откорректированный вес на место
            w(i,j,:) = w_cur;
            % "совесть" сети
            y=((I==i) && (J==j));
            p(i,j) = p(I,J)+0.5*(y-p(I,J));
            b(i,j) = S*(1.0/(w_num^2)-p(I,J));
        end;
    end;
    
    % Строим график каждые 100 итераций
    if(mod(n,100)==0)
        Plot_net(m,w);
    end;
      
    % Проверка критерия останова
    if (sqrt(delta_sum(1)+delta_sum(2))<1e-3)
        break;
    end;
    
    n=n+1; % увеличение номера итерации (дискретное время процесса)
    
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
sprintf('Потребовалось %d итераций',n);
