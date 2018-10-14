clear all;
close all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\SOM-Toolbox-master'));

n = 100;
d1 = rand(n,3)+1;
d2 = rand(n,3);
d = [d1;d2];

%% 
hold on;
plot(d1(:,1),d1(:,2),'.r');
plot(d2(:,1),d2(:,2),'.b');

%% 
d = table2array(Chainlink);
hold on;
plot3(d(1:end/2,1),d(1:end/2,2),d(1:end/2,3),'.r');
plot3(d(end/2+1:end,1),d(end/2+1:end,2),d(end/2+1:end,3),'.b');

%% 
epoches = [10 20 40 80 160];
threshes = [100];
rests = [100];
errors = [];
ind = 1;
for i=1:length(epoches)
    for j=1:length(rests)
        for k=1:length(threshes)
            sD = som_data_struct(d);
            sM = som_randinit(sD);
            sM = som_prototrain(sM,sD, epoches(1,i));
            [qe te] = som_quality(sM,sD);
        end
        %          ind=ind+2;
    end
    %      ind=ind+2;
    errors(i,ind)=qe;
    errors(i,ind+1)=te;
end
e1=errors;
save e1.mat e1;

close all
load e1.mat e1;
load e2.mat e2;
hold on
epoches = [10,20,40,80,160];
plot(epoches,e1(:,1),'-r');
plot(epoches,e2(:,1),'-b');
