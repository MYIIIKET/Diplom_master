clear all;
close all
clc;

creaturesNumber = 1000;
range = 100;
speed = 0.1;

positions = (range*rand(creaturesNumber, 2));
directions = zeros(creaturesNumber, 2);
targets = (range*rand(creaturesNumber, 2));
classes = zeros(creaturesNumber, 1);
classes(1,1) = 1;

a = [0,1];

distances = util.getDistances(positions, positions);

distances2 = util.getDistances(positions, targets);



e = 1;
r = 5;
ind = util.getNeighbors(e, r, distances);

f = figure;
hold on;
plot(positions(:,1), positions(:,2), '.');
plot(targets(:,1), targets(:,2), 'o');
pause();
for i=1:1000
    distance = diag(distances2);

    directions(:,1) = (targets(:,1) - positions(:,1))./distance;
    directions(:,2) = (targets(:,2) - positions(:,2))./distance;
    
    positions(:,1) = positions(:,1)+speed.*directions(:,1);
    positions(:,2) = positions(:,2)+speed.*directions(:,2);

    
    distances = util.getDistances(positions, positions);
    distances2 = util.getDistances(positions, targets);
    clf(f);
    hold on;
    plot(positions(:,1), positions(:,2), '.');
    plot(targets(:,1), targets(:,2), 'o');
    pause(0);
    hold off;
end


% hold on
% plot(positions(e,1), positions(e,2), 'go');
% plot(positions(:,1), positions(:,2), '.');
% plot(positions(ind(:,2),1), positions(ind(:,2),2), 'o');