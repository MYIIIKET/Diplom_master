clear all;
close all;
clc;

creaturesNumber = 100;
range = 100;
field = zeros(range);
speed = 5;

creatures(creaturesNumber,1) = creature;
position = round(range*rand(creaturesNumber, 2));
target = round(200*rand(creaturesNumber, 2)-50);

f = figure;
hold on;
for i=1:creaturesNumber
    creatures(i).x = position(i,1);
    creatures(i).y = position(i,2);
    creatures(i).speed = speed;
    creatures(i).class = class.getRandom;
    creatures(i).draw(f);
end
plot(target(:,1), target(:,2), 'bo');
pause();
clf(f);

hold on;
for i=1:100
    clf(f);
    hold on;
    plot(target(:,1), target(:,2), 'bo');
    
    for j=1:creaturesNumber
        creatures(j).moveTo(target(j,:));
        creatures(j).draw(f);
    end
    pause(0.5);
end
