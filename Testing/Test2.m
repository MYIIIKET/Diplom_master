clear all;
close all;
clc;

creaturesNumber = 100;
speed = 5;

radius = 5;
border = 100;
range = 100;
creatures = struct([]);


position = round(range*rand(creaturesNumber, 2));

target = round(200*rand(creaturesNumber, 2)-50);


f = figure;
hold on;
for j=1:creaturesNumber
    creatures(j).x = position(j,1);
    creatures(j).y = position(j,2);
    creatures(j).speed = speed;
    plot(creatures(j).x, creatures(j).y, 'r.');
    plot(target(:,1), target(:,2), 'bo');
end
pause();
hold off;

direction = zeros(1, 2);
for i=1:100
    for j=1:creaturesNumber
        
        distanceX = (target(j,1) - creatures(j).x)^2;
        distanceY = (target(j,2) - creatures(j).y)^2;
        distance = sqrt(distanceX + distanceY);
        
        directionX = (target(j,1) - creatures(j).x)/distance;
        directionY = (target(j,2) - creatures(j).y)/distance;
        
        
        if abs(distance)<speed
            creatures(j).x = target(j,1);
            creatures(j).y = target(j,2);
            creatures(j).speed = 0;
        else
            creatures(j).x = creatures(j).x+creatures(j).speed*directionX;
            creatures(j).y = creatures(j).y+creatures(j).speed*directionY;
        end
    end
    
    clf(f);
    hold on;
    for j=1:creaturesNumber
        plot(creatures(j).x, creatures(j).y, 'r.');
        plot(target(:,1), target(:,2), 'bo');
    end
    hold off;
    pause(0.00001);
end



