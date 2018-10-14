function res = getY(sz)
res = zeros(sz,3);

a = 0.5;
b = 1;
r = (b-a).*rand(sz,1) + a;
res(:,1)=r;

a = 0.5;
b = 1;
r = (b-a).*rand(sz,1) + a;
res(:,2)=r;

a = 0;
b = 0.25;
r = (b-a).*rand(sz,1) + a;
res(:,3)=r;
end