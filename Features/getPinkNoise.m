function y=getPinkNoise(x, dB)
points = length(x);
if rem(points, 2)
    M = points+1;
else
    M = points;
end
wNoise = randn(1, M);
X = fft(wNoise);
NumUniquePts = M/2 + 1;     % number of the unique fft points
n = 1:NumUniquePts;         % vector with frequency indexes 

X = X(1:NumUniquePts);      
X = X./sqrt(n);

X = [X conj(X(end-1:-1:2))];
y = real(ifft(X));


y = y(1, 1:points);
y = y - mean(y);
y = y/std(y, 1);
y=getSigma(x,dB).*y;
y = x + y;
end