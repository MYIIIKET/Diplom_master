function H = getTriangleFilters(fs, f_low, f_high, M, K)
f_min = 0;
f_max = 0.5*fs;
hz2mel = @(hz)(1127*log(1+hz/700));
mel2hz = @(mel)(700*exp(mel/1127)-700);
f = linspace( f_min, f_max, K );
c = mel2hz( hz2mel(f_low)+[0:M+1]*((hz2mel(f_high)-hz2mel(f_low))/(M+1)) );
H = zeros( M, K );

for m = 1:M
    k = f>=c(m)&f<=c(m+1);
    H(m,k) = 2*(f(k)-c(m)) / ((c(m+2)-c(m))*(c(m+1)-c(m)));
    k = f>=c(m+1)&f<=c(m+2);
    H(m,k) = 2*(c(m+2)-f(k)) / ((c(m+2)-c(m))*(c(m+2)-c(m+1)));
end

H = H./repmat(max(H,[],2),1,K);
end