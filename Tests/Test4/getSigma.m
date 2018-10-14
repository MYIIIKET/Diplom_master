function sigma = getSigma(U, delta, iters)
sigma = 0;
sigma_arr = zeros(iters,1);
sigma_arr(1,1)=sigma;
for i=1:iters
    urp = sqrt((((sigma-min(U(:)))/(max(U(:))-min(U(:))))^2)+(1-ecdf(sigma)).^2);
    sigma_arr(i,1)=min(urp);
    sigma=sigma+delta;
end
sigma = max(sigma_arr);
end