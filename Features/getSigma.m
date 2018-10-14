function sigma=getSigma(x, SNR)
% noise generation
Ps = 10*log10(std(x).^2);       % signal power, dBV^2
Pn = Ps - SNR;                  % noise power, dBV^2
Pn = 10^(Pn/10);                % noise power, V^2
sigma = sqrt(Pn);               % noise RMS, V
end