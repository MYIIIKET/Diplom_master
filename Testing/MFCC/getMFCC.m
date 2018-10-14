function c = getMFCC(x, fs)
y = filter([1 -0.97], 1, x);

%% Framing
hamm = hamming(fs*0.04);
[frames, ~, frameSize] = getFrames(y, fs, 0.04, 0.01, hamm);

%% FFT
fftFrames = abs(fft(frames'));
%% PS
psFrames = (fftFrames.^2)./frameSize;

%% Triangle
H = getTriangleFilters(fs, 300, 8000, 40, frameSize);
filterFrames = H*psFrames(1:frameSize,:);

%% DCT
dctm = @(N,M)(sqrt(2.0/M)*cos(repmat([0:N-1].',1,M).*repmat(pi*([1:M]-0.5)/M,N,1)));
dct = dctm(13, 40);
filterFrames(:,all(filterFrames==0))=[];
c = dct*log(filterFrames);
end