function frames = getCustomMFCC(x)
%% Framing 
B = [1 -0.095];
y = filter(B, 1, x);

%% Framing
size = length(y);
frameDuration = 40/1000;
frameSize = frameDuration*fs;
frameStep = frameSize/4;
frameNumber = round(size/frameStep - frameSize/frameStep + 1);
frames = zeros(frameNumber, frameSize);

a=1;
b=frameSize;
for i=1:frameNumber
    d = a:b;
    frames(i,:) = y(d);
    a = a+frameStep;
    b = b+frameStep;
end

%% Hamming
h = hamming(frameSize);
for i=1:frameNumber
    frames(i,:) = frames(i,:).*h';
end

%% FFT transform
framesFFT = zeros(frameNumber, 257);
for i=1:frameNumber
    res = fft(frames(i,:), 512);
    framesFFT(i,:) = res(1:257);
end
frames = framesFFT;

%% PSD estimate
window = hamming(257);
periodograms = zeros(frameNumber, 257);
for i=1:frameNumber
    psd = periodogram(frames(i,:),window);
    periodograms(i,:) = psd(1:257, 1);
end
frames = periodograms;

%% Triag filter bank
f_min = 0;
f_low = 300;
f_high = 8000;
f_max = 0.5*fs;
M = 26;
K = 257;
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

%% Apply trifbank
results = zeros(frameNumber,M);
for i=1:frameNumber
    results(i,:) = frames(i,:)*H';
%     plot(results(i,:));
%     pause();
end
frames = results;

%% DCT
results = zeros(frameNumber, 13);
for i=1:frameNumber
    res = dct(frames(i,:));
    results(i,:) = res(:,1:13);
end
frames = results;
end