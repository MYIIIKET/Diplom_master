clear all;
close all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\rastamat'));
addpath(genpath('C:\Users\MYlll\Desktop\Tools\mfcc\mfcc'));

[x fs] = audioread('fsp1train5.wav');

B = [1 -0.97];
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


%% Compute power spectrum
result = zeros(frameNumber, 513);
for i=1:frameNumber
    result(i,:) = periodogram(fftFrames(i,:),[],length(fftFrames(i,:)),fs);
end
frames = result;

MAG = abs(fft(frames,1024,1));

figure;
plot(frames(:,1));
figure;
plot(MAG(:,1));
%% Triag filter bank
f_min = 0;
f_low = 300;
f_high = 8000;
f_max = 0.5*fs;
M = 26;
K = 513;
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
end
frames = results;
%% DCT
cNumber = 13;
results = zeros(frameNumber, cNumber);
for i=1:26
    res = dct(frames(i,:));
    results(i,:) = res(1,1:cNumber);
end
CC = results';
%% test

[cepstra,aspectrum,pspectrum] = melfcc(x, fs,...
    'minfreq', 300, 'maxfreq', 8000, 'numcep', 13,...
'nbands', 40, 'wintime', 0.04, 'hoptime', 0.01);

[ MFCCs, FBEs, frames ] = mfcc( x, fs, 40, 10, 0.97, @hamming, [300 8000], 26, 12+1, 22 );


% ceplifter = @( N, L )( 1+0.5*L*sin(pi*[0:N-1]/L) );
% lifter = ceplifter( 13, 22 );
% CC = diag(lifter)*CC;

t = 42;
figure;
plot(CC(:,t));
                 
figure;
plot(cepstra(:,t));
                                