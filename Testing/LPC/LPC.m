clear all;
clc;

[x fs] = audioread('fsp2train1.wav');

y = filter([1 -0.97], 1, x);

hamm = hamming(fs*0.04);
[frames, framesNumber, frameSize] = getFrames(y, fs, 0.04, 0.01, hamm);

[lp,g] = lpc(frames', 27);
aFrames = zeros(framesNumber, 28);
for i=1:framesNumber
    aFrames(i,:) = g(i,1)./abs(fft(lp(i,:))).^2;
end
cFrames = zeros(framesNumber, 28);
for i=1:framesNumber
    cFrames(i,:) = -ifft(log(aFrames(i,:)));
end
% cFrames(any(isnan(cFrames), 2), :) = [];
% find(isnan(cFrames))