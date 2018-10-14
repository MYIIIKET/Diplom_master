clear all;
close all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\rastamat'));
addpath(genpath('C:\Users\MYlll\Desktop\Tools\mfcc\mfcc'));

[x fs] = audioread('fsp2train1.wav');

y = filter([1 -0.97], 1, x);


%% Framing
hamm = hamming(fs*0.04);
[frames, framesNumber, frameSize] = getFrames(y, fs, 0.04, 0.01, hamm);

%% FFT
fftFrames = abs(fft(frames'));
%% PS
psFrames = (fftFrames.^2)./frameSize;

%% Triangle
H = getTriangleFilters(fs, 300, 8000, 40, frameSize);
filterFrames = H*psFrames(1:frameSize,:);

%% DCT
find(filterFrames(6,:) == 0)
% filterFrames(:,all(filterFrames==0))=[];
c = -dct(log(filterFrames(1:13,:)));

%% Deltas
d = deltas(c);
dd = deltas(deltas(c,5),5);
c = [c;d;dd];