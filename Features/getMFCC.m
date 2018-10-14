function c = getMFCC(x, fs)
y = filter([1 -0.97], 1, x);

%% Framing
hamm = hamming(fs*0.0004);
[frames, ~, frameSize] = getFrames(y, fs, 0.0004, 0.0002, hamm);

%% FFT
fftFrames = abs(fft(frames'));
%% PS
psFrames = (fftFrames.^2)./frameSize;

%% Triangle
H = getTriangleFilters(fs, 300, 8000, 40, frameSize);
filterFrames = H*psFrames(1:frameSize,:);

%% DCT
filterFrames(:,all(filterFrames==0))=[];
c = -dct(log(filterFrames(1:13,:)));

%% Deltas
% d = deltas(c);
% dd = deltas(deltas(c,5),5);
% c = [c;d;dd];
end