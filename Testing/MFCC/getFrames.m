function [frames, frameNumber, frameSize] = getFrames(y, fs, frameDuration, frameStep, wnd)
size = length(y);
frameSize = frameDuration*fs;
stepSize = frameStep*fs;
frameNumber = round(size/stepSize - frameSize/stepSize + 1);
frames = zeros(frameNumber, frameSize);

a=1;
b=frameSize;
for i=1:frameNumber
    d = a:b;
    frames(i,:) = y(d).*wnd;
    a = a+stepSize;
    b = b+stepSize;
end
end