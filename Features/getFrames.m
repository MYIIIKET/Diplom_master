function [frames, frameNumber, frameSize] = getFrames(y, fs, frameDuration, frameStep, wnd)
size = length(y);
frameSize = floor(frameDuration*fs);
stepSize = ceil(frameStep*fs);
frameNumber = round(size/stepSize - frameSize/stepSize + 1);
frames = zeros(frameNumber, frameSize);

a=1;
b=frameSize;
for i=1:frameNumber
    d = a:b;
    if d(end) > size
        break;
    end
    frames(i,:) = y(d).*wnd(1);
    a = a+stepSize;
    b = b+stepSize;
end
end