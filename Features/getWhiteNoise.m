function y=getWhiteNoise(x, dB)
y=getSigma(x,dB).*wgn(1,length(x),0);
y=y+x;
end