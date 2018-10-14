clear all;
clc;

f = 0:44100;
m = 1125*log(1+f/700);

plot(f,m);