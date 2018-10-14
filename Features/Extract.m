clear all;
clc;
addpath(genpath('C:\Users\MYlll\Desktop\Tools\rastamat'));

%% 
path = 'G:\grant\wav';

pattern = '.wav';
allFiles = dir(path);
fileNames = {allFiles.name};
mkdir mfcc1;

for i=3:length(fileNames)
    name = fileNames{i};
    [~,~,ext] = fileparts(name);
    if strcmp(pattern, ext) == 1
        [x, fs] = audioread(strcat(path,'\',name));
        data = melfcc(x, fs, 'dither', 1);
        name = strrep(name, pattern, '');
        save(strcat('mfcc1\',strcat(name,'.mat')),'data');
    end
end
% plot(res)
% sound(res,fs)
% f1_test = melfcc(res, fs, 'dither', 1);
% f2 = rastaplp(res, fs);
% res = rastaplp(res, fs, 0, 12);
% f1_test = lpc2cep(getLPC(res,fs)',28);
% del = deltas(f1_test);
% ddel = deltas(deltas(f1_test,5),5);
% f1_test = [f1_test;del;ddel];
