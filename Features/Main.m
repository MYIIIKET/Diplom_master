clear all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\rastamat'));
addpath(genpath('C:\Users\MYlll\Desktop\Features'));
% addpath(genpath('C:\Users\MYlll\Desktop\Tools\AuditoryToolbox'));
% addpath(genpath('C:\Users\MYlll\Desktop\Tools\drtoolbox'));
% addpath(genpath('C:\Users\MYlll\Desktop\Features\SpeakerRecognition'));

cName = {'_mfcc'};
noise = {'_orig'};
op = {'_add_5dB', '_sub_5dB', '_0dB'};

% op = {};


pattern = '.wav';
replacement = '';

allFiles = dir('raw');
cd raw
allDirs = {allFiles.name};
for i=3:24
    cd(allDirs{i});
    names = dir('.');
    files = {names.name};
    n = length(names);
    for j=3:n
        name = files{j};
        [~,~,ext] = fileparts(name);
        if strcmp(pattern, ext) == 1
            [basex, fs] = audioread(name);
            name = regexprep(name,pattern,replacement);
            varName = name;
            for t1 = 1:length(cName)
                x=basex;
                varName = strcat(varName,cName{t1});
                for t2 = 1:length(noise)
                    x=basex;
                    varName = strcat(varName,noise{t2});
                    for t3 = 1:length(op)
                        x=basex;
                        if regexp(varName, '_orig')
                            if regexp(varName, '_data')
                                variable.(varName) = x;
                            end
                            if regexp(varName, '_mfcc')
%                                 mfcc = melfcc(x, fs, 'minfreq', 300, 'maxfreq', 8000, 'numcep', 13, 'nbands', 40, 'fbtype', 'fcmel', 'dcttype', 1, 'usecmp', 1, 'wintime', 0.04, 'hoptime', 0.01, 'preemph', 0.97, 'dither', 1);
                                mfcc = [];
                                mfcc = rastaplp(x, fs, 0, 12);
                                del = deltas(mfcc);
                                ddel = deltas(deltas(mfcc,5),5);
                                mfcc = [mfcc;del;ddel];
%                                 mfcc = melfcc(x, fs, 'minfreq', 0, 'maxfreq', 8000, 'numcep', 20, 'nbands', 40, 'fbtype', 'fcmel', 'dcttype', 1, 'usecmp', 1, 'wintime', 0.032, 'hoptime', 0.016, 'preemph', 0, 'dither', 1);
                                variable.(varName) = mfcc;
                            end
                            if regexp(varName, '_lpc')
                                y = powspec(x, fs, 0.04, 0.01, 1);
                                aspectrum = audspec(y, fs, 21, 'bark', 300, 8000, 0, 1);
                                lpc = dolpc(aspectrum,20);
                                variable.(varName) = lpc;
                            end
                            if regexp(varName, '_plp')
                                plp = melfcc(x, fs, 'minfreq', 300, 'maxfreq', 8000, 'numcep', 21,...
                                    'nbands', 20, 'fbtype', 'bark', 'dcttype', 1, 'usecmp', 1,...
                                    'wintime', 0.04, 'hoptime', 0.01, 'preemph', 0, 'dither', 1,...
                                    'minfreq', 300, 'maxfreq', 8000);
                                variable.(varName) = plp;
                            end
                        end
                        x=basex;
                        varName = strcat(varName,op{t3});
                        if regexp(varName, '_white')
                            if regexp(varName, '_add_5dB')
                                x=getWhiteNoise(x', +5); x=x';
                            end
                            if regexp(varName, '_sub_5dB')
                                x=getWhiteNoise(x', -5); x=x';
                            end
                            if regexp(varName, '_0dB')
                                x=getWhiteNoise(x', 0); x=x';
                            end
                            if regexp(varName, '_data')
                                variable.(varName) = x;
                            end
%                             if regexp(varName, '_mfcc')
%                                 mfcc = melfcc(x, fs, 'minfreq', 300, 'maxfreq', 8000, 'numcep', 21,...
%                                     'nbands', 20, 'fbtype', 'fcmel', 'dcttype', 1, 'usecmp', 1,...
%                                     'wintime', 0.04, 'hoptime', 0.01, 'preemph', 0, 'dither', 1,...
%                                     'minfreq', 300, 'maxfreq', 8000);
%                                 variable.(varName) = mfcc;
%                             end
                            if regexp(varName, '_lpc')
                                y = powspec(x, fs, 0.04, 0.01, 1);
                                aspectrum = audspec(y, fs, 21, 'bark', 300, 8000, 0, 1);
                                lpc = dolpc(aspectrum,20);
                                variable.(varName) = lpc;
                            end
                            if regexp(varName, '_plp')
                                plp = melfcc(x, fs, 'minfreq', 300, 'maxfreq', 8000, 'numcep', 21,...
                                    'nbands', 20, 'fbtype', 'bark', 'dcttype', 1, 'usecmp', 1,...
                                    'wintime', 0.04, 'hoptime', 0.01, 'preemph', 0, 'dither', 1,...
                                    'minfreq', 300, 'maxfreq', 8000);
                                variable.(varName) = plp;
                            end
                        end
                        x=basex;
                        if regexp(varName, '_pink')
                            if regexp(varName, '_add_5dB')
                                x=getPinkNoise(x', +5); x=x';
                            end
                            if regexp(varName, '_sub_5dB')
                                x=getPinkNoise(x', -5); x=x';
                            end
                            if regexp(varName, '_0dB')
                                x=getPinkNoise(x', 0); x=x';
                            end
                            if regexp(varName, '_data')
                                variable.(varName) = x;
                            end
%                             if regexp(varName, '_mfcc')
%                                 mfcc = melfcc(x, fs, 'minfreq', 300, 'maxfreq', 8000, 'numcep', 21,...
%                                     'nbands', 20, 'fbtype', 'fcmel', 'dcttype', 1, 'usecmp', 1,...
%                                     'wintime', 0.04, 'hoptime', 0.01, 'preemph', 0, 'dither', 1,...
%                                     'minfreq', 300, 'maxfreq', 8000);
%                                 variable.(varName) = mfcc;
%                             end
                            if regexp(varName, '_lpc')
                                y = powspec(x, fs, 0.04, 0.01, 1);
                                aspectrum = audspec(y, fs, 21, 'bark', 300, 8000, 0, 1);
                                lpc = dolpc(aspectrum,20);
                                variable.(varName) = lpc;
                            end
                            if regexp(varName, '_plp')
                                plp = melfcc(x, fs, 'minfreq', 300, 'maxfreq', 8000, 'numcep', 21,...
                                    'nbands', 20, 'fbtype', 'bark', 'dcttype', 1, 'usecmp', 1,...
                                    'wintime', 0.04, 'hoptime', 0.01, 'preemph', 0, 'dither', 1,...
                                    'minfreq', 300, 'maxfreq', 8000);
                                variable.(varName) = plp;
                            end
                        end
                        varName = regexprep(varName, op{t3},'');
                    end
                    varName = regexprep(varName, noise{t2},'');
                end
                varName = regexprep(varName, cName{t1},'');
            end
            j=3;
            names = dir('.');
            n = length(names);
        end
    end
    cd ..
end
cd ..


Data.Data = variable;
mkdir flatData
cd flatData

fields = fieldnames(Data.Data);
variable = Data.Data;
for i=1:numel(fields)
    res = variable.(fields{i});
    save(strcat(fields{i},'.mat'), 'res');
end
cd ../



