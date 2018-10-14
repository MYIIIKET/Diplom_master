clear all;
clc;

addpath(genpath('C:\Users\MYlll\Desktop\Tools\rastamat'));

F1 = 'C:\Users\MYlll\Desktop\Features\dev-clean\1673\143396\1673-143396-0004.flac'
F2 = 'C:\Users\MYlll\Desktop\Features\dev-clean\2035\147960\2035-147960-0010.flac'
[x, fs] = audioread(F2);
mfcc = getMFCC(x, fs);
s.('F2') = mfcc;
save('F2.mat', '-struct', 's');


% for i=1:length(dirs)
%     D = [];
%     allPath = dir(dirs{i});
%     allNames = {allPath.name};
%     allFolders = {allPath.folder};
%     for j=3:length(allNames)
%         [x, fs] = audioread(strcat(allFolders{j},'\',allNames{j}));
%         mfcc = getMFCC(x, fs);
% %         mfcc = melfcc(x, fs, 'minfreq', 300, 'maxfreq', 8000, 'numcep', 13, 'nbands', 40, 'fbtype', 'fcmel', 'dcttype', 1, 'usecmp', 1, 'wintime', 0.04, 'hoptime', 0.01, 'preemph', 0.97, 'dither', 1);
%         D = [D'; mfcc']';
%     end
%     d = 'cello';
%     p = strcat(d,'\');
%     mkdir(d);
%     n = strcat(strcat(d,'_'),folders{i});
%     s = struct;
%     s.(n) = D;
%     save(strcat(p,'\',n,'.mat'), '-struct', 's');
% end