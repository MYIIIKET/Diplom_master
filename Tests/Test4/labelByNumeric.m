function sM = labelByNumeric(sM, labels)
prefix = 'cluster_';
for i=1:length(labels)
    sM.labels{i} = strcat(prefix,int2str(labels(i)));
end
end