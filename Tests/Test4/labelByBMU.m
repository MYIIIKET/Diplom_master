function sM = labelByBMU(sM, data_struct)
prefix = 'cluster_';
for i=1:length(data_struct.data)
    bmus = som_bmus(sM,data_struct.data{i});
    for j=1:length(bmus)
        sM.labels{bmus(j)} = strcat(prefix,int2str(i));
    end
end
end