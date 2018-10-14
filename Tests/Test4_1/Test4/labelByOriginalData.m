function sD_train = labelByOriginalData(sD, data_struct, cluster_number)
prefix = 'cluster_';
sD_train = sD;
len = length(sD.labels(:,1));
a = 1;
b = 0;
for i=1:cluster_number
    b = b + length(data_struct.data{i});
    sD_train = som_label(sD_train,'add',[a:b]',strcat(prefix,int2str(i)));
    a = b + 1;
end
end