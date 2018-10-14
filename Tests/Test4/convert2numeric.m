function labels = convert2numeric(sM, cluster_number)
label_size = length(sM.labels);
labels = zeros(label_size,1);
% labels = som_label2num(sM);
for i=1:label_size
    if strcmp(sM.labels{i},'cluster_1')
        labels(i) = 1;
    elseif strcmp(sM.labels{i},'cluster_2')
            labels(i) = 2;
    else
        labels(i) = 0;
    end
end
% if max(labels)>cluster_number
%     labels(labels(:)==max(labels))=0;
% end
end