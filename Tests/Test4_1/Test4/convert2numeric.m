function labels = convert2numeric(sM, cluster_number)
labels = som_label2num(sM);
if max(labels)>cluster_number
    labels(labels(:)==max(labels))=0;
end
end