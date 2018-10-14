function res = getProjection(arg, data, dim)
if strcmp(arg, 'sammon')
    res = sammon(data', dim, 10^(-5), '3');
%     res = compute_mapping(data', 'Sammon', dim);
    return;
end
if strcmp(arg, 'pcaproj')
    res = pcaproj(data', dim);
    return;
end
if strcmp(arg, 'cca')
    pca = pcaproj(data', dim);
    res = cca(data', pca, 100);
    return;
end
end