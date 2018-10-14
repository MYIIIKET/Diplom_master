function showHits(sM, H)
figure;
colormap(gray);
som_show(sM,'umat','all')

hitFields = fieldnames(H.hits);
colorFields = fieldnames(H.colors);
fieldNumber = length(hitFields);

for i=1:fieldNumber
    som_show_add('hit',H.hits.(hitFields{i}), 'Subplot', 1,...
        'Markercolor', H.colors.(colorFields{i}));
end
end