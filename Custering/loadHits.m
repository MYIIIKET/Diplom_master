function sD = loadHits(sM, sD)

colors = [  [0 0 1];
            [0 1 0];
            [0 1 1];
            [1 0 0];
            [1 0 1];
            [1 1 0];
            [0 0 0.5];
            [0 0.5 0];
            [0 0.5 0.5];
            [0.5 0 0];
            [0.5 0 0.5];
            [0.5 0.5 0]];

H = struct;

lab = unique(sD.labels(:,1));
for i=1:length(lab)
    d=find(strcmp(sD.labels(:,1),lab{i}))';
    h = som_hits(sM,sD.data(d,:));
    color = colors(i,:);
    H.('hits').(strcat('h', num2str(i))) = h;
    H.('colors').(strcat('color',num2str(i))) = color;
end

sD.H = H;
end