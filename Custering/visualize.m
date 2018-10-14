function visualize(sM, sD)

[q,t] = som_quality(sM,sD);
adm = som_distortion(sM,sD);
fprintf('QE = %f; TE = %f; ADM = %f;\n', q, t, adm);

figure;
som_show(sM,'umat','all');

% showHits(sM, H);

% if length(sM.codebook)>2
%     return;
% end
figure;
som_grid(sM,'Coord',sM.codebook(:,[1,2,3]),...
	 'Markersize',7,'Linecolor','k');
figure;
som_grid(sM,'Coord',sM.codebook(:,[1,2]),...
	 'Markersize',7,'Linecolor','k');
figure; hold on;
som_grid(sM,'Coord',sM.codebook(:,[1,2,3]),...
	 'Markersize',7,'Linecolor','k');
plot3(sD.data(:,1),sD.data(:,2),sD.data(:,3), 'r.');
% bmu = som_bmus(sM);
% co = sM.codebook(bmu,:);
% plot(co(1),co(2),'g+', 'Markersize', 10);
hold off;
end