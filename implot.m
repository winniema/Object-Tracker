function[] = implot(frame, vid, timepos_r, timepos_g)
figure; imagesc(read(vid, frame)); axis image;
hold on
xlim([0 1280])
ylim([0 720])
scatter(timepos_r(:,2), timepos_r(:,3), 'm+')
scatter(timepos_g(:,2), timepos_g(:,3), 'c+')
set(gca,'xaxislocation','top','yaxislocation','left','ydir','reverse')
hold off
end