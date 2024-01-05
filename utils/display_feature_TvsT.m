function display_feature_TvsT_mean(var1,newfig,rev_color)


color_test1 = cbrewer('seq','YlGnBu',100); % colormap for corr 
color_test2 = color_test1(end:-1:1,:);

for mouse = 1:size(var1,2)
    if newfig
        figure;
    end
    h = imagesc(var1{mouse}); 
    if newfig
        title(['Mouse ' num2str(mouse)]);
    end
    if rev_color
        colormap(gca,color_test2); 
    else
        colormap(gca,color_test1); 
    end
    hold on
    set(gca,'TickDir','out'); xlabel('Time'); ylabel('Time'); set(gca,'YDir','normal')
    colorbar; set(h,'AlphaData',~isnan(var1{mouse}));  box off
    pbaspect([1 1 1]);
    tickmarks=find(double(~isnan(nanmean(var1{mouse},2))));
%     xticks([tickmarks]); yticks([tickmarks]);
    
    set(gca, 'XTick', [1:3:tickmarks(end)], 'XTickLabel', [1:3:tickmarks(end)]); % axis is too dense otherwise
    set(gca, 'YTick', [1:3:tickmarks(end)], 'YTickLabel', [1:3:tickmarks(end)])
    
    set(gcf,'color','w'); set(gca,'YDir','normal'); pbaspect([1 1 1]);
end

end