function [time_vs_session time_vs_session_all] = display_feature_dTvsdS_mean(var1,varT,varS,newfig,rev_color,showplot)


color_test1 = cbrewer('seq','YlGnBu',100); % colormap for corr 
color_test2 = color_test1(end:-1:1,:);

for mouse = 1:length(var1)
    clear time_vs_session_matrixes
    % Time from onset first day
    for i = 1:max(varS{mouse}(:))+1 % cylce through all possible interleaved session numbers
        for ii = 1:max(varT{mouse})+1 % cylce through all possible time intervals
            comparison_type_index = (varS{mouse} == i-1).*(varT{mouse}== ii-1);
            comparison_type_index = triu(comparison_type_index);
            comparison_type_index(comparison_type_index==0) = NaN;
            time_vs_session_matrixes{i,ii} = var1{mouse}.*comparison_type_index;
        end
    end
    
    time_vs_session_all{mouse} = cellfun(@(v) (v(~isnan(v(:)))),time_vs_session_matrixes,'UniformOutput',false);
    time_vs_session{mouse} = cellfun(@(v) nanmean(v(:)),time_vs_session_matrixes,'UniformOutput',true);
    time_vs_session{mouse} = real(time_vs_session{mouse}); % for some reason it sometimes adds a 0 imaginary part part to the real number
    
    if showplot
        if newfig
            figure;
        end
        h = imagesc(time_vs_session{mouse});
        if newfig
            title(['Mouse ' num2str(mouse)]);
        end
        if rev_color
            colormap(gca,color_test2);
        else
            colormap(gca,color_test1);
        end
        hold on
        set(gca,'TickDir','out'); xlabel('Delta Time'); ylabel('Delta Sessions'); set(gca,'YDir','normal')
        colorbar; set(h,'AlphaData',~isnan(time_vs_session{mouse}));  box off
        pbaspect([1 1 1]);
        tickmarksy=find(double(~isnan(nanmean(time_vs_session{mouse},2))));
        tickmarksx=find(double(~isnan(nanmean(time_vs_session{mouse},1))));
        
        set(gca, 'XTick', [0:3:tickmarksx(end)], 'XTickLabel', [0:3:tickmarksx(end)]-1); % axis is too dense otherwise
        set(gca, 'YTick', [0:3:tickmarksy(end)], 'YTickLabel', [0:3:tickmarksy(end)]-1)
        
        set(gcf,'color','w'); set(gca,'YDir','normal'); pbaspect([1 1 1]);
    end
end