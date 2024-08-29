%% calculate pixel-wise correlations of roi crops
data_stability = roi_correlations(data_stability,2,'stability');

%% plot settings
set(groot,'DefaultAxesTickDir', 'out')
set(groot,'DefaultAxesTickDirMode', 'manual');
set(groot,'DefaultAxesFontWeight', 'normal');
set(groot,'DefaultAxesFontName', 'Arial');
set(groot,'DefaultAxesFontSizeMode', 'manual')
set(groot,'DefaultAxesFontSize', 12)
set(groot,'DefaultAxesTitleFontSizeMultiplier',1.2)
set(groot,'DefaultAxesTickLength', [0.01 0.01]);
set(groot,'DefaultFigureColor','w')
set(groot,'DefaultLineLineWidth',2)
set(groot,'DefaultAxesLineWidth',1)
set(groot,'DefaultHistogramFaceColor',[0.5 0.5 0.5])
set(groot,'DefaultScatterMarkerEdgeColor',[0.5 0.5 0.5])
set(groot,'DefaultScatterMarkerFaceColor',[0.5 0.5 0.5])
set(groot,'DefaultPolaraxesFontSizeMode', 'manual')
set(groot,'DefaultPolaraxesFontWeight', 'normal');
set(groot,'DefaultAxesFontName', 'Arial');
set(groot,'DefaultPolaraxesFontSize',10)
set(groot,'DefaultPolaraxesTitleFontSizeMultiplier',1.2);
set(groot,'DefaultTextFontName', 'Arial');
set(groot,'DefaultLegendFontName', 'Arial');
set(groot,'DefaultColorbarFontName', 'Arial');
set(groot,'DefaultColorbarTickDirection','out');
set(groot,'DefaultFigureRenderer','painter');

savefigs = 1;
save_fig_loc = 'D:\stability_figures';
fig_number = 0;
close all
clear fig_hand
%% ROI correlations stability data
roi_corr_threshold = [-1 0.1 0.25 0.5];
data = data_stability;
mouse_group =[1,2,3,4,5];
% data = data_orientation_deprivation;
% mouse_group = [1:6 9 10 11 14 16 17 18];

% plot distribution of min roi correlations with threshold
clear correlations collect_rois
cell_counter = 0;
for mouse_n = mouse_group
    for cell_n = 1:length(data{mouse_n})
        cell_counter = cell_counter+1;
        roi_corrs = data{mouse_n}(cell_n).ROI_correlations;
        roi_corrs(roi_corrs==0)=nan;
        cell_min_cor = nanmin(roi_corrs(:));
        correlations(cell_counter) = cell_min_cor;
        
        [x,y] = find(roi_corrs==cell_min_cor);
        idx = [x(1),y(1)];
        if isfield(data{mouse_n}(cell_n),'suite2P_template')
            day_1_roi = data{mouse_n}(cell_n).suite2P_template{idx(1)};
            day_2_roi = data{mouse_n}(cell_n).suite2P_template{idx(2)};
        else
            day_1_roi = data{mouse_n}(cell_n).roi_template{idx(1)};
            day_2_roi = data{mouse_n}(cell_n).roi_template{idx(2)};
        end
        
        roi_crop_size = size(day_1_roi,1);
        % take average
        day_1_roia = day_1_roi(1:roi_crop_size,1:roi_crop_size,:);
        day_2_roia = day_2_roi(1:roi_crop_size,1:roi_crop_size,:);
        try
            day_1_roib = day_1_roi(1:roi_crop_size,roi_crop_size+1:roi_crop_size*2,:);
        catch
            day_1_roib = nan(size(day_1_roia));
        end
        try
            day_2_roib = day_2_roi(1:roi_crop_size,roi_crop_size+1:roi_crop_size*2,:);
        catch
            day_2_roib = nan(size(day_2_roia));
        end
        day_1_roi = nanmean(cat(4,day_1_roia,day_1_roib),4);
        day_2_roi = nanmean(cat(4,day_2_roia,day_2_roib),4);
        
        collect_rois{cell_counter} = cat(2,day_1_roi,day_2_roi);
    end
end
for cell_n = 1:cell_counter
   temp = collect_rois{cell_n};
   if size(temp,1)==25 && size(temp,2)==50
       temp(:,:,1) = rescale(temp(:,:,1));
       temp(:,:,2) = rescale(temp(:,:,2));
       temp(:,:,3) = zeros(size(temp(:,:,1)));
       collect_rois{cell_counter} = temp;
   elseif size(temp,1)==41 && size(temp,2)==82
       temp = imresize(temp,[25 50]);
       collect_rois{cell_counter} = temp;
   else
       collect_rois{cell_counter} = nan(25,50,3);
   end
end

correlations_all = cellfun(@(x) {x(:).ROI_correlations},data(mouse_group),'UniformOutput',0);
correlations_all = cat(2,correlations_all{:});
correlations_all = cellfun(@(x) x(:)',correlations_all,'UniformOutput',0);
correlations_all = cat(2,correlations_all{:});
correlations_all(correlations_all==1) = [];

fig_number = fig_number+1;
fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'ROI corr distribution','position',[680,558,946,413]);
subplot(1,2,1)
histogram(correlations_all,'BinWidth',0.1); xlim([-1 1]); hold on
plot([roi_corr_threshold(2), roi_corr_threshold(2)],ylim,'--r')
plot([roi_corr_threshold(3), roi_corr_threshold(3)],ylim,'--r')
plot([roi_corr_threshold(4), roi_corr_threshold(4)],ylim,'--r')
title([num2str(sum(correlations_all<roi_corr_threshold(2))) ...
    '/' num2str(sum(correlations_all<roi_corr_threshold(3))) ...
    '/' num2str(sum(correlations_all<roi_corr_threshold(4))) ...
    '/' num2str(length(correlations_all)) ])
box off
xlabel('minimum ROI corr per cell')
subplot(1,2,2)
histogram(correlations,'BinWidth',0.1); xlim([-1 1]); hold on
plot([roi_corr_threshold(2), roi_corr_threshold(2)],ylim,'--r')
plot([roi_corr_threshold(3), roi_corr_threshold(3)],ylim,'--r')
plot([roi_corr_threshold(4), roi_corr_threshold(4)],ylim,'--r')
title([num2str(sum(correlations<roi_corr_threshold(2))) ...
    '/' num2str(sum(correlations<roi_corr_threshold(3))) ...
    '/' num2str(sum(correlations<roi_corr_threshold(4))) ...
    '/' num2str(length(correlations)) ])
box off
xlabel('minimum ROI corr per cell')

%% figure 1d&h with roi corr cutoffs: stability with roi corr threshold
mouse_group = [1,2,3,4,5];
roi_corr_threshold = [-1 0.1 0.25 0.5];

for th = 1:length(roi_corr_threshold)
    
    fig_number = fig_number+1;
    fig_hand(fig_number) = figure;
    
    set(fig_hand(fig_number), 'Name', ['PO stability roicor over th ' char(th+64)], 'Position', [204,233,1473,598]);
    
    % PO at tp1 vs tp2 and tp20
    allPO_d1a = [];
    allPO_d1b = [];
    allPO_d2 = [];
    allissig_d2 = [];
    allPO_d20 = [];
    allissig_d20 = [];
    allcells_a = [];
    allcells_b = [];
    allmice_a = [];
    allmice_b = [];
    cell_counter = 0;
    for mouse_n = mouse_group
        dT1idx = data_stability{(mouse_n)}(1).comp_type_mat_days==1;
        dT20idx = data_stability{(mouse_n)}(1).comp_type_mat_days==20;
        dT1idx = tril(dT1idx,-1);
        dT20idx = tril(dT20idx,-1);
        dT1idx = find(dT1idx);
        dT20idx = find(dT20idx);
        
        for celln = 1:length(data_stability{(mouse_n)})
            if all(all(data_stability{mouse_n}(celln).ROI_correlations>roi_corr_threshold(th)))
                cell_counter = cell_counter+1;
                
                PO_d1_temp = data_stability{(mouse_n)}(celln).PO(1);
                PO_d2_temp = data_stability{(mouse_n)}(celln).POdif_all(dT1idx);
                PO_d20_temp = data_stability{(mouse_n)}(celln).POdif_all(dT20idx);
                
                allPO_d2 = cat(1,allPO_d2,...
                    wrapTo180((PO_d1_temp + PO_d2_temp).*2)./2);
                allissig_d2 = cat(1,allissig_d2,...
                    data_stability{(mouse_n)}(celln).POdif_issig(dT1idx));
                
                allPO_d20 = cat(1,allPO_d20,...
                    wrapTo180((PO_d1_temp + PO_d20_temp).*2)./2);
                allissig_d20 = cat(1,allissig_d20,...
                    data_stability{(mouse_n)}(celln).POdif_issig(dT20idx));
                
                allPO_d1a = cat(1,allPO_d1a, ones(length(dT1idx),1)*PO_d1_temp);
                allPO_d1b = cat(1,allPO_d1b, ones(length(dT20idx),1)*PO_d1_temp);
                allcells_a = cat(1,allcells_a, ones(length(dT1idx),1)*cell_counter);
                allcells_b = cat(1,allcells_b, ones(length(dT20idx),1)*cell_counter);
                allmice_a = cat(1,allmice_a, ones(length(dT1idx),1)*mouse_n); % to check from which mice the cells are comming from
                allmice_b = cat(1,allmice_b, ones(length(dT20idx),1)*mouse_n);
            end
            
        end
    end
    
    subplot(2,3,1)
    scatter(allPO_d1a(~allissig_d2),allPO_d2(~allissig_d2),'b','filled','MarkerFaceAlpha',0.2); hold on
    scatter(allPO_d1a(find(allissig_d2)),allPO_d2(find(allissig_d2)),'r','filled','MarkerFaceAlpha',0.2);
    set(gca,'TickDir','out')
    axis square
    xlim([-90 90]);ylim([-90 90])
    xticks(-90:30:90);yticks(-90:30:90)
    xlabel('PO day n')
    ylabel('PO day n+1')
    n_comp = sum(~isnan(allPO_d1a)&~isnan(allPO_d2));
    n_cells = length(unique(allcells_a(~isnan(allPO_d1a)&~isnan(allPO_d2))));
    n_mice = length(unique(allmice_a(~isnan(allPO_d1a)&~isnan(allPO_d2))));
    [r_cc, p_cc]=circ_corrcc_withUniCorrection(wrapToPi(deg2rad(allPO_d1a.*2)),wrapToPi(deg2rad(allPO_d2.*2)),0);
    title(['n: ' num2str(n_comp) '/' num2str(n_cells) '/' num2str(n_mice)...
        ', r: ' num2str(r_cc,3) '(' num2str(p_cc,3) ')'])
    
    subplot(2,3,2)
    scatter(allPO_d1b(~allissig_d20),allPO_d20(~allissig_d20),'b','filled','MarkerFaceAlpha',0.2); hold on
    scatter(allPO_d1b(find(allissig_d20)),allPO_d20(find(allissig_d20)),'r','filled','MarkerFaceAlpha',0.2)
    set(gca,'TickDir','out')
    axis square
    xlim([-90 90]);ylim([-90 90])
    xticks([-90:30:90]);yticks([-90:30:90])
    xlabel('PO day n')
    ylabel('PO day n+20')
    n_comp = sum(~isnan(allPO_d1b)&~isnan(allPO_d20));
    n_cells = length(unique(allcells_b(~isnan(allPO_d1b)&~isnan(allPO_d20))));
    n_mice = length(unique(allmice_b(~isnan(allPO_d1b)&~isnan(allPO_d20))));
    [r_cc, p_cc]=circ_corrcc_withUniCorrection(wrapToPi(deg2rad(allPO_d1b.*2)),wrapToPi(deg2rad(allPO_d20.*2)),0);
    title(['n: ' num2str(n_comp) '/' num2str(n_cells) '/' num2str(n_mice)...
        ', r: ' num2str(r_cc,3) '(' num2str(p_cc,3) ')'])
    
    [~,dT_bin_labels,~,~,~,~,out]=...
        eval_PO_stability(data_stability(mouse_group),'dT_binsize',2,...
        'roi_corr_threshold',roi_corr_threshold(th));
    
    all_days = cellfun(@(x) x(1).day, data_stability(mouse_group),'UniformOutput',0);
    max_days = max(cellfun(@(x) max(x),all_days));
    min_days = min(cellfun(@(x) max(x),all_days));
    over_n1_days = sort(cellfun(@(x) max(x),all_days));over_n1_days=over_n1_days(end-1);
    n_perint = cellfun(@(x) length(unique(x(7,:))),out(cellfun(@(x) ~isempty(x),out)));
    
    % dPO sig%
    subplot(2,3,3)
    all_sig = nan(length(mouse_group),length(out)-1);
    for mouse_n = mouse_group
        clear sigfrac_animal
        sigfrac_animal(1)=nan;
        for dT_n = 2:length(out)
            temp = out{dT_n};
            temp = temp(4,temp(7,:)==mouse_n); % 4th row is sig, 6th row is animal
            sigfrac_animal(dT_n) = sum(temp)*100/length(temp);
        end
        plot(sigfrac_animal,'color',[0.5 0.5 0.5],'LineWidth',1); hold on
        all_sig(mouse_n,1:length(sigfrac_animal)) = sigfrac_animal;
    end
    all_sig(:,~[n_perint>3])=nan;
    all_mean = nanmean(all_sig,1);
    all_std = nanstd(all_sig,1);
    all_n = sum(~isnan(all_sig),1);
    all_sem = all_std./sqrt(all_n);
    days_all = find(~isnan(all_mean));
    all_mean(isnan(all_mean)) = [];
    all_std(isnan(all_std)) = [];
    all_sem(isnan(all_sem)) = [];
    errorbar(days_all,all_mean,all_sem,'-k','LineWidth',2)
    ylim([0 60])
    box off
    ylabel('fraction of sig changes [%]')
    xlabel('time interval [days]')
    xticks(1:length(dT_bin_labels))
    xlim([1 length(dT_bin_labels)])
    xticklabels(dT_bin_labels)
    set(gca,'XTickLabelRotation',45)
    % stats (anova + dunnets)
    day_temp = reshape(repmat(1:size(all_sig(:,2:end-2),2),size(all_sig(:,2:end-2),1),1),[],1);
    all_sig_temp = reshape(all_sig(:,2:end-2),[],1);
    nans_idx = find(isnan(all_sig_temp));
    day_temp(nans_idx)= [];
    all_sig_temp(nans_idx)= [];
    [p,tbl,stats] = anova1(all_sig_temp,day_temp,'off');
    [results] = dunnett(stats);
    title(['anova1: p' num2str(p,3) ', F' num2str(tbl{2,5},4)])
    for i=1:length(results)
        if results(i)<0.05
            scatter(days_all(i),50,'*','MarkerEdgeColor','k')
        end
    end
    
    % |dPO| cumprob
    subplot(2,3,4)
    colors = cool(size(out,2));
    for dT_bin = 2:length(out)
        if n_perint(dT_bin)>3
            dPO_all{dT_bin} = out{dT_bin}(3,:);
            binsize = 10;
            line(dT_bin) = cdfplot(abs(dPO_all{dT_bin})); hold on
            set(line(dT_bin),'color',colors(dT_bin,:),'LineWidth',2)
        end
    end
    xlabel('|PO| change (°)'); ylabel('cumulative prob'); title('')
    axis square
    xlim([0 90])
    set(gca,'TickDir','out')
    cax = colorbar; colormap('cool');
    temp_ticks = [0:1/(length(dT_bin_labels)-1):1];
    cax.Ticks = temp_ticks(n_perint>3); cax.TickLabels = dT_bin_labels(n_perint>3);
    cax.Title.String = 'dT (days)'; cax.TickDirection = 'out';
    set(gca,'FontSize',11);
    
    % median |dPO|
    subplot(2,3,5); cla
    all_sig = nan(length(mouse_group),length(out)-1);
    clear dPO_med dPO_CI_temp dPO_sig_med dPO_sig_CI_temp
    for dT_bin = 1:length(out)
        if n_perint(dT_bin)>3
            dPO_all{dT_bin} = abs(out{dT_bin}(3,:));
            dPO_med(dT_bin) = nanmedian(dPO_all{dT_bin});
            nanmedian_boot = bootstrp(500,@nanmedian,dPO_all{dT_bin});
            dPO_CI_temp{dT_bin} = [prctile(nanmedian_boot,2.5); prctile(nanmedian_boot,97.5)];
            
            dPO_sig_all{dT_bin} = abs(out{dT_bin}(3,find(out{dT_bin}(4,:)==1)));
            dPO_sig_med(dT_bin) = nanmedian(dPO_sig_all{dT_bin});
            nanmedian_boot = bootstrp(500,@nanmedian,dPO_sig_all{dT_bin});
            dPO_sig_CI_temp{dT_bin} = [prctile(nanmedian_boot,2.5); prctile(nanmedian_boot,97.5)];
        else
            dPO_med(dT_bin) = nan;
            dPO_CI_temp{dT_bin} = [nan,nan]';
            dPO_sig_med(dT_bin) = nan;
            dPO_sig_CI_temp{dT_bin} = [nan,nan]';
        end
    end
    dPO_CI = cat(2,dPO_CI_temp{:})-dPO_med;
    errorbar(0:length(dPO_med)-1,dPO_med,dPO_CI(1,:),dPO_CI(2,:),'-k','LineWidth',2); hold on
    dPO_sig_CI = cat(2,dPO_sig_CI_temp{:})-dPO_sig_med;
    errorbar(0:length(dPO_sig_med)-1,dPO_sig_med,dPO_sig_CI(1,:),dPO_sig_CI(2,:),'--r','LineWidth',2)
    ylim([0 30])
    box off
    ylabel('med |dPO| [°]')
    xlabel('time interval [days]')
    xticks([0:length(dT_bin_labels)])
    xticklabels(dT_bin_labels)
    xlim([0 length(dT_bin_labels)-1])
    set(gca,'XTickLabelRotation',45)
    % stats
    clear combinedMatrix
    combinedMatrix = NaN(length(dPO_all), max(cellfun(@(x) length(x),dPO_all)));
    for i = 1:length(dPO_all)
        combinedMatrix(i, 1:length(dPO_all{i})) = dPO_all{i};
    end
    combinedMatrix(1,:)=[];
    day_temp = repmat(1:size(combinedMatrix,1),size(combinedMatrix,2),1)';
    day_temp = reshape(day_temp,[],1);
    all_dPO_temp = reshape(combinedMatrix,[],1);
    nans_idx = find(isnan(all_dPO_temp));
    day_temp(nans_idx)= [];
    all_dPO_temp(nans_idx)= [];
    [p_all,ANOVATAB_all] = kruskalwallis(all_dPO_temp',day_temp','off');
    p_ranksum = [nan,nan];
    for i = 3:length(dPO_all)
        p_ranksum(i) = ranksum(dPO_all{2}, dPO_all{i});
        if p_ranksum(i)<0.05/(length(dPO_all)-2)
            scatter(i-1,9,'*','MarkerEdgeColor','k')
        end
    end
    
    clear combinedMatrix
    combinedMatrix = NaN(length(dPO_sig_all), max(cellfun(@(x) length(x),dPO_sig_all)));
    for i = 1:length(dPO_sig_all)
        combinedMatrix(i, 1:length(dPO_sig_all{i})) = dPO_sig_all{i};
    end
    combinedMatrix(1,:)=[];
    day_temp = repmat(1:size(combinedMatrix,1),size(combinedMatrix,2),1)';
    day_temp = reshape(day_temp,[],1);
    all_dPO_temp = reshape(combinedMatrix,[],1);
    nans_idx = find(isnan(all_dPO_temp));
    day_temp(nans_idx)= [];
    all_dPO_temp(nans_idx)= [];
    [p_sig,ANOVATAB_sig] = kruskalwallis(all_dPO_temp',day_temp','off');
    p_ranksum = [nan,nan];
    for i = 3:length(dPO_all)
        p_ranksum(i) = ranksum(dPO_all{2}, dPO_all{i});
        if p_ranksum(i)<0.05/(length(dPO_all)-2)
            scatter(i-1,25,'*','MarkerEdgeColor','r')
        end
    end
    
    title({['Kruskal-Wallis on all: p' num2str(p_all,3) ', X' num2str(ANOVATAB_all{2,5},3)];
        ['Kruskal-Wallis on sig: ' num2str(p_sig,3) ', X' num2str(ANOVATAB_sig{2,5},3)]})
end

%% save figs

if savefigs
    mkdir([save_fig_loc '\rebuttal fig\'])
    for i=1:length(fig_hand)
        saveas(fig_hand(i),[save_fig_loc '\rebuttal fig\' fig_hand(i).Name],'svg')
        saveas(fig_hand(i),[save_fig_loc '\rebuttal fig\' fig_hand(i).Name],'fig');
    end
end