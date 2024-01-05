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

fig_number = 0;

%% 
mouse_group = [1,2,3,4,5,6];%
example_mouse = 5;

% cor decay
fig_number = fig_number+1;
fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'pairwise signal corr','position',[1060,740,343,248]);
[~,~,dSvals_PS_cor,PS_cor_FIT_dS] = calc_stability(Master_ROI_prop_stability(example_mouse),'PS_cor',1,fig_number); % pairwise singal corr (implement)
axis square

% cor decay
fig_number = fig_number+1;
fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'PS TC ER corr','position',[1031,528,875,420]);
[dTvals_PS_cor,PS_cor_FIT_dT,~,~,dTvals_PS_cor_allmice,PS_cor_FIT_dT_allmice,~,~,F_equation_PSC] = calc_stability(Master_ROI_prop_stability(mouse_group),'PS_cor',0); % pairwise singal corr (implement)
[dTvals_TC_cor,TC_cor_FIT_dT,~,~,dTvals_TC_cor_allmice,TC_cor_FIT_dT_allmice,~,~,F_equation_TC] = calc_stability(Master_ROI_prop_stability(mouse_group),'TC_cor',0); % tuning curve corr (check code)
[dTvals_ER_cor,ER_cor_FIT_dT,~,~,dTvals_ER_cor_allmice,ER_cor_FIT_dT_allmice,~,~,F_equation_ER] = calc_stability(Master_ROI_prop_stability(mouse_group),'ER_cor',0); % ensemble rate corr  (check code)
% [dTvals_PV_cor,PV_cor_FIT_dT,dSvals_PV_cor,PV_cor_FIT_dS] = calc_stability(Master_ROI_prop_stability(mouse_group),'PV_cor',0); % ensemble rate corr  (check code)

subplot(1,3,1); cla
for i = 1:length(dTvals_PS_cor_allmice)
   plot(dTvals_PS_cor_allmice{i},PS_cor_FIT_dT_allmice{i},'color',[0.5 0.5 0.5]); hold on
end
plot(dTvals_PS_cor,PS_cor_FIT_dT,'r'); hold on
box off
ylim([0, 1]); xlim([0, 30]);
xlabel('time interval [days]')
ylabel('PS cor.')
axis square
title(F_equation_PSC,'Interpreter','none','FontSize',10)

subplot(1,3,2); cla
for i = 1:length(dTvals_PS_cor_allmice)
   plot(dTvals_TC_cor_allmice{i},TC_cor_FIT_dT_allmice{i},'color',[0.5 0.5 0.5]); hold on
end
plot(dTvals_TC_cor,TC_cor_FIT_dT,'r'); hold on
box off
ylim([0, 1]); xlim([0, 30]);
xlabel('time interval [days]')
ylabel('TC cor.')
axis square
title(F_equation_TC,'Interpreter','none','FontSize',10)

subplot(1,3,3); cla
for i = 1:length(dTvals_PS_cor_allmice)
   plot(dTvals_ER_cor_allmice{i},ER_cor_FIT_dT_allmice{i},'color',[0.5 0.5 0.5]); hold on
end
plot(dTvals_ER_cor,ER_cor_FIT_dT,'r'); hold on
box off
ylim([0, 1]); xlim([0, 30]);
xlabel('time interval [days]')
ylabel('ER cor.')
axis square
title(F_equation_ER,'Interpreter','none','FontSize',10)

%% figure 1d-h: stability
mouse_group = [1,2,3,4,5,6];

fig_number = fig_number+1;
fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'PO stability quantification', 'Position', [204,233,1473,598]);

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
    dT1idx = Master_ROI_prop_stability{(mouse_n)}(1).comp_type_mat_days==1;
    dT20idx = Master_ROI_prop_stability{(mouse_n)}(1).comp_type_mat_days==20;
    dT1idx = tril(dT1idx,-1);
    dT20idx = tril(dT20idx,-1);
    dT1idx = find(dT1idx);
    dT20idx = find(dT20idx);
    
    for celln = 1:length(Master_ROI_prop_stability{(mouse_n)})
        cell_counter = cell_counter+1;
        
        PO_d1_temp = Master_ROI_prop_stability{(mouse_n)}(celln).PO(1);
        PO_d2_temp = Master_ROI_prop_stability{(mouse_n)}(celln).POdif_all(dT1idx);
        PO_d20_temp = Master_ROI_prop_stability{(mouse_n)}(celln).POdif_all(dT20idx);
        
        allPO_d2 = cat(1,allPO_d2,...
            wrapTo180((PO_d1_temp + PO_d2_temp).*2)./2);
        allissig_d2 = cat(1,allissig_d2,...
            Master_ROI_prop_stability{(mouse_n)}(celln).POdif_issig(dT1idx));
        
        allPO_d20 = cat(1,allPO_d20,...
            wrapTo180((PO_d1_temp + PO_d20_temp).*2)./2);
        allissig_d20 = cat(1,allissig_d20,...
            Master_ROI_prop_stability{(mouse_n)}(celln).POdif_issig(dT20idx));
        
        allPO_d1a = cat(1,allPO_d1a, ones(length(dT1idx),1)*PO_d1_temp);
        allPO_d1b = cat(1,allPO_d1b, ones(length(dT20idx),1)*PO_d1_temp);
        allcells_a = cat(1,allcells_a, ones(length(dT1idx),1)*cell_counter);
        allcells_b = cat(1,allcells_b, ones(length(dT20idx),1)*cell_counter);
        allmice_a = cat(1,allmice_a, ones(length(dT1idx),1)*mouse_n); % to check from which mice the cells are comming from
        allmice_b = cat(1,allmice_b, ones(length(dT20idx),1)*mouse_n);
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
    eval_PO_stability(Master_ROI_prop_stability(mouse_group),'dT_binsize',2);

all_days = cellfun(@(x) x(1).day, Master_ROI_prop_stability(mouse_group),'UniformOutput',0);
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

%% sup fig: Running controls
fig_number = 0; clear fig_hand; %close all;% start figure count for handles
mouse_group = [1,2,3,4,5];% 6 does not have running extracted?

% Approach 1
clear PO_dif_deg PO_issig PO_still PO_run
clear vect_PO_still vect_POCI_dist_still vect_POCI_still
clear vect_PO_run vect_POCI_dist_run vect_POCI_run PO_issig
running_threshold = 2;
counter_c=0;
neuron_counter = 0;
for mouse_n = 1:5
    for cell_n = 1:length(Master_ROI_prop_stability{mouse_n})
        included_neuron = 0;
        vis_resp = Master_ROI_prop_stability{mouse_n}(cell_n).visually_responsive(:);
        tuned = Master_ROI_prop_stability{mouse_n}(cell_n).Stat_tuned(:);
        
        for tp_n=find(vis_resp&tuned,1,'first')
            if vis_resp(tp_n) & tuned(tp_n)
                running_temp = squeeze(mean(Master_ROI_prop_stability{mouse_n}(1).all_Running_traces_stim{tp_n},1));
                activ_temp = squeeze(mean(Master_ROI_prop_stability{mouse_n}(cell_n).all_traces_stim{tp_n},1));
                running_idx = running_temp>running_threshold;
                activ_run = activ_temp';
                activ_stil = activ_temp';
                activ_stil(running_idx) = nan;
                activ_run(~running_idx) = nan;
                for stim_n = 1:12
                    [~,I] = sort(isnan(activ_stil(stim_n,:)));
                    activ_stil(stim_n,:) = activ_stil(stim_n,I);
                    [~,I] = sort(isnan(activ_run(stim_n,:)));
                    activ_run(stim_n,:) = activ_run(stim_n,I);
                end
                ntrials_still=min(sum(~isnan(activ_stil),2));
                ntrials_run=min(sum(~isnan(activ_run),2));
                
                ratio_trials=min([ntrials_still,ntrials_run])/max([ntrials_still,ntrials_run]);
                
                if ratio_trials>0.20 && min([ntrials_still,ntrials_run])>6
                    included_neuron=1;
                    counter_c=counter_c+1;
                    [vect_PO_temp1, vect_POCI_dist_temp1, vect_POCI_temp1] = ...
                        vect_PO_boot(activ_stil,0,0);
                    [vect_PO_temp2, vect_POCI_dist_temp2, vect_POCI_temp2] = ...
                        vect_PO_boot(activ_run,0,0);
                    
                    PO_rad_still = vect_PO_temp1;
                    PO_still(counter_c) = wrapTo180(rad2deg(vect_PO_temp1'));
                    POCI_still{counter_c}=wrapTo180(rad2deg(vect_POCI_temp1'));
                    PO_CI_lb_still = wrapToPi(deg2rad(POCI_still{counter_c}(1)));
                    PO_CI_ub_still = wrapToPi(deg2rad(POCI_still{counter_c}(2)));
                    
                    PO_rad_run = vect_PO_temp2;
                    PO_run(counter_c) = wrapTo180(rad2deg(vect_PO_temp2'));
                    POCI_run{counter_c}=wrapTo180(rad2deg(vect_POCI_temp2'));
                    PO_CI_lb_run = wrapToPi(deg2rad(POCI_run{counter_c}(1)));
                    PO_CI_ub_run = wrapToPi(deg2rad(POCI_run{counter_c}(2)));
                    
                    PO_dif_deg(counter_c) = 2*wrapTo180(rad2deg(wrapToPi(circ_dist(PO_rad_still.*2,PO_rad_run.*2)./2))./2);
                    
                    PO_issig(counter_c) = ~all(circ_dist([PO_rad_still,PO_rad_still].*2,[PO_CI_ub_run PO_CI_lb_run].*2)>0==[0 1]) && ...
                        ~all(circ_dist([PO_rad_run,PO_rad_run].*2,[PO_CI_ub_still PO_CI_lb_still].*2)>0==[0 1]);
                    
                end
            end
        end
        if included_neuron
            neuron_counter=neuron_counter+1;
        end
    end
end

still_vs_run = cat(1,PO_still,PO_run,PO_dif_deg,PO_issig);
sig_frac_still_vs_run = sum(still_vs_run(4,:)==1)/size(still_vs_run,2);

fig_number = fig_number+1;
fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'running vs still', 'Position', [1241,558,560, 420]);
subplot(1,1,1); cla; hold on
[~,dT_bin_labels,~,~,~,dPOall_perint_binned,out]=...
    eval_PO_stability(Master_ROI_prop_stability(mouse_group),'dT_binsize',2); % include all mice for across day comparison
short_TPn = 2;
long_TPn = 11;
int_short = cdfplot(abs(wrapTo180(out{short_TPn}(3,:))));
int_long = cdfplot(abs(wrapTo180(out{long_TPn}(3,:))));
int_runvsStill = cdfplot(abs(wrapTo180(still_vs_run(3,:))));
set(int_short,'Color','c')
set(int_long,'Color','m')
set(int_runvsStill,'Color','k')
axis square
xlim([0 90])
xlabel('|dPO|')
ylabel('cumulative prob.')
title([])
[p_short,~,tbl_short] = ranksum(int_short.XData,int_runvsStill.XData);
[p_long,~,tbl_long] = ranksum(int_long.XData,int_runvsStill.XData);
legend({[num2str(length(out{short_TPn}(3,:))) 'n, '],...
    [num2str(length(out{long_TPn}(3,:))) 'n, '],...
    [num2str(length(still_vs_run(3,:))) 'n, ' num2str(neuron_counter) 'neurons']},...
    'Location','SouthEast')

title({[dT_bin_labels{short_TPn} ': p' num2str(p_short,3) ', U' num2str(tbl_short.ranksum,3)];
    [dT_bin_labels{long_TPn} ': p' num2str(p_long,3) ', U' num2str(tbl_long.ranksum,3)]});

% Approach 2
[~,dT_bin_labels,~,~,~,~,out]=...
    eval_PO_stability(Master_ROI_prop_stability(mouse_group),'dT_binsize',2,'max_interval',20);

clear out_nonsig dPOabs  RMI dSpeed speedmod AMI dPupil dPupil pupilmod issig dT_idx
for i = 2:length(dT_bin_labels)
    dPOabs{i} = abs(out{i}(3,:))';
    speedmod{i} = abs(out{i}(5,:).*out{i}(6,:))';
    pupilmod{i} = abs(out{i}(8,:).*out{i}(9,:))';
    issig{i} = out{i}(4,:)';
    dT_idx{i} = repmat((i-1),size(out{i},2),1);
    
    % remove all comparisons where nans are present in any of the metrics
    idx_nan = any(isnan(cat(2,dPOabs{i},speedmod{i},pupilmod{i})),2);
    speedmod{i}(idx_nan) = [];
    pupilmod{i}(idx_nan) = [];
    issig{i}(idx_nan)    = [];
    dT_idx{i}(idx_nan)   = [];
end
dT = cellfun(@(x) x(1),dT_idx(2:end));

speedmod_all = cat(1,speedmod{:});
speedmod_all = speedmod_all(~isnan(speedmod_all));
SMI_thr_run = prctile(speedmod_all(~isnan(speedmod_all)),50); % same as median btw
pupilmod_all = cat(1,pupilmod{:});
pupilmod_all = pupilmod_all(~isnan(pupilmod_all));
PMI_thr_pupil = prctile(pupilmod_all,50); % same as median btw

fig_number = fig_number+1;
fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'behavioural modulation index based exclusion', 'Position', [28,421,1095,540]);

subplot(2,4,1:2); cla
histogram(speedmod_all,'BinWidth',SMI_thr_run); box off
xlim([0 inf])
vline(SMI_thr_run,'--r')
xlabel('Running modulation')
ylabel('count')

subplot(2,4,3:4); cla
histogram(pupilmod_all,'BinWidth',PMI_thr_pupil); box off
xlim([0 inf])
vline(PMI_thr_pupil,'--b')
xlabel('Arousal modulation')
ylabel('count')

subplot(2,4,5:6); cla; hold on
perSig = @(x) 100*sum(x)/length(x);
perSig_all = cellfun(@(x) perSig(x),issig(2:end));
perSig_norun = cellfun(@(x,y) perSig(x(find(y<SMI_thr_run))),issig(2:end),speedmod(2:end));
perSig_nopupil = cellfun(@(x,y) perSig(x(find(y<PMI_thr_pupil))),issig(2:end),pupilmod(2:end));
boot95_l = @(x) prctile(bootstrp(500,@(x) 100*sum(x)/length(x),x),2.5);
boot95_u = @(x) prctile(bootstrp(500,@(x) 100*sum(x)/length(x),x),97.5);

l95_all = cellfun(@(x) boot95_l(x),issig(2:end));
l95_norun = cellfun(@(x,y) boot95_l(x(find(y<SMI_thr_run))),issig(2:end),speedmod(2:end));
l95_nopupil = cellfun(@(x,y) boot95_l(x(find(y<PMI_thr_pupil))),issig(2:end),pupilmod(2:end));
u95_all = cellfun(@(x) boot95_u(x),issig(2:end));
u95_norun = cellfun(@(x,y) boot95_u(x(find(y<SMI_thr_run))),issig(2:end),speedmod(2:end));
u95_nopupil = cellfun(@(x,y) boot95_u(x(find(y<PMI_thr_pupil))),issig(2:end),pupilmod(2:end));
l95_all = l95_all-perSig_all;
l95_norun = l95_norun-perSig_norun;
l95_nopupil = l95_nopupil-perSig_nopupil;
u95_all = u95_all-perSig_all;
u95_norun = u95_norun-perSig_norun;
u95_nopupil = u95_nopupil-perSig_nopupil;

errorbar(1:length(perSig_all),perSig_all,l95_all,u95_all,[],[],'k')
errorbar(1:length(perSig_norun)-0.1,perSig_norun,l95_norun,u95_norun,[],[],'r')
errorbar(1:length(perSig_nopupil)+0.1,perSig_nopupil,l95_nopupil,u95_nopupil,[],[],'b')
xticks(1:1:length(dT_bin_labels))
xticklabels(dT_bin_labels(1:1:end))
xtickangle(45)
xlabel('time interval [days]')
ylabel('%sig changes')

subplot(2,4,7:8); cla; hold on
median_all = cellfun(@(x) nanmedian(x),dPOabs(2:end));
median_norun = cellfun(@(x,y) nanmedian(x(find(abs(y)<SMI_thr_run))),dPOabs(2:end),speedmod(2:end));
median_nopupil = cellfun(@(x,y) nanmedian(x(find(abs(y)<PMI_thr_pupil))),dPOabs(2:end),pupilmod(2:end));
boot95_l = @(x) prctile(bootstrp(500,@nanmedian,x),2.5);
boot95_u = @(x) prctile(bootstrp(500,@nanmedian,x),97.5);

l95_all = cellfun(@(x) boot95_l(x),dPOabs(2:end));
l95_norun = cellfun(@(x,y) boot95_l(x(find(abs(y)<SMI_thr_run))),dPOabs(2:end),speedmod(2:end));
l95_nopupil = cellfun(@(x,y) boot95_l(x(find(abs(y)<PMI_thr_pupil))),dPOabs(2:end),pupilmod(2:end));
u95_all = cellfun(@(x) boot95_u(x),dPOabs(2:end));
u95_norun = cellfun(@(x,y) boot95_u(x(find(abs(y)<SMI_thr_run))),dPOabs(2:end),speedmod(2:end));
u95_nopupil = cellfun(@(x,y) boot95_u(x(find(abs(y)<PMI_thr_pupil))),dPOabs(2:end),pupilmod(2:end));
l95_all = l95_all-median_all;
l95_norun = l95_norun-median_norun;
l95_nopupil = l95_nopupil-median_nopupil;
u95_all = u95_all-median_all;
u95_norun = u95_norun-median_norun;
u95_nopupil = u95_nopupil-median_nopupil;

errorbar(1:length(median_all),median_all,l95_all,u95_all,[],[],'k')
errorbar(1:length(median_norun)-0.1,median_norun,l95_norun,u95_norun,[],[],'r')
errorbar(1:length(median_nopupil)+0.1,median_nopupil,l95_nopupil,u95_nopupil,[],[],'b')
xticks(1:1:length(dT_bin_labels))
xticklabels(dT_bin_labels(1:1:end))
xtickangle(45)
ylim([0 8])
xlabel('time interval [days]')
ylabel('drift magnitude')
