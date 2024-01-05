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

%% stripe rearing population effect
mouse_group = [9 10 11 14 16 17 18];

within_bin_ratios = 0;
cross_animal = 1;

mintp = 3;
SR_idx = 3;

time_point_labels = {'BL1', 'BL2', 'SR1'};
all_cell_filters = {'resp_anytime',['resp_and_tuned_from_BL2toSR' num2str(SR_idx-2)],['resp_and_tuned_at_BL2xorSR' num2str(SR_idx-2)]};

fig_number = fig_number+1;
fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'SR28d population effect', 'Position', [1,31,1920,973]);

clear PO_bin_frac_sig_ofTuned
for i=1:length(all_cell_filters)
    figure(fig_number) % select correct figure
    apply_cell_filter = all_cell_filters{i};
    clear PO_to_PrefO_count_per_mouse PO_to_PrefO_count PO_bin_centers
    
    for mouse_n = mouse_group
        if ~isempty(Master_ROI_prop_StripeRearing{mouse_n})
            [PO_to_PrefO_count_per_mouse_bin90{mouse_n},~] = StripeRearing_stability(Master_ROI_prop_StripeRearing(mouse_n),time_point_labels(1:mintp),...
                'comparison_type','BL2-All','Binning_interval' , 90,...
                'cell_filter',apply_cell_filter,'use_only_sig_changes',0);
            [PO_to_PrefO_count_per_mouse_bin30{mouse_n},~] = StripeRearing_stability(data_orientation_deprivation(mouse_n),time_point_labels(1:mintp),...
                'comparison_type','BL2-All','Binning_interval' , 30,...
                'cell_filter',apply_cell_filter,'use_only_sig_changes',0);
            
            if i == 1 % get
                [~,~,PO_bin_frac_sig_ofTuned{mouse_n}] = StripeRearing_stability(data_orientation_deprivation(mouse_n),time_point_labels(1:mintp),...
                    'comparison_type','BL2-All','Binning_interval' , 90,...
                    'cell_filter',apply_cell_filter,'use_only_sig_changes',0);
            end
        end
    end
    
    idx_BL2 = find(strcmp(time_point_labels,'BL2'));
    change_PermOri_per_mouse = cellfun(@(x) (x(SR_idx,2) - x(idx_BL2,2)),PO_to_PrefO_count_per_mouse_bin90(mouse_group));
    change_OrthOri_per_mouse = cellfun(@(x) (x(SR_idx,1) - x(idx_BL2,1)),PO_to_PrefO_count_per_mouse_bin90(mouse_group));
    
    cellcount_BL2=num2cell(cellfun(@(x)  x(idx_BL2,1) + x(idx_BL2,2), PO_to_PrefO_count_per_mouse_bin90(mouse_group)));
    percent_change_PermOri_per_mouse = cellfun(@(x,y) (x(SR_idx,2) - x(idx_BL2,2))/(y),PO_to_PrefO_count_per_mouse_bin90(mouse_group),cellcount_BL2);
    percent_change_OrthOri_per_mouse = cellfun(@(x,y) (x(SR_idx,1) - x(idx_BL2,1))/(y),PO_to_PrefO_count_per_mouse_bin90(mouse_group),cellcount_BL2);
    
    figure(fig_number) % select correct figure
    subplot(2,3,i); cla
    colors = [0.5 0.5 0.5;0 0 0; 1 0 0];
    for tp = [2 SR_idx] %1:mintp
        all_dist = cellfun(@(x) x(tp,:)/sum(x(2,:)), PO_to_PrefO_count_per_mouse_bin30(mouse_group),'UniformOutput', false);
        
        mean_vals = mean(cat(1,all_dist{:}).*100,1);
        errer_val = std(cat(1,all_dist{:}).*100)./sqrt(length(mouse_group));
        x_val = -90:30:60;
        errorbar(x_val,mean_vals,errer_val,'color',colors(tp,:),'LineWidth',2); hold on
        plot(x_val,mean_vals,'color',colors(tp,:),'LineWidth',2)
        
        set(gca,'TickDir','out'); box off; set(gca,'TickDir','out')
        xticks(-90:30:60); xlim([-100 70]);ylim([0 inf])
        xlabel([{'pref. Ori. (PO) diff.'};{'to Permitted Ori'}])
        if i==1
            ylabel({'% relative to' ;'total count at BL2'})
        end
    end
    pbaspect([1.5,0.75,0.5])
    ylim([0 30])
    if i==3
        ylim([0 50])
    end
    title(apply_cell_filter,'Interpreter','none')
    
    figure(fig_number);  % select correct figure
    ax2(i) = subplot(2,3,i+3);
    for mouse_n = 1:length(mouse_group)
        plot([2,1],100.*[percent_change_OrthOri_per_mouse(mouse_n), percent_change_PermOri_per_mouse(mouse_n)],'color',[0.5 0.5 0.5]); hold on
    end
    plot([2,1],100.*[mean(percent_change_OrthOri_per_mouse), mean(percent_change_PermOri_per_mouse)],'k','LineWidth',5); hold on
    xticks([1,2]);xticklabels({'0' '90'})
    xlim([0.5 2.5]);
    if i~=3
        ylim([-40 40])
    end
    plot(xlim,[0 0],'--k')
    box off
    set(gca,'TickDir','out')
    xlabel([{'pref. Ori. (PO) diff.'};{'to Permitted Ori'}])
    ylabel({'% change in cell count'; 'from BL2 to SR1'})
    ylim([-1*max(abs(ylim)) max(abs(ylim))])
    [~,p_90,~,STATS_90]=ttest(percent_change_OrthOri_per_mouse);
    [~,p_0,~,STATS_0]=ttest(percent_change_PermOri_per_mouse);
    [~,p_paired,~,STATS_paired]=ttest(percent_change_PermOri_per_mouse-percent_change_OrthOri_per_mouse);
    title({['t(' num2str(STATS_paired.df) ')' num2str(STATS_paired.tstat,3) ',p=' num2str(p_paired,'%.3f')];...
        ['t(' num2str(STATS_0.df) ')' num2str(STATS_0.tstat,3) ',p=' num2str(p_0,'%.3f') ...
        '; t(' num2str(STATS_90.df) ')' num2str(STATS_90.tstat,3) ',p=' num2str(p_90,'%.3f')]...
        },'Interpreter','none')
    pbaspect([0.75,1.3,0.5])
    ax = gca;
    ax.TitleFontSizeMultiplier = 1;
end
linkaxes(ax2(1:2))

all_POs = cellfun(@(x) [x.PO],data_orientation_deprivation([mouse_group]),'UniformOutput', false);
all_POs = cat(2,all_POs{:});
clear all_vt all_vr
for mouse_n = mouse_group
    temp = arrayfun(@(x) x.visually_responsive',data_orientation_deprivation{mouse_n},'UniformOutput', 0);
    all_vr{mouse_n} = cat(2,temp{:});
    temp = arrayfun(@(x) x.Stat_tuned',data_orientation_deprivation{mouse_n},'UniformOutput', 0);
    all_vt{mouse_n} = cat(1,temp{:})';
end
all_vr = logical(cat(2,all_vr{:}));
all_vt = cat(2,all_vt{:});
all_POs(~all_vr | ~all_vt)=nan;

all_permittedOri = cellfun(@(x) [x.permitted_ori],data_orientation_deprivation(mouse_group),'UniformOutput', false);
all_POissig = cellfun(@(x) arrayfun(@(y) y.POdif_issig(2,SR_idx),x,'UniformOutput', false), data_orientation_deprivation(mouse_group),'UniformOutput', false);
all_permittedOri = cat(2,all_permittedOri{:});
all_POissig = cat(2,all_POissig{:}); all_POissig = cat(2,all_POissig{:});
all_POs_relToPerm = wrapTo180([all_POs-all_permittedOri].*2)./2;

%% scatter plots

% 28d Stripe rearing
mouse_group = [9 10 11 14 16 17 18];


time_point_labels = {'BL1', 'BL2', 'SR1'};
[~, ~, ~, comp_type_label,~,~,~,~,data_out] = StripeRearing_stability(data_orientation_deprivation(mouse_group),time_point_labels,...
    'comparison_type','BL2-All','Binning_interval' , 30,...
    'cell_filter','resp_anytime','use_only_sig_changes',0);
for IT_bin = 1:length(data_out)
    temp = isnan(data_out{IT_bin});
    hasnan = any(temp(1:2,:),1);
    data_out{IT_bin}(:,hasnan)=[];
end
data_out = data_out{strcmp(comp_type_label,'BL2-SR1')};


issig=logical(data_out(4,:));
dPOperm = abs(data_out(1,:))-abs(data_out(2,:));

fig_number = fig_number+1;
fig_hand(fig_number) = figure;

set(fig_hand(fig_number), 'Name', 'conv_n_drift_28d_SR', 'Position', ...
    [900,661,887,333],'Renderer','painters')

subplot(1,3,1)
p1 = scatter(data_out(1,~issig),data_out(2,~issig),'b','filled'); hold on
p2 = scatter(data_out(1,issig),data_out(2,issig),'r','filled');
alpha(p1,0.2);
alpha(p2,0.2);
xlim([-90,90]);
ylim([-90,90]);
plot([-90 90],[-90 90],'--k')
xticks(-90:30:90);
yticks(-90:30:90);
axis square
xlabel('PO1-perm')
ylabel('PO2-perm')
notnan = find(~isnan(data_out(1,:))&~isnan(data_out(2,:)));
[r,p]=circ_corrcc_withUniCorrection(wrapToPi(deg2rad(data_out(1,notnan)').*2),wrapToPi(deg2rad(data_out(2,notnan)').*2),1);
title({['n=' num2str(size(data_out,2)) ', r=' num2str(r,3) ', p=' num2str(p,3)]},'Interpreter','none')

subplot(1,3,2); cla
p1 = patch([0 90 90 0],[-90 -90 0 -90],[0.8 0.8 0.8],'EdgeColor','none');
p2 = patch([0 90 0 0],[0 90 90 0],[0.8 0.8 0.8],'EdgeColor','none'); hold on
alpha(p1,0.2);
alpha(p2,0.2);
p1 = scatter(abs(data_out(1,~issig)),dPOperm(~issig),'b','filled'); hold on
p2 = scatter(abs(data_out(1,issig)),dPOperm(issig),'r','filled');
alpha(p1,0.2);
alpha(p2,0.2);
xlim([0,90]);
ylim([-90,90]);
plot([0 90],[0 0],'--k')
xticks([0:30:90]);
yticks([-90:30:90]);
axis square
xlabel('|PO1-perm|')
ylabel('d|PO-perm|')

subplot(1,3,3)
bin_width = 30;
bin_centers = [0:bin_width:90-(bin_width)]+bin_width/2;
p1 = scatter(abs(data_out(1,~issig)),abs(data_out(3,~issig)),'b','filled'); hold on
p2 = scatter(abs(data_out(1,issig)),abs(data_out(3,issig)),'r','filled');
alpha(p1,0.2);
alpha(p2,0.2);
xlim([0,90]);
ylim([0,90]);
xticks([0:30:90]);
yticks([0:30:90]);
axis square
xlabel('|PO1-perm|')
ylabel('|dPO|')
[r,p]=corr(abs(data_out(1,:))',abs(data_out(3,:))','Type','Spearman');
title({[num2str(sum(issig)) '/' num2str(length(issig)) ' (sig/all)'];
    ['r=' num2str(r,2) ', p=' num2str(p,2)]},'FontSize',10)


%% drift and convergence as a function of time for all conditions
% stability data
mouse_group = [1,2,3,4,5,6];%
all_days = cellfun(@(x) x(1).day, data_stability(mouse_group),'UniformOutput',0);

[~,dT_bin_labels,~,~,~,dPOall_perint,out] = ...
    eval_PO_stability(data_stability(mouse_group),'max_interval',28);
n_perint = cellfun(@(x) length(unique(x(7,:))),out(cellfun(@(x) ~isempty(x),out)));

clear dPOsig_AUC dPOsig_AUC_CI dPOall_AUC dPOall_AUC_CI dPOsig_med dPOsig_med_CI dPOall_med dPOall_med_CI
cdf_trapz=@(x) trapz(sort(abs(x)),0:1/(length(x)-1):1);
for dTbin = 2:length(dT_bin_labels)
    if ~isempty(dPOall_perint{dTbin})
        temp = dPOall_perint{dTbin}; temp(isnan(dPOall_perint{dTbin}))=[];
        dPOall_AUC(dTbin) = cdf_trapz(temp);
        dPOall_AUC_CI(dTbin,:) = bootci(500,@(x) cdf_trapz(x),temp);
        dPOall_med(dTbin) = nanmedian(abs(temp));
        dPOall_med_CI(dTbin,:) = bootci(500,@(x) nanmedian(abs(x)),temp);
    else
        dPOall_AUC(dTbin) = nan;
        dPOall_AUC_CI(dTbin,:)  = [nan nan];
        dPOall_med(dTbin) = nan;
        dPOall_med_CI(dTbin,:)  = [nan nan];
    end
end
fig_number = fig_number+1;
fig_hand(fig_number) = figure(206); cla
set(fig_hand(fig_number), 'Name', 'med of dPO', 'Position', [12,376,570,260]);
temp_nonan = find(~isnan(dPOall_med(2:end)))+1;
errorbar(find(n_perint>3)-1,dPOall_med(n_perint>3),...
    dPOall_med_CI(n_perint>3,1)-dPOall_med(n_perint>3)',...
    dPOall_med_CI(n_perint>3,2)-dPOall_med(n_perint>3)','-k'); hold on
scatter(find(n_perint>3)-1,dPOall_med(n_perint>3),'k','filled'); hold on

xlabel('dT (days)'); ylabel('median |dPO| [°]')
set(gca,'TickDir','out')
ylim([0 inf])
box off

% 7 day stripe rearing
mouse_group = [1 2 3 4 5 6 7 8];
maxtp = 6;

time_point_labels = {'BL1', 'BL2', 'SR1', 'SR2', 'SR3', 'SR4'};
[~, ~, ~, comp_type_label,~,~,~,dPOall_perint] = StripeRearing_stability(data_orientation_deprivation(mouse_group),time_point_labels,...
    'comparison_type','BL2-All','Binning_interval' , 30,...
    'cell_filter','resp_anytime','use_only_sig_changes',0);

clear dPOsig_AUC dPOsig_AUC_CI dPOall_AUC dPOall_AUC_CI dPOsig_med dPOsig_med_CI dPOall_med dPOall_med_CI
cdf_trapz=@(x) trapz(sort(abs(x)),0:1/(length(x)-1):1);
for dTbin = 2:length(comp_type_label)
    if ~isempty(dPOall_perint{dTbin})
        temp = dPOall_perint{dTbin}; temp(isnan(dPOall_perint{dTbin}))=[];
        dPOall_AUC(dTbin) = cdf_trapz(temp);
        dPOall_AUC_CI(dTbin,:) = bootci(500,@(x) cdf_trapz(x),temp);
        dPOall_med(dTbin) = nanmedian(abs(temp));
        dPOall_med_CI(dTbin,:) = bootci(500,@(x) nanmedian(abs(x)),temp);
    else
        dPOall_AUC(dTbin) = nan;
        dPOall_AUC_CI(dTbin,:)  = [nan nan];
        dPOall_med(dTbin) = nan;
        dPOall_med_CI(dTbin,:)  = [nan nan];
    end
end
figure(206);
scatter(7:7:7*(maxtp-2),dPOall_med(2:end),'r','filled'); hold on
errorbar(7:7:7*(maxtp-2),dPOall_med(2:end),dPOall_med_CI(2:end,1)-dPOall_med(2:end)',dPOall_med_CI(2:end,2)-dPOall_med(2:end)','-r')

% 28d Stripe rearing
mouse_group = [9 10 11 14 16 17 18];

time_point_labels = {'BL1', 'BL2', 'SR1'};
[~, ~, ~, comp_type_label,~,~,~,dPOall_perint] = StripeRearing_stability(data_orientation_deprivation(mouse_group),time_point_labels,...
    'comparison_type','BL2-All','Binning_interval', 30,...
    'cell_filter','resp_anytime','use_only_sig_changes',0);

clear dPOsig_AUC dPOsig_AUC_CI dPOall_AUC dPOall_AUC_CI dPOsig_med dPOsig_med_CI dPOall_med dPOall_med_CI
cdf_trapz=@(x) trapz(sort(abs(x)),0:1/(length(x)-1):1);
for dTbin = 2:length(comp_type_label)
    if ~isempty(dPOall_perint{dTbin})
        temp = dPOall_perint{dTbin}; temp(isnan(dPOall_perint{dTbin}))=[];
        dPOall_AUC(dTbin) = cdf_trapz(temp);
        dPOall_AUC_CI(dTbin,:) = bootci(500,@(x) cdf_trapz(x),temp);
        dPOall_med(dTbin) = nanmedian(abs(temp));
        dPOall_med_CI(dTbin,:) = bootci(500,@(x) nanmedian(abs(x)),temp);
    else
        dPOall_AUC(dTbin) = nan;
        dPOall_AUC_CI(dTbin,:)  = [nan nan];
        dPOall_med(dTbin) = nan;
        dPOall_med_CI(dTbin,:)  = [nan nan];
    end
end
figure(206);
scatter(28,dPOall_med(2),'g','filled'); hold on
plot(28,dPOall_med(2),'--g'); hold on
errorbar(28,dPOall_med(2),dPOall_med_CI(2,1)-dPOall_med(2)',dPOall_med_CI(2,2)-dPOall_med(2)','.g')

% now look at convergence
mouse_group = [9 10 11 14 16 17 18];
mean_permitted = mean(cellfun(@(x) x(1).permitted_ori,data_orientation_deprivation(mouse_group)));

% stability data
mouse_group = [1,2,3,4,5,6];%
all_days = cellfun(@(x) x(1).day, data_stability(mouse_group),'UniformOutput',0);
min_days = min(cellfun(@(x) max(x),all_days));
max_days = max(cellfun(@(x) max(x),all_days));
[~,dT_bin_labels,conTOdiv,conTOdiv_CI] = ...
    eval_PO_stability(data_stability(mouse_group),'PO_relative',mean_permitted,...
    'use_only_sig_changes',0,'max_interval',max_days-1);

fig_number = fig_number+1;
fig_hand(fig_number) = figure(203);
idx_short = find(strcmp(dT_bin_labels,'20'));
set(fig_hand(fig_number), 'Name', 'convergence all', 'Position', [9,34,570,254])
errorbar(find(n_perint>3)-1,[conTOdiv{n_perint>3}],cellfun(@(x) x(1),conTOdiv_CI(n_perint>3))-[conTOdiv{n_perint>3}],...
    cellfun(@(x) x(2),conTOdiv_CI(n_perint>3))-[conTOdiv{n_perint>3}],'color','k'); hold on
scatter(find(n_perint>3)-1,[conTOdiv{n_perint>3}],'k','filled');

% 7 day stripe rearing
mouse_group = [1 2 3 4 5 6 7 8];
maxtp = 6; %max(cellfun(@(x) size(x(1).date,1),data_orientation_deprivation(mouse_group))); %
mintp = 6; %min(cellfun(@(x) size(x(1).date,1),data_orientation_deprivation(mouse_group))); %

time_point_labels = {'BL1', 'BL2', 'SR1', 'SR2', 'SR3', 'SR4'};
[~, ~, ~, ~,conTOdiv,conTOdiv_CI] = StripeRearing_stability(data_orientation_deprivation(mouse_group),time_point_labels,...
    'comparison_type','BL2-All','Binning_interval' , 30,...
    'cell_filter','resp_anytime','use_only_sig_changes',0);

figure(203); set(gcf,'color','w');
errorbar(7:7:7*(maxtp-2),[conTOdiv{2:end}],cellfun(@(x) x(1),conTOdiv_CI(2:end))-[conTOdiv{2:end}],...
    cellfun(@(x) x(2),conTOdiv_CI(2:end))-[conTOdiv{2:end}],'r')
scatter(7:7:7*(maxtp-2),[conTOdiv{2:end}],'r','filled'); hold on

% 28d Stripe rearing
mouse_group = [9 10 11 14 16 17 18];
maxtp = 3;
mintp = 3;

time_point_labels = {'BL1', 'BL2', 'SR1'};
[~, ~, ~, ~,conTOdiv,conTOdiv_CI] = StripeRearing_stability(data_orientation_deprivation(mouse_group),time_point_labels,...
    'comparison_type','BL2-All','Binning_interval' , 30,...
    'cell_filter','resp_anytime','use_only_sig_changes',0);
figure(203);
errorbar(28,[conTOdiv{2:end}],cellfun(@(x) x(1),conTOdiv_CI(2:end))-[conTOdiv{2:end}],...
    cellfun(@(x) x(2),conTOdiv_CI(2:end))-[conTOdiv{2:end}],'.g'); hold on
plot(28,[conTOdiv{2:end}],'--g'); hold on
scatter(28,[conTOdiv{2:end}],'g','filled'); hold on

xlim([0,30])
ylim([-2 5])
xlabel('SR length or dT (days)')
ylabel('median d|rPO| [°]')
set(gca,'TickDir','out')
plot(xlim,[0 0],'--k')
box off

%% Magnitude and direction shuffle & optimization
fig_number = fig_number+1;
fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'magnitude vs direction shuffle', 'Position', [790, 682, 560, 292])

for SRdays = 1:2
    if SRdays==2
        mouse_group = [9 10 11 14 16 17 18];
        tp_color = [0.4660 0.6740 0.1880];
    elseif SRdays==1
        mouse_group = [1 2 3 4 5 6 7 8]; % switch to 7 day interval
        tp_color = [0.8500 0.3250 0.0980];
    end
    maxtp = 3;
    mintp = 3;
    time_point_labels = {'BL1', 'BL2', 'SR1'};
    
    % single cell recovery
    [~, ~, ~, comp_type_label,conTOdiv,conTOdiv_CI,~,~,data_out_all] = StripeRearing_stability(data_orientation_deprivation(mouse_group),time_point_labels,...
        'comparison_type','Sequential','Binning_interval' , 30,...
        'cell_filter','resp_anytime','use_only_sig_changes',0);
    
    dPOperm = data_out_all{2}(1,:); % rPO start
    dPO_orig = data_out_all{2}(3,:); % d(rPO)
    idx_nan = isnan(abs(dPO_orig)); % remove nans
    dPOperm(idx_nan)=[]; % remove nans
    dPO_orig(idx_nan)=[]; % remove nans
    issig = data_out_all{2}(4,:); % sig dPOs
    issig(idx_nan)=[]; % remove nans
    
    avg_func=@nanmedian; % avg metric to use (can switch to mean)
    
    % can switch to only using significant chagnes
    onlysig = 0;
    if onlysig
        idx_0 = ~logical(issig);
        dPOperm(idx_0)=[];
        dPO_orig(idx_0)=[];
        issig(idx_0)=[];
    end
    
    % Dir: original vs optimized vs shuffled
    dPO = dPO_orig;
    
    % original
    dPO_abs = abs(dPO);
    dPOperm2 = wrapTo180((dPOperm+dPO)*2)/2; % equivalent too end rPO: data_out_all{2}(2,:);
    dabsrPO = wrapTo180((abs(dPOperm)-abs(dPOperm2))*2)/2;  % d(abs(rPO)) normal data
    
    % randomize
    dPO_abs_rand = dPO_abs.*sign(dPO(randperm(length(dPO))));
    dPOperm2_rand = wrapTo180((dPOperm+dPO_abs_rand)*2)/2;
    dabsrPO_rand = wrapTo180((abs(dPOperm)-abs(dPOperm2_rand))*2)/2;
    
    % optimize
    dPO_abs_opt = dPO_abs.*(-1*sign(dPOperm));
    dPOperm2_opt = wrapTo180((dPOperm+dPO_abs_opt)*2)/2;
    dabsrPO_opt = wrapTo180((abs(dPOperm)-abs(dPOperm2_opt))*2)/2;
    
    subplot(1,2,2)
    dabsrPO_ci = bootci(1000,avg_func,dabsrPO)-avg_func(dabsrPO);
    dabsrPO_rand_ci = bootci(1000,avg_func,dabsrPO_rand)-avg_func(dabsrPO_rand);
    dabsrPO_opt_ci = bootci(1000,avg_func,dabsrPO_opt)-avg_func(dabsrPO_opt);
    errorbar([1,2,3],[avg_func(dabsrPO),avg_func(dabsrPO_rand),avg_func(dabsrPO_opt)],...
        [dabsrPO_ci(1),dabsrPO_rand_ci(1),dabsrPO_opt_ci(1)],...
        [dabsrPO_ci(2),dabsrPO_rand_ci(2),dabsrPO_opt_ci(2)],'Color',tp_color,'LineWidth',2); hold on
    scatter([1,2,3],[median(dabsrPO),median(dabsrPO_rand),median(dabsrPO_opt)],'MarkerFaceColor',tp_color,'MarkerEdgeColor',tp_color)
    box off
    xlim([0.5 3.5]);plot(xlim,[0 0],'--k')
    ylim([-2 8])
    xticks([1 2 3])
    xticklabels({'data', 'rand dPO dir','opt dPO dir'}); xtickangle(45)
    ylabel(['median d|rPO| [°]'])
    [p,~,STATS]=signrank(dabsrPO,dabsrPO_rand);
    title(['z' num2str(STATS.zval,3)  ', p' num2str(p,'%.3f')])

    % Magnitude: original vs optimized vs shuffled
    
    % original
    dPO = dPO_orig;
    dPO_abs = abs(dPO);
    dPOperm2 = wrapTo180((dPOperm+dPO)*2)/2; % equivalent too end rPO: data_out_all{2}(2,:);
    dabsrPO = wrapTo180((abs(dPOperm)-abs(dPOperm2))*2)/2;  % d(abs(rPO)) normal data
    
    % randomize
    dPO = dPO_orig;
    dPO = abs(dPO(randperm(length(dPO)))).*sign(dPO); % shuffle drift magnitude but not direction
    dPO_abs_rand = abs(dPO);
    dPOperm2_rand = wrapTo180((dPOperm+dPO)*2)/2; % equivalent too end rPO: data_out_all{2}(2,:);
    dabsrPO_rand = wrapTo180((abs(dPOperm)-abs(dPOperm2_rand))*2)/2;  % d(abs(rPO)) convergence randomized mag
    
    % optimize
    dPO = dPO_orig;
    [dPO_abs,I] = sort(abs(dPO)); % sort based on dPO size
    [temp,I]=sort(abs(dPOperm)); % sorting of start |rPO|
    dPO_abs(I) = dPO_abs; % sort based on |rPO|
    dPO = dPO_abs.*sign(dPO);
    dPO_abs_opt = abs(dPO);
    dPOperm2_opt = wrapTo180((dPOperm+dPO)*2)/2;
    dabsrPO_opt = wrapTo180((abs(dPOperm)-abs(dPOperm2_opt))*2)/2;
    
    subplot(1,2,1); 
    dabsrPO_ci = bootci(1000,avg_func,dabsrPO)-avg_func(dabsrPO);
    dabsrPO_rand_ci = bootci(1000,avg_func,dabsrPO_rand)-avg_func(dabsrPO_rand);
    dabsrPO_opt_ci = bootci(1000,avg_func,dabsrPO_opt)-avg_func(dabsrPO_opt);
    errorbar([1,2,3],[avg_func(dabsrPO),avg_func(dabsrPO_rand),avg_func(dabsrPO_opt)],...
        [dabsrPO_ci(1),dabsrPO_rand_ci(1),dabsrPO_opt_ci(1)],...
        [dabsrPO_ci(2),dabsrPO_rand_ci(2),dabsrPO_opt_ci(2)],'Color',tp_color,'LineWidth',2); hold on
    scatter([1,2,3],[median(dabsrPO),median(dabsrPO_rand),median(dabsrPO_opt)],'MarkerFaceColor',tp_color,'MarkerEdgeColor',tp_color)
    box off
    xlim([0.5 3.5]);plot(xlim,[0 0],'--k')
    ylim([-2 8])
    xticks([1 2 3])
    xticklabels({'data', 'rand dPO mag','opt dPO mag'}); xtickangle(45)
    ylabel(['median d|rPO| [°]'])
    [p,~,STATS]=signrank(dabsrPO,dabsrPO_rand);
    title(['z' num2str(STATS.zval,3)  ', p' num2str(p,'%.3f')])
end

%% recovery from Stripe rearing
mouse_group = [9 10 11 14 16 17 18]; % 28d interval
maxtp = 4;
mintp = 4;

idx_BL2 = 2;
SR_idx = 3;
RE_idx = 4;

time_point_labels = {'BL1', 'BL2', 'SR1', 'RE1'};

% population effect
fig_number = fig_number+1;
fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'recovery population effect', 'Position', [155,149,1645,769]);

clear PO_bin_frac_sig_ofTuned
ax = [];
for mouse_n = mouse_group
    if ~isempty(data_orientation_deprivation{mouse_n})
        [PO_to_PrefO_count_per_mouse_bin90{mouse_n},~] = StripeRearing_stability(data_orientation_deprivation(mouse_n),time_point_labels(1:mintp),...
            'comparison_type','BL2-All','Binning_interval' , 90,...
            'cell_filter','resp_anytime','use_only_sig_changes',0);
        [PO_to_PrefO_count_per_mouse_bin30{mouse_n},~] = StripeRearing_stability(data_orientation_deprivation(mouse_n),time_point_labels(1:mintp),...
            'comparison_type','BL2-All','Binning_interval' , 30,...
            'cell_filter','resp_anytime','use_only_sig_changes',0);
        
        if i == 1 % get
            [~,~,PO_bin_frac_sig_ofTuned{mouse_n}] = StripeRearing_stability(data_orientation_deprivation(mouse_n),time_point_labels(1:mintp),...
                'comparison_type','BL2-All','Binning_interval' , 90,...
                'cell_filter','resp_anytime','use_only_sig_changes',0);
        end
    end
end
[PO_to_PrefO_count,PO_bin_centers] = StripeRearing_stability(data_orientation_deprivation(mouse_group),time_point_labels(1:mintp),...
    'comparison_type','BL2-All','Binning_interval', 90,...
    'cell_filter','resp_anytime','use_only_sig_changes',0);

changeBLSR_PermOri_per_mouse = cellfun(@(x) (x(SR_idx,2) - x(idx_BL2,2)),PO_to_PrefO_count_per_mouse_bin90(mouse_group));
changeBLSR_OrthOri_per_mouse = cellfun(@(x) (x(SR_idx,1) - x(idx_BL2,1)),PO_to_PrefO_count_per_mouse_bin90(mouse_group));
changeSRRE_PermOri_per_mouse = cellfun(@(x) (x(RE_idx,2) - x(SR_idx,2)),PO_to_PrefO_count_per_mouse_bin90(mouse_group));
changeSRRE_OrthOri_per_mouse = cellfun(@(x) (x(RE_idx,1) - x(SR_idx,1)),PO_to_PrefO_count_per_mouse_bin90(mouse_group));
changeBLRE_PermOri_per_mouse = cellfun(@(x) (x(RE_idx,2) - x(idx_BL2,2)),PO_to_PrefO_count_per_mouse_bin90(mouse_group));
changeBLRE_OrthOri_per_mouse = cellfun(@(x) (x(RE_idx,1) - x(idx_BL2,1)),PO_to_PrefO_count_per_mouse_bin90(mouse_group));

cellcount_BL2=num2cell(cellfun(@(x)  x(idx_BL2,1) + x(idx_BL2,2), PO_to_PrefO_count_per_mouse_bin90(mouse_group)));
cellcount_SR=num2cell(cellfun(@(x)  x(SR_idx,1) + x(SR_idx,2), PO_to_PrefO_count_per_mouse_bin90(mouse_group)));
percent_changeBLSR_PermOri_per_mouse = cellfun(@(x,y) (x(SR_idx,2) - x(idx_BL2,2))/(y),PO_to_PrefO_count_per_mouse_bin90(mouse_group),cellcount_BL2);
percent_changeBLSR_OrthOri_per_mouse = cellfun(@(x,y) (x(SR_idx,1) - x(idx_BL2,1))/(y),PO_to_PrefO_count_per_mouse_bin90(mouse_group),cellcount_BL2);
percent_changeSRRE_PermOri_per_mouse = cellfun(@(x,y) (x(RE_idx,2) - x(SR_idx,2))/(y),PO_to_PrefO_count_per_mouse_bin90(mouse_group),cellcount_SR);
percent_changeSRRE_OrthOri_per_mouse = cellfun(@(x,y) (x(RE_idx,1) - x(SR_idx,1))/(y),PO_to_PrefO_count_per_mouse_bin90(mouse_group),cellcount_SR);
percent_changeBLRE_PermOri_per_mouse = cellfun(@(x,y) (x(RE_idx,2) - x(idx_BL2,2))/(y),PO_to_PrefO_count_per_mouse_bin90(mouse_group),cellcount_BL2);
percent_changeBLRE_OrthOri_per_mouse = cellfun(@(x,y) (x(RE_idx,1) - x(idx_BL2,1))/(y),PO_to_PrefO_count_per_mouse_bin90(mouse_group),cellcount_BL2);

subplot(2,3,1:3); cla
if maxtp == 3
    colors = [0.5 0.5 0.5;0 0 0; 1 0 0];
elseif maxtp == 4
    colors = [0.5 0.5 0.5;0 0 0; 1 0 0;1 0.6 0];
elseif maxtp > 4
    colors = [0.5 0.5 0.5; 0 0 0; flip(autumn(maxtp-2))];
end

for tp = [idx_BL2 SR_idx RE_idx]
    all_dist = cellfun(@(x) x(tp,:)/sum(x(2,:)), PO_to_PrefO_count_per_mouse_bin30(mouse_group),'UniformOutput', false);
    
    errorbar([-90:30:60],mean(cat(1,all_dist{:})).*100,std(cat(1,all_dist{:}).*100)./sqrt(length(mouse_group)),'color',colors(tp,:),'LineWidth',2); hold on
    plot([-90:30:60],mean(cat(1,all_dist{:}).*100,1),'color',colors(tp,:),'LineWidth',2)
    
    set(gca,'TickDir','out'); box off; set(gca,'TickDir','out')
    xticks([-90:30:60]); xlim([-100 70]);ylim([0 inf])
    xlabel([{'pref. Ori. (PO) diff.'};{'to Permitted Ori'}])
    if i==1
        ylabel({['% relative to'] ;['total count at BL2']})
    end
end
pbaspect([1.5,0.75,0.5])
ylim([0 35])
title('resp_anytime','Interpreter','none')

ax_temp = subplot(2,3,4); cla; ax = cat(1,ax,ax_temp);
for mouse_n = 1:length(mouse_group)
    plot([2,1],100.*[percent_changeBLSR_OrthOri_per_mouse(mouse_n), percent_changeBLSR_PermOri_per_mouse(mouse_n)],'color',[0.5 0.5 0.5]); hold on
end
plot([2,1],100.*[mean(percent_changeBLSR_OrthOri_per_mouse), mean(percent_changeBLSR_PermOri_per_mouse)],'k','LineWidth',5); hold on
xticks([1,2]);xticklabels({'0' '90'})
xlim([0.5 2.5]); ylim([-40 40])
plot(xlim,[0 0],'--k')
box off
set(gca,'TickDir','out')
xlabel([{'pref. Ori. (PO) diff.'};{'to Permitted Ori'}])
ylabel({'% change in cell count'; 'from BL to SR'})
ylim([-1*max(abs(ylim)) max(abs(ylim))])
pbaspect([0.75,1.3,0.5])
%stats
[~,p_90,~,STATS_90]=ttest([percent_changeBLSR_OrthOri_per_mouse]);
[~,p_0,~,STATS_0]=ttest([percent_changeBLSR_PermOri_per_mouse]);
[~,p_paired,~,STATS_paired]=ttest([percent_changeBLSR_PermOri_per_mouse-percent_changeBLSR_OrthOri_per_mouse]);
title({['t(' num2str(STATS_paired.df) ')' num2str(STATS_paired.tstat,3) ',p=' num2str(p_paired,'%.3f')];...
    ['t(' num2str(STATS_0.df) ')' num2str(STATS_0.tstat,3) ',p=' num2str(p_0,'%.3f') ...
    '; t(' num2str(STATS_90.df) ')' num2str(STATS_90.tstat,3) ',p=' num2str(p_90,'%.3f')]...
    },'Interpreter','none','FontSize',10)

ax_temp = subplot(2,3,5); cla; ax = cat(1,ax,ax_temp);
for mouse_n = 1:length(mouse_group)
    plot([2,1],100.*[percent_changeSRRE_OrthOri_per_mouse(mouse_n), percent_changeSRRE_PermOri_per_mouse(mouse_n)],'color',[0.5 0.5 0.5]); hold on
end
plot([2,1],100.*[mean(percent_changeSRRE_OrthOri_per_mouse), mean(percent_changeSRRE_PermOri_per_mouse)],'k','LineWidth',5); hold on
xticks([1,2]);xticklabels({'0' '90'})
xlim([0.5 2.5]); ylim([-40 40])
plot(xlim,[0 0],'--k')
box off
set(gca,'TickDir','out')
xlabel([{'pref. Ori. (PO) diff.'};{'to Permitted Ori'}])
ylabel({'% change in cell count'; 'from SR to RE'})
ylim([-1*max(abs(ylim)) max(abs(ylim))])
pbaspect([0.75,1.3,0.5])
%stats
[~,p_90,~,STATS_90]=ttest([percent_changeSRRE_OrthOri_per_mouse]);
[~,p_0,~,STATS_0]=ttest([percent_changeSRRE_PermOri_per_mouse]);
[~,p_paired,~,STATS_paired]=ttest([percent_changeSRRE_PermOri_per_mouse-percent_changeSRRE_OrthOri_per_mouse]);
title({['t(' num2str(STATS_paired.df) ')' num2str(STATS_paired.tstat,3) ',p=' num2str(p_paired,'%.3f')];...
    ['t(' num2str(STATS_0.df) ')' num2str(STATS_0.tstat,3) ',p=' num2str(p_0,'%.3f') ...
    '; t(' num2str(STATS_90.df) ')' num2str(STATS_90.tstat,3) ',p=' num2str(p_90,'%.3f')]...
    },'Interpreter','none','FontSize',10)

ax_temp = subplot(2,3,6); cla; ax = cat(1,ax,ax_temp);
for mouse_n = 1:length(mouse_group)
    plot([2,1],100.*[percent_changeBLRE_OrthOri_per_mouse(mouse_n), percent_changeBLRE_PermOri_per_mouse(mouse_n)],'color',[0.5 0.5 0.5]); hold on
end
plot([2,1],100.*[mean(percent_changeBLRE_OrthOri_per_mouse), mean(percent_changeBLRE_PermOri_per_mouse)],'k','LineWidth',5); hold on
xticks([1,2]);xticklabels({'0' '90'})
xlim([0.5 2.5]); ylim([-40 40])
plot(xlim,[0 0],'--k')
box off
set(gca,'TickDir','out')
xlabel([{'pref. Ori. (PO) diff.'};{'to Permitted Ori'}])
ylabel({'% change in cell count'; 'from BL to RE'})
ylim([-1*max(abs(ylim)) max(abs(ylim))])
pbaspect([0.75,1.3,0.5])
%stats
[~,p_90,~,STATS_90]=ttest([percent_changeBLRE_OrthOri_per_mouse]);
[~,p_0,~,STATS_0]=ttest([percent_changeBLRE_PermOri_per_mouse]);
[~,p_paired,~,STATS_paired]=ttest([percent_changeBLRE_PermOri_per_mouse-percent_changeBLRE_OrthOri_per_mouse]);
title({['t(' num2str(STATS_paired.df) ')' num2str(STATS_paired.tstat,3) ',p=' num2str(p_paired,'%.3f')];...
    ['t(' num2str(STATS_0.df) ')' num2str(STATS_0.tstat,3) ',p=' num2str(p_0,'%.3f') ...
    '; t(' num2str(STATS_90.df) ')' num2str(STATS_90.tstat,3) ',p=' num2str(p_90,'%.3f')]...
    },'Interpreter','none','FontSize',10)

linkaxes(ax)