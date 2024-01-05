function[vect_PO, vect_POCI_dist, vect_POCI, CI_under_45, OS, OS_CI, OSnn, OSnn_CI, vect_PD, pref_stim, pref_stim_resp, DS] = ...
    vect_PO_boot(mean_resp_dir_rep,stim_dir_correction,showplot)
% in order to keep rep number equal across stims, if there is a nan in a
% rep for any stim, that whole rep is removed 
mean_resp_dir_rep(:,find(any(isnan(mean_resp_dir_rep),1))) = nan;

% first stimulus isdir 0 which is right to left, angles increasing clockwise
% dirs = deg2rad([360/size(mean_resp_dir_rep,1):360/size(mean_resp_dir_rep,1):360]);
% dirs = circshift(dirs,-5); 
dirs = deg2rad(0:360/size(mean_resp_dir_rep,1):360-360/size(mean_resp_dir_rep,1));

if length(stim_dir_correction)==1
    dirs = wrapToPi(dirs-deg2rad(round(stim_dir_correction))); % convert to relative direction (instead of screen)
    dirs_all = repmat(dirs,size(mean_resp_dir_rep,2),1)';
elseif length(stim_dir_correction)==2 % different angles for the two aquisistions
    dirs_all = repmat(dirs,size(mean_resp_dir_rep,2),1)';
    dirs_all(:,1:16) = wrapToPi(dirs_all(:,1:16)-deg2rad(round(stim_dir_correction(1)))); % not sure if this is compensating in the right direction, check!
    if ~isnan(stim_dir_correction(2))
        dirs_all(:,17:32) = wrapToPi(dirs_all(:,17:32)-deg2rad(round(stim_dir_correction(2)))); % not sure if this is compensating in the right direction, check!
    end
    dirs = wrapToPi(dirs-deg2rad(round(nanmean(stim_dir_correction)))); % convert to relative direction (instead of screen)
end
dirs_all = dirs_all(:);

oris = round(wrapToPi(dirs.*2)./2,14); % if i dont round then the two directions with equivalent ori are slightly different
oris_all = repmat(oris,size(mean_resp_dir_rep,2),1)';
oris_all = oris_all(:);
% should make oris & dirs outputs to be saved in main data structure

% two alternatives! baseline subract or set neg responses to 0
% 1)
mean_resp_dir_rep_all = mean_resp_dir_rep(:)-min(mean_resp_dir_rep(:)); % min subtracted 
% 2) not a fan of this one
% mean_resp_dir_rep_all_nonneg = mean_resp_dir_rep(:);
% mean_resp_dir_rep_all_nonneg(mean_resp_dir_rep_all_nonneg(:)<0) = 0; % neg values and nans set to 0

idx_nan = find(isnan(mean_resp_dir_rep_all));
% if ~isempty(idx_nan) & ~all(diff(hist(oris_all,unique(oris_all)))==0)
%    error('nans in responseve vector detected, check for unequal stimulus representation') 
% end
mean_resp_dir_rep_all(idx_nan) = []; % this is dangerouse if there is an eneven dist of nans across stims
oris_all(idx_nan) = [];
dirs_all(idx_nan) = [];

[X,Y]=pol2cart(oris_all*2,mean_resp_dir_rep_all);
[ori_mean]=cart2pol(mean(X),mean(Y));
ori_mean = wrapToPi(ori_mean./2);

%% PO
[X,Y]=pol2cart((oris_all-ori_mean)*2,mean_resp_dir_rep_all);
% [ori_mean0]=cart2pol(mean(X),mean(Y));

boots_on = @(x,y) wrapToPi(cart2pol(mean(x),mean(y))./2); % gives 
[STATS_PO]= bootstrp(1000,boots_on,X(:),Y(:)); % output is in actual angle so needs to be *2 for circ stats
ori_bootmean = circ_mean(STATS_PO.*2)./2;
nouniform_boot = circ_rtest(STATS_PO.*2)<0.05;
% ori_CI = prctile(abs(STATS_PO),95).*[-1; 1]+ori_bmean;
ori_CI = [prctile(STATS_PO,5); prctile(STATS_PO,95)];

vect_PO = wrapToPi([ori_bootmean+ori_mean]*2)/2;
vect_POCI = wrapToPi([ori_CI+ori_mean]*2)/2;
vect_POCI_dist = circ_dist(vect_POCI'.*2,[vect_PO vect_PO].*2)./2;
% CI_under_45 = rad2deg(abs(vect_POCI_dist(1)))<45;
CI_under_45 = rad2deg(mean(abs(vect_POCI_dist)))<45;

%% OS (orientation selectivity)

% multiple methods

% OSI formula
OS_fun = @(theta,R) abs( sum(mean(R,1).*exp(2*1i*mean(theta,1))) / sum(mean(R,1)) ) ; % from Mazurek et al., 2014
% OSI formula with PO subtraction to make sure PO->OS correlation is not an artefact
% OS_fun = @(theta,R) abs( sum(mean(R,1).*exp(2*1i* (mean(theta,1)-vect_PO))  ) / sum(mean(R,1)) ) ; % from Mazurek et al., 2014

% min response subracted
OS = OS_fun(rad2deg(dirs),mean_resp_dir_rep'-min(mean_resp_dir_rep(:)));
[STATS_OS] = bootstrp(1000,OS_fun,...
    repmat(rad2deg(dirs),size(mean_resp_dir_rep,2),1),...
    mean_resp_dir_rep'-min(mean_resp_dir_rep(:)));
OS_bootmean = mean(STATS_OS);
OS_CI = [prctile(STATS_OS,5); prctile(STATS_OS,95)];

% 2 OSnonneg (setting all neg resp to 0)
mean_resp_dir_rep_nonneg = mean_resp_dir_rep;
mean_resp_dir_rep_nonneg(mean_resp_dir_rep_nonneg<0) = 0;
OSnn = OS_fun(rad2deg(dirs),mean_resp_dir_rep_nonneg');
[STATS_OSnn] = bootstrp(1000,OS_fun,...
    repmat(rad2deg(dirs),size(mean_resp_dir_rep_nonneg,2),1),...
    mean_resp_dir_rep_nonneg'-min(mean_resp_dir_rep_nonneg(:)));
OSnn_bootmean = mean(STATS_OSnn);
OSnn_CI = [prctile(STATS_OSnn,5); prctile(STATS_OSnn,95)];

%% PD & DS - did not boostrap this to get CI for now
% there are several posibilities of how to calculate PD and DS. here i have
% gone for the most similar to previous studies iv found

% % based on closest stimulus to the vect based PO
pref_Ori_stims_n = find(abs(oris-vect_PO) == min(abs(oris-vect_PO)));
if length(pref_Ori_stims_n)<2
    pref_Ori_stims_n = [pref_Ori_stims_n, rem((pref_Ori_stims_n + length(dirs)/2),length(dirs))];
    if any(pref_Ori_stims_n==0)
        pref_Ori_stims_n(pref_Ori_stims_n==0)=length(dirs);
    end
end
pref_Ori_stims = dirs(pref_Ori_stims_n);

% based on the stimulus with the largest response amplitude
% [~,pref_Ori_stims_n] = max(nanmean(mean_resp_dir_rep-min(mean_resp_dir_rep(:)),2));
% pref_Ori_stims_n = mod([pref_Ori_stims_n,pref_Ori_stims_n+size(mean_resp_dir_rep,1)/2],size(mean_resp_dir_rep,1));
% pref_Ori_stims_n(pref_Ori_stims_n==0)=size(mean_resp_dir_rep,1);
% pref_Ori_stims = dirs(pref_Ori_stims_n);


idx_a = dirs_all==pref_Ori_stims(1);
idx_b = dirs_all==pref_Ori_stims(2);
resp_temp = [mean(mean_resp_dir_rep_all(idx_a)),mean(mean_resp_dir_rep_all(idx_b))]; % mean resp to prefered dir and null direction

[~,idx] = max(resp_temp);
pref_stim = pref_Ori_stims(idx);
pref_stim_resp = resp_temp(idx);

if abs(wrapToPi(pref_stim))>pi/2
    vect_PD = wrapToPi(vect_PO+pi);
else
    vect_PD = vect_PO;
end

DS = (abs(diff(resp_temp)))/(sum(resp_temp)); % difference over sum of prefered direction vs null
if DS>1
    DS=1;
elseif DS<0;
    DS=0;
end

%%
if showplot==1
    figure(101); hold off; 
    histogram(STATS_OS); hold on; ylabel('count'); xlabel('OS')
    plot([1 1].*OS_bootmean,ylim,'b'); 
    plot([1 1].*OS,ylim,'r'); 
    plot([1 1].*OS_CI(1),ylim,'--b'); 
    plot([1 1].*OS_CI(2),ylim,'--b')
    
    histogram(STATS_OSnn); hold on; ylabel('count'); xlabel('OS')
    plot([1 1].*OSnn_bootmean,ylim,'b'); 
    plot([1 1].*OSnn,ylim,'r'); 
    plot([1 1].*OSnn_CI(1),ylim,'--b'); 
    plot([1 1].*OSnn_CI(2),ylim,'--b')
    xlim([0 1])
    set(gcf,'color','w');
    set(gca,'FontSize',15);
    
    figure(102); hold off
    a = polarscatter(oris_all.*2,(mean_resp_dir_rep_all(:))./(max((mean_resp_dir_rep_all(:)))).*2); hold on
    a.MarkerEdgeColor = [0.65 0.65 0.65];
    a.MarkerFaceColor = [0.65 0.65 0.65];
    a.MarkerEdgeAlpha = 0.3;
    a.MarkerFaceAlpha = 0.3;
%     polarplot([ori_mean ori_mean].*2, [rlim],'b')
    polarplot([vect_PO vect_PO].*2, [rlim],'k')
    polarplot([ori_CI(1)+ori_mean ori_CI(1)+ori_mean].*2, [rlim],'--k')
    polarplot([ori_CI(2)+ori_mean ori_CI(2)+ori_mean].*2, [rlim],'--k')
    temp = nanmean(reshape(nanmean(mean_resp_dir_rep-min(mean_resp_dir_rep(:)),2),6,2),2)./max((mean_resp_dir_rep_all(:))).*2;
    polarplot(oris([1:6,1]).*2,temp([1:6,1]),'-k')
%     polarhistogram(STATS.*2+ori_mean.*2,'BinWidth',2*pi/12,'FaceColor','b','Normalization','probability')
    polarhistogram(STATS_PO.*2+ori_mean.*2,'BinWidth',2*pi/24,'FaceColor','k','Normalization','probability')
    thetaticks(sort(rad2deg(oris(1:end/2).*2)))
    thetaticklabels(num2cell(wrapTo180(sort(rad2deg(oris(1:end/2))))))
    set(gcf,'color','w');
    set(gca,'FontSize',15);
    
    figure(103); hold off
    a = polarscatter(dirs_all,(mean_resp_dir_rep_all(:))./(max((mean_resp_dir_rep_all(:))))); hold on
    a.MarkerEdgeColor = [0.65 0.65 0.65];
    a.MarkerFaceColor = [0.65 0.65 0.65];
    a.MarkerEdgeAlpha = 0.3;
    a.MarkerFaceAlpha = 0.3;
    polarplot([pref_stim pref_stim], [rlim],'k')
    polarplot([vect_PD vect_PD], [rlim],'b')
    temp = nanmean(mean_resp_dir_rep-min(mean_resp_dir_rep(:)),2)./max((mean_resp_dir_rep_all(:)));
    polarplot(dirs([1:12,1]),temp([1:12,1]),'-k')
    thetaticks(sort(rad2deg(dirs)))
    thetaticklabels(num2cell(sort(rad2deg(dirs))))
    set(gcf,'color','w');
    set(gca,'FontSize',15);
    
    figure(104); cla
    plot(rad2deg(dirs),mean(mean_resp_dir_rep,2)); hold on
    xlabel('stim. dir.'); ylabel('resp amp')
    set(gca,'TickDir','out')
    set(gcf,'color','w');
    set(gca,'FontSize',15);    
    plot(repmat(rad2deg(vect_PO),1,2),ylim,'k')
    plot(repmat(rad2deg(vect_PO+pi),1,2),ylim,'--k')
end
end
