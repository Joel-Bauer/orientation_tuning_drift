function [frac_sig_bothtuned,...
    dT_bin_labels,...
    ratio_conTOdiv,...
    ratio_conTOdiv_CI,...
    dPOsig_perint_binned,...
    dPOall_perint_binned,...
    output] = ...
    eval_PO_stability(data,varargin)

p = inputParser();
p.KeepUnmatched = false;
p.CaseSensitive = false;
p.StructExpand  = false;
p.PartialMatching = true;

addParameter(p, 'dT_binsize' , 1)
addParameter(p, 'cell_filter' , 'resp_anytime')
addParameter(p, 'startPO_sig_frac' , 0)
addParameter(p, 'max_interval' , max(cellfun(@(x) max(x(1).day), data)))
addParameter(p, 'PO_relative' , [])
addParameter(p, 'use_only_sig_changes', 1)

parse(p, varargin{:})
names = fieldnames(p.Results);
for i=1:length(names)
    eval([names{i} '=p.Results.' names{i} ]);
end

% get running
for mouse_n = 1:length(data)
    try
        mouse_running{mouse_n}.mouse = data{mouse_n}(1).mouse;
        mouse_running{mouse_n}.all_Running_traces_stim = data{mouse_n}(1).all_Running_traces_stim;
        mouse_running{mouse_n}.all_Running_traces_poststim = data{mouse_n}(1).all_Running_traces_poststim;
    catch
        mouse_running{mouse_n}.all_Running_traces_stim = [];
        mouse_running{mouse_n}.all_Running_traces_poststim = [];
    end
end

% get pupilsize
for mouse_n = 1:length(data)
    try
        mouse_pupilsize{mouse_n}.mouse = data{mouse_n}(1).mouse;
        mouse_pupilsize{mouse_n}.all_Pupilsize_traces_stim = data{mouse_n}(1).all_Pupilsize_traces_stim;
        mouse_pupilsize{mouse_n}.all_Pupilsize_traces_poststim = data{mouse_n}(1).all_Pupilsize_traces_poststim;
    catch
        mouse_pupilsize{mouse_n}.all_Pupilsize_traces_stim = [];
        mouse_pupilsize{mouse_n}.all_Pupilsize_traces_poststim = [];
    end
end

% subselect cells based on filter
if strcmp(cell_filter,'continously_resp')
    for mouse_n = 1:length(data)
        idx_filter = cellfun(@(x) all(x),{data{mouse_n}.visually_responsive});
        data{mouse_n} = data{mouse_n}(find(idx_filter));
    end
elseif strcmp(cell_filter,'continously_tuned')
    for mouse_n = 1:length(data)
        idx_filter = cellfun(@(x) all(x),{data{mouse_n}.Stat_tuned});
        data{mouse_n} = data{mouse_n}(find(idx_filter));
    end
elseif strcmp(cell_filter,'continously_resp_and_tuned')
    for mouse_n = 1:length(data)
        idx_filter = cellfun(@(x) all(x),{data{mouse_n}.visually_responsive});
        data{mouse_n} = data{mouse_n}(find(idx_filter));
    end
    for mouse_n = 1:length(data)
        idx_filter = cellfun(@(x) all(x),{data{mouse_n}.Stat_tuned});
        data{mouse_n} = data{mouse_n}(find(idx_filter));
    end
elseif strcmp(cell_filter,'resp_anytime')
    for mouse_n = 1:length(data)
        idx_filter = cellfun(@(x) any(x),{data{mouse_n}.visually_responsive});
        data{mouse_n} = data{mouse_n}(find(idx_filter));
    end
elseif strcmp(cell_filter,'tuned_anytime')
    for mouse_n = 1:length(data)
        idx_filter = cellfun(@(x) any(x),{data{mouse_n}.Stat_tuned});
        data{mouse_n} = data{mouse_n}(find(idx_filter));
    end
end

% calc running modulation index per cell (run_mod has a wider range so use that for now)
clear  run_mod run_mod_tpmeans pupil_mod pupil_mod_tpmeans 
for mouse_n = 1:length(data)
    if ~isempty(mouse_running{mouse_n}.all_Running_traces_stim)
        running_all = cellfun(@(x,y) cat(1,x,y),mouse_running{mouse_n}.all_Running_traces_stim, mouse_running{mouse_n}.all_Running_traces_poststim,'UniformOutput',0)';
        pupilsize_all = cellfun(@(x,y) cat(1,x,y),mouse_pupilsize{mouse_n}.all_Pupilsize_traces_stim, mouse_pupilsize{mouse_n}.all_Pupilsize_traces_poststim,'UniformOutput',0)';
        
        mean_speed{mouse_n} = cellfun(@(x) mean(x(:)),running_all);
        mean_pupilsize{mouse_n} = cellfun(@(x) mean(x(:)),pupilsize_all);

        mean_dSpeed{mouse_n} = mean_speed{mouse_n}'-mean_speed{mouse_n}; % mean_speed_dif{mouse_n}(1,2) change in mean running speed from tp1 to tp2
        mean_dPupil{mouse_n} = mean_pupilsize{mouse_n}'-mean_pupilsize{mouse_n}; 
        mean_dBehav{mouse_n} = nanmean(cat(3,mean_dSpeed{mouse_n}, mean_dPupil{mouse_n}),3);
        
        for cell_n = 1:length(data{mouse_n})
            mouse_resp_all = cellfun(@(x,y) cat(1,x,y),data{mouse_n}(cell_n).all_traces_stim,data{mouse_n}(cell_n).all_traces_poststim,'UniformOutput',0);
            for tp = 1:length(mouse_resp_all)
                temp_run = mean(running_all{tp},1);
                temp_pupil = mean(pupilsize_all{tp},1);
                temp_activity = mean(mouse_resp_all{tp},1);
                run_mod{mouse_n}(cell_n,tp) = corr(temp_run(:),temp_activity(:));
                pupil_mod{mouse_n}(cell_n,tp) = corr(temp_pupil(:),temp_activity(:));
            end
            run_mod_tpmeans{mouse_n}(cell_n,:,:) = repmat(max(run_mod{mouse_n}(cell_n,:)),length(mouse_resp_all),length(mouse_resp_all));
            pupil_mod_tpmeans{mouse_n}(cell_n,:,:) = repmat(max(pupil_mod{mouse_n}(cell_n,:)),length(mouse_resp_all),length(mouse_resp_all));
            behavioural_mod_tpmeans{mouse_n}(cell_n,:,:) = nanmean(cat(3,squeeze(run_mod_tpmeans{mouse_n}(cell_n,:,:)),squeeze(pupil_mod_tpmeans{mouse_n}(cell_n,:,:))),3);
        end
        dSpeed_Rmod{mouse_n} = permute(run_mod_tpmeans{mouse_n},[2,3,1]).*mean_dSpeed{mouse_n};
        dPupil_Pmod{mouse_n} = permute(pupil_mod_tpmeans{mouse_n},[2,3,1]).*mean_dPupil{mouse_n};
        dBehav_Bmod{mouse_n} = permute(behavioural_mod_tpmeans{mouse_n},[2,3,1]).*mean_dBehav{mouse_n};
    else
        mean_speed{mouse_n}        =nan(length(data{mouse_n}(1).day),1);
        mean_dSpeed{mouse_n}       =nan(length(data{mouse_n}(1).day),length(data{mouse_n}(1).day));
        mean_pupilsize{mouse_n}    =nan(length(data{mouse_n}(1).day),1);
        mean_dPupil{mouse_n}=nan(length(data{mouse_n}(1).day),length(data{mouse_n}(1).day));
        mean_dBehav{mouse_n}=nan(length(data{mouse_n}(1).day),length(data{mouse_n}(1).day));
        run_mod_tpmeans{mouse_n}  =nan(length(data{mouse_n}),length(data{mouse_n}(1).day),length(data{mouse_n}(1).day));
        pupil_mod_tpmeans{mouse_n}=nan(length(data{mouse_n}),length(data{mouse_n}(1).day),length(data{mouse_n}(1).day));
        behavioural_mod_tpmeans{mouse_n}=nan(length(data{mouse_n}),length(data{mouse_n}(1).day),length(data{mouse_n}(1).day));
        dSpeed_Rmod{mouse_n}=nan(length(data{mouse_n}(1).day),length(data{mouse_n}(1).day),length(data{mouse_n}));
        dPupil_Pmod{mouse_n}=nan(length(data{mouse_n}(1).day),length(data{mouse_n}(1).day),length(data{mouse_n}));
        dBehav_Bmod{mouse_n}=nan(length(data{mouse_n}(1).day),length(data{mouse_n}(1).day),length(data{mouse_n}));
    end
end

% decide on interval type grouping
dT_binstart = [0,1:dT_binsize:max_interval];
dT_binend = [0,dT_binsize:dT_binsize:max_interval];
if dT_binend(end)<max_interval
    dT_binend(end+1)=max_interval;
end
dT_bins = [dT_binstart; dT_binend];

for bin = 1:size(dT_bins,2)
    if bin==1
        dT_bin_mice(bin) = sum( cellfun(@(x) max(histcounts(x(1).day,'BinWidth',0.5))>1,data));
    else
        dT_bin_mice(bin) = sum(cellfun(@(x) any(ismember(x(1).comp_type_mat_days(:),dT_bins(1,bin):dT_bins(2,bin))), data));
    end
end

clear dT_bin_labels
dT_bin_labels{1} = '0';
for dT_bin = 2:size(dT_bins,2)
    if dT_binsize == 1
        dT_bin_labels{dT_bin} = [num2str(dT_bins(1,dT_bin))];
    else
        dT_bin_labels{dT_bin} = [num2str(dT_bins(1,dT_bin)) '-' num2str(dT_bins(2,dT_bin))];
    end
end

mouse_n_perint = cell(1,max_interval+1); % variable to collect mouse number for each cell changes for each time interval
dPOsig_perint = cell(1,max_interval+1); % variable to collect all PO changes in based on time interval
dPOall_perint = cell(1,max_interval+1); % variable to collect all PO changes in based on time interval
issig_perint = cell(1,max_interval+1);  % variable to collect all index of is change is sig
start_PO_perint_all = cell(1,max_interval+1);  % variable to collect initial PO of sig changes
start_PO_perint = cell(1,max_interval+1);  % variable to collect initial PO of sig changes
end_PO_perint_all = cell(1,max_interval+1);  % variable to collect initial PO of sig changes
end_PO_perint = cell(1,max_interval+1);  % variable to collect initial PO of sig changes
dPOmin_perint = cell(1,max_interval+1);  % variable to collect all min detectable dPO
both_tuned_perint = cell(1,max_interval+1);  % variable to collect all cell tuned at both tps
both_visresp_perint = cell(1,max_interval+1);  % variable to collect all cell tuned at both tps
cells_per_perint_all = cell(1,max_interval+1);
cells_per_perint_sig = cell(1,max_interval+1);
run_mod_perint = cell(1,max_interval+1);
dSpeed_perint = cell(1,max_interval+1);
dSpeed_Rmod_perint = cell(1,max_interval+1);
pupil_mod_perint = cell(1,max_interval+1);
dPupil_perint = cell(1,max_interval+1);
dPupil_Pmod_perint = cell(1,max_interval+1);
behav_mod_perint = cell(1,max_interval+1);
dBehav_perint = cell(1,max_interval+1);
dBehav_Bmod_perint = cell(1,max_interval+1);

mouse_n_perint_binned = cell(1,size(dT_bins,2));
dPOsig_perint_binned = cell(1,size(dT_bins,2)); % variable to collect all PO changes in based on time interval
dPOall_perint_binned = cell(1,size(dT_bins,2)); % variable to collect all PO changes in based on time interval
issig_perint_binned = cell(1,size(dT_bins,2));  % variable to collect all index of is change is sig (might not need)
start_PO_perint_binned_all = cell(1,size(dT_bins,2));
start_PO_perint_binned = cell(1,size(dT_bins,2));
end_PO_perint_binned_all = cell(1,size(dT_bins,2));
end_PO_perint_binned = cell(1,size(dT_bins,2));
dPOmin_perint_binned = cell(1,size(dT_bins,2));
both_tuned_perint_binned = cell(1,size(dT_bins,2));
both_visresp_perint_binned = cell(1,size(dT_bins,2));
cells_per_perint_binned = cell(1,size(dT_bins,2));
run_mod_perint_binned = cell(1,size(dT_bins,2));
dSpeed_perint_binned = cell(1,size(dT_bins,2));
pupil_mod_perint_binned = cell(1,size(dT_bins,2));
dPupil_perint_binned = cell(1,size(dT_bins,2));
behav_mod_perint_binned = cell(1,size(dT_bins,2));
dBehav_perint_binned = cell(1,size(dT_bins,2));

for mouse_n = 1:length(data)
    % sig dPO change matrix
    POdif_issig = single(cat(3,data{mouse_n}.POdif_issig));
    for i = 1:size(POdif_issig,1)
        POdif_issig(i,i,:)=nan;
    end
    
    % matrix of sig dPO magnitudes
    POdif_sig = cat(3,data{mouse_n}.POdif_sig).*-1;% the sign of the circ_dist function is oposite to as it caculates x-y not y-x
    POdif_all = cat(3,data{mouse_n}.POdif_all).*-1;% the sign of the circ_dist function is oposite to as it caculates x-y not y-x
    POdif_mouse_n = ones(size(POdif_all))*mouse_n;
    
    % interval type matrix
    int_mat = data{mouse_n}(1).comp_type_mat_days;
    
    % min detectable dPO magnitude & cell is tuned at both tps
    clear start_PO_all end_PO_all dPO_min_min both_tunned both_visresp start_PO_CI_all end_PO_CI_all
    for celln = 1:length(data{mouse_n})
        PO_all = data{mouse_n}(celln).PO;
        PO_CI_all = cellfun(@(x) abs(circ_dist(deg2rad(x(2)),deg2rad(x(1)))),cat(2,data{mouse_n}(celln).PO_CI_norm)); % CI size
        for tp1 = 1:length(PO_all)
            for tp2 = 1:length(PO_all)
                tp_start = min([tp1,tp2]); % makes matrix symetrical, change direction is calcualted consistently
                tp_end = max([tp1,tp2]);
                start_PO_all(tp1,tp2,celln) = PO_all(tp_start);
                end_PO_all(tp1,tp2,celln) = PO_all(tp_end);
                start_PO_CI_all(tp1,tp2,celln) = PO_CI_all(tp_start);
                end_PO_CI_all(tp1,tp2,celln) = PO_CI_all(tp_end);
            end
        end
        
        % matrix of min detectable dPO magnitude
        PO_CI_all_temp = data{mouse_n}(celln).PO_CI_norm;
        PO_CI_all_temp = cellfun(@(x) wrapToPi(deg2rad(x)),PO_CI_all_temp,'UniformOutput',0);
        dPO_min_all = cellfun(@(x) rad2deg(mean(abs(x))), PO_CI_all_temp);
        for tp1 = 1:length(dPO_min_all)
            for tp2 = 1:length(dPO_min_all)
                dPO_min_min(tp1,tp2,celln) = max(dPO_min_all([tp1,tp2]),[],'includenan'); % in this case if one is a nan there should not be a min dPO
            end
        end
        
        % matrix of cell is tunned at both tps
        stat_tuned = data{mouse_n}(celln).Stat_tuned;
        for tp1 = 1:length(dPO_min_all)
            for tp2 = 1:length(dPO_min_all)
                both_tunned(tp1,tp2,celln) = all(stat_tuned([tp1,tp2]));
            end
        end
        
        % matrix of cell is tunned at both tps
        stat_visresp = data{mouse_n}(celln).visually_responsive;
        for tp1 = 1:length(dPO_min_all)
            for tp2 = 1:length(dPO_min_all)
                both_visresp(tp1,tp2,celln) = all(stat_visresp([tp1,tp2]));
            end
        end
        
    end
    
    start_PO=start_PO_all;
    end_PO=end_PO_all;
    start_PO_CI_sig = start_PO_CI_all;
    end_PO_CI_sig = end_PO_CI_all;
    
    % only take significant changes
    start_PO(find(POdif_issig~=1))=nan; % sig only
    end_PO(find(POdif_issig~=1))=nan; % sig only
    start_PO_CI_sig(find(POdif_issig~=1))=nan; % sig only
    end_PO_CI_sig(find(POdif_issig~=1))=nan; % sig only
    
    for dT = 1:max_interval+1
        idx_dT = int_mat==dT-1; % find interval type (first dT looks for same day comparisons, hence the -1)
        idx_dT = tril(idx_dT,-1); % lower half excluding diag
        [x,y] = find(idx_dT);
        
        % cells id mat
        cell_idx_mat_all = zeros(size(POdif_all));
        for i = 1:size(cell_idx_mat_all,3); cell_idx_mat_all(:,:,i) = mouse_n*10000+i;end
        cell_idx_mat_sig = cell_idx_mat_all;
        cell_idx_mat_all(isnan(POdif_all))=nan;
        cell_idx_mat_sig(isnan(POdif_sig))=nan;
        
        for i=1:length(x)
            % merge sig change
            temp_POdif_mouse_n = squeeze(POdif_mouse_n(x(i),y(i),:));
            mouse_n_perint{dT} = cat(1,mouse_n_perint{dT},temp_POdif_mouse_n);
            
            % merge sig change
            temp_sigdPOs = squeeze(POdif_sig(x(i),y(i),:));
            dPOsig_perint{dT} = cat(1,dPOsig_perint{dT},temp_sigdPOs);
            
            % merge all change
            temp_alldPOs = squeeze(POdif_all(x(i),y(i),:));
            dPOall_perint{dT} = cat(1,dPOall_perint{dT},temp_alldPOs);
            
            % merge sig dPO magnitudes
            temp_issigs = squeeze(POdif_issig(x(i),y(i),:)); % collect indicese from interval length
            issig_perint{dT} = cat(1,issig_perint{dT},temp_issigs);
            
            % merge start PO of sig changes
            temp_start_PO = squeeze(start_PO_all(x(i),y(i),:)); % collect indicese from interval length
            start_PO_perint_all{dT} = cat(1,start_PO_perint_all{dT},temp_start_PO);
            
            % merge start PO of sig changes
            temp_start_PO = squeeze(start_PO(x(i),y(i),:)); % collect indicese from interval length
            start_PO_perint{dT} = cat(1,start_PO_perint{dT},temp_start_PO);
            
            % merge end PO of sig changes
            temp_end_PO = squeeze(end_PO_all(x(i),y(i),:)); % collect indicese from interval length
            end_PO_perint_all{dT} = cat(1,end_PO_perint_all{dT},temp_end_PO);
            
            % merge end PO of sig changes
            temp_end_PO = squeeze(end_PO(x(i),y(i),:)); % collect indicese from interval length
            end_PO_perint{dT} = cat(1,end_PO_perint{dT},temp_end_PO);
            
            % merge min detectable dPO
            temp_dPOmin = squeeze(dPO_min_min(x(i),y(i),:));
            dPOmin_perint{dT} = cat(1,dPOmin_perint{dT},temp_dPOmin);
            
            % comparisons where cells are tunned at both tps
            temp_both_tuned = squeeze(both_tunned(x(i),y(i),:));
            both_tuned_perint{dT} = cat(1,both_tuned_perint{dT},temp_both_tuned);
            
            % comparisons where cells are visually responsive at both tps
            temp_both_visresp = squeeze(both_visresp(x(i),y(i),:));
            both_visresp_perint{dT} = cat(1,both_visresp_perint{dT},temp_both_visresp);
            
            % collect cell ids (not true cell id)
            cells_per_perint_all{dT} = cat(1,cells_per_perint_all{dT},squeeze(cell_idx_mat_all(x(i),y(i),:)));
            cells_per_perint_sig{dT} = cat(1,cells_per_perint_sig{dT},squeeze(cell_idx_mat_sig(x(i),y(i),:)));
            
            % merge running vars
            temp_run_mod = squeeze(run_mod_tpmeans{mouse_n}(:,x(i),y(i))); % collect indicese from interval length
            temp_dSpeed = repmat(mean_dSpeed{mouse_n}(x(i),y(i)),size(start_PO,3),1); % collect indicese from interval length
            temp_dSpeed_Rmod = squeeze(dSpeed_Rmod{mouse_n}(x(i),y(i),:)); % collect indicese from interval length
            run_mod_perint{dT} = cat(1,run_mod_perint{dT},temp_run_mod);
            dSpeed_perint{dT} = cat(1,dSpeed_perint{dT},temp_dSpeed);
            dSpeed_Rmod_perint{dT} = cat(1,dSpeed_Rmod_perint{dT},temp_dSpeed_Rmod);
            
            % merge pupil vars
            temp_pupil_mod = squeeze(pupil_mod_tpmeans{mouse_n}(:,x(i),y(i))); % collect indicese from interval length
            temp_dPupil = repmat(mean_dPupil{mouse_n}(x(i),y(i)),size(start_PO,3),1); % collect indicese from interval length
            temp_dPupil_Pmod = squeeze(dPupil_Pmod{mouse_n}(x(i),y(i),:)); % collect indicese from interval length
            pupil_mod_perint{dT} = cat(1,pupil_mod_perint{dT},temp_pupil_mod);
            dPupil_perint{dT} = cat(1,dPupil_perint{dT},temp_dPupil);
            dPupil_Pmod_perint{dT} = cat(1,dPupil_Pmod_perint{dT},temp_dPupil_Pmod);
            
            % merge behav vars
            temp_behav_mod = squeeze(behavioural_mod_tpmeans{mouse_n}(:,x(i),y(i))); % collect indicese from interval length
            temp_dBehav = repmat(mean_dBehav{mouse_n}(x(i),y(i)),size(start_PO,3),1); % collect indicese from interval length
            temp_dBehav_Bmod = squeeze(dBehav_Bmod{mouse_n}(x(i),y(i),:)); % collect indicese from interval length
            behav_mod_perint{dT} = cat(1,behav_mod_perint{dT},temp_behav_mod);
            dBehav_perint{dT} = cat(1,dBehav_perint{dT},temp_dBehav);
            dBehav_Bmod_perint{dT} = cat(1,dBehav_Bmod_perint{dT},temp_dBehav_Bmod);            
        end
    end
end

% combine multiple time intervals into bins
for dT_bin = 1:size(dT_bins,2)
    mouse_n_perint_binned{dT_bin}        = cat(1,mouse_n_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1});
    dPOsig_perint_binned{dT_bin}        = cat(1,dPOsig_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1});
    dPOall_perint_binned{dT_bin}        = cat(1,dPOall_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1});
    issig_perint_binned{dT_bin}         = cat(1,issig_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1});
    start_PO_perint_binned_all{dT_bin}  = cat(1,start_PO_perint_all{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1});
    start_PO_perint_binned{dT_bin}      = cat(1,start_PO_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1});
    end_PO_perint_binned_all{dT_bin}    = cat(1,end_PO_perint_all{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1});
    end_PO_perint_binned{dT_bin}        = cat(1,end_PO_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1});
    dPOmin_perint_binned{dT_bin}        = cat(1,dPOmin_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1});
    both_tuned_perint_binned{dT_bin}    = cat(1,both_tuned_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1});
    both_visresp_perint_binned{dT_bin}  = cat(1,both_visresp_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1});
    cells_per_perint_binned{dT_bin}  = cat(1,cells_per_perint_all {[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1});
    
    run_mod_perint_binned{dT_bin} = cat(1,run_mod_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1}); % running modulation index for cell (mean across tps)
    dSpeed_perint_binned{dT_bin} = cat(1,dSpeed_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1}); % mean speed diff across sessions
    pupil_mod_perint_binned{dT_bin} = cat(1,pupil_mod_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1}); % running modulation index for cell (mean across tps)
    dPupil_perint_binned{dT_bin} = cat(1,dPupil_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1}); % mean speed diff across sessions
    behav_mod_perint_binned{dT_bin} = cat(1,behav_mod_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1}); % running modulation index for cell (mean across tps)
    dBehav_perint_binned{dT_bin} = cat(1,dBehav_perint{[dT_bins(1,dT_bin):dT_bins(2,dT_bin)]+1}); % mean speed diff across sessions
    
    if isempty(issig_perint_binned{dT_bin})
        mouse_n_perint_binned{dT_bin} = [];
        dPOsig_perint_binned{dT_bin} = [];
        dPOall_perint_binned{dT_bin} = [];
        issig_perint_binned{dT_bin} = [];
        start_PO_perint_binned_all{dT_bin} = [];
        start_PO_perint_binned{dT_bin} = [];
        end_PO_perint_binned_all{dT_bin} = [];
        end_PO_perint_binned{dT_bin} = [];
        dPOmin_perint_binned{dT_bin} = [];
        both_tuned_perint_binned{dT_bin} = [];
        both_visresp_perint_binned{dT_bin} = [];
        cells_per_perint_binned{dT_bin}  = [];
        
        run_mod_perint_binned{dT_bin} = []; % running modulation index for cell (mean across tps)
        dSpeed_perint_binned{dT_bin} = []; % mean speed diff across sessions
        pupil_mod_perint_binned{dT_bin} = []; % running modulation index for cell (mean across tps)
        dPupil_perint_binned{dT_bin} = []; % mean speed diff across sessions
        behav_mod_perint_binned{dT_bin} = []; % running modulation index for cell (mean across tps)
        dBehav_perint_binned{dT_bin} = []; % mean speed diff across sessions
    end
end

both_visrespANDtuned_perint_binned = cellfun(@(a,b) a.*b, both_tuned_perint_binned,both_visresp_perint_binned,'UniformOutput', false);

% calculate fractions
numb_sig = cellfun(@(x) sum(~isnan(x)),dPOsig_perint_binned,'UniformOutput', false); numb_sig = cat(2,numb_sig{:});
numb_tuned = cellfun(@(x) sum(~isnan(x)),dPOall_perint_binned,'UniformOutput', false); numb_tuned = cat(2,numb_tuned{:}); % only considering those that are conc tuned and vis resp (othewise they have no PO to compare)
numb_all = cellfun(@(x) length(x),dPOall_perint_binned,'UniformOutput', false); numb_all = cat(2,numb_all{:});
frac_sig_all = numb_sig./numb_all;
frac_sig_bothtuned = numb_sig./numb_tuned;
frac_both_tuned = numb_tuned./numb_all;

numb_both_tuned = cellfun(@(x) sum(x),both_visresp_perint_binned,'UniformOutput', false); numb_both_tuned = cat(2,numb_both_tuned{:});

% combine into 1D variables for boxplot/scatter
dPOsig_vect = [];
dPOall_vect = [];
dT_bin_dPOsig_vect = [];
dT_bin_dPOall_vect = [];
mindPO_vect = [];
dT_bin_mindPO_vect = [];
for dT_bin = 1:size(dT_bins,2)
    if isempty(dPOsig_perint_binned{dT_bin})
        dPOsig_vect = cat(1,dPOsig_vect,nan);
        dT_bin_dPOsig_vect = cat(1,dT_bin_dPOsig_vect,dT_bin);
    else
        dPOsig_vect = cat(1,dPOsig_vect,dPOsig_perint_binned{dT_bin});
        dT_bin_dPOsig_vect = cat(1,dT_bin_dPOsig_vect,ones(length(dPOsig_perint_binned{dT_bin}),1)*dT_bin);
    end
    
    if isempty(dPOall_perint_binned{dT_bin})
        dPOall_vect = cat(1,dPOall_vect,nan);
        dT_bin_dPOall_vect = cat(1,dT_bin_dPOall_vect,dT_bin);
    else
        dPOall_vect = cat(1,dPOall_vect,dPOall_perint_binned{dT_bin});
        dT_bin_dPOall_vect = cat(1,dT_bin_dPOall_vect,ones(length(dPOall_perint_binned{dT_bin}),1)*dT_bin);
    end
    
    if isempty(dPOmin_perint_binned{dT_bin})
        mindPO_vect = cat(1,mindPO_vect,nan);
    else
        mindPO_vect = cat(1,mindPO_vect,dPOmin_perint_binned{dT_bin});
    end
end

% combine to make 2D hist of starting PO and dT
PO_bins = -90:30:90;
clear PO_bin_labels
for PO_bin = 1:size(PO_bins,2)-1
    PO_bin_labels{PO_bin} = [num2str(PO_bins(PO_bin)) '-' num2str(PO_bins(PO_bin+1))];
end
for dT_bin = 1:size(dT_bins,2)
    PO_idx{dT_bin} = discretize(start_PO_perint_binned{dT_bin},PO_bins); % PO bin index
    for idx=1:length(PO_bins)-1
        % fraction sig dPO
        temp =issig_perint_binned{dT_bin}(find(PO_idx{dT_bin}==idx));
        PO_bin_frac_sig(dT_bin,idx) = sum(temp)/length(temp);
        
        % magnitudes of sig dPO
        dPOsig_dT_PO{dT_bin,idx} = dPOsig_perint_binned{dT_bin}(find(PO_idx{dT_bin}==idx));
        
        % min detectable dPO
        dPOmin_dT_PO{dT_bin,idx} = dPOmin_perint_binned{dT_bin}(find(PO_idx{dT_bin}==idx));
        
        % fraction concurently tuned comparisons (cells are tuned for both tps)
        temp = both_tuned_perint_binned{dT_bin}(find(PO_idx{dT_bin}==idx));
        both_tuned_dT_PO(dT_bin,idx) = sum(temp)/length(temp);
    end
end

ratio_ = @(x) nanmedian(x);

for dT_bin = 1:size(dT_bins,2)
    theta_start_all = wrapTo180(([start_PO_perint_binned_all{dT_bin}]')).*2;
    theta_end_all = wrapTo180(([end_PO_perint_binned_all{dT_bin}]')).*2;
    theta_end_all(isnan(theta_start_all(:)))=[];
    theta_start_all(isnan(theta_start_all(:)))=[];
    
    if ~isempty(PO_relative)
        theta_start_all = wrapTo180([theta_start_all-PO_relative.*2]);
        theta_end_all = wrapTo180([theta_end_all-PO_relative.*2]);
    end
    
    theta_change_all = abs(theta_start_all)-abs(theta_end_all); % testing
    realtive_theta_change_all{dT_bin} = theta_change_all;
    
    if length(realtive_theta_change_all{dT_bin})>3
        ratio_boot = bootstrp(1000,ratio_,realtive_theta_change_all{dT_bin}./2);
        ratio_conTOdiv{dT_bin} = ratio_(realtive_theta_change_all{dT_bin}./2);
        ratio_conTOdiv_CI{dT_bin} = [prctile(ratio_boot,5); prctile(ratio_boot,95)];
    else
        ratio_conTOdiv{dT_bin} = nan;
        ratio_conTOdiv_CI{dT_bin} = [nan nan];
    end
end

% output for making all possible change plots
for dT_bin = 1:size(dT_bins,2)
    if ~isempty(start_PO_perint_binned_all{dT_bin})
        output{dT_bin}(1,:)=start_PO_perint_binned_all{dT_bin}; % this is normalized to permitted unless input specifies otherwise
        output{dT_bin}(2,:)=end_PO_perint_binned_all{dT_bin};
        output{dT_bin}(3,:)=dPOall_perint_binned{dT_bin}; % change relative to true starting PO (not relative to permitted ori)
        output{dT_bin}(4,:)=issig_perint_binned{dT_bin};
        
        output{dT_bin}(5,:)= run_mod_perint_binned{dT_bin};% behavioural modulation index for cell (mean across tps)
        output{dT_bin}(6,:)= dSpeed_perint_binned{dT_bin};% mean behaviour (running and pupilsize) diff across sessions
        output{dT_bin}(8,:)= pupil_mod_perint_binned{dT_bin};% behavioural modulation index for cell (mean across tps)
        output{dT_bin}(9,:)= dPupil_perint_binned{dT_bin};% mean behaviour (running and pupilsize) diff across sessions

        output{dT_bin}(7,:)=mouse_n_perint_binned{dT_bin};
        
        temp = isnan(output{dT_bin});
        hasnan = any(temp(1:2,:),1);
        
        output{dT_bin}(:,hasnan)=[];
    else
        output{dT_bin} = [];
    end
end

end


function [PO_dif, PO_issig, PO_sig] = PO_diff_matrix(data, showplots)
%% PO from vector average with significance based on bootstrapped 95% conf
PO_all = wrapToPi(deg2rad([data.PO]));
PO_CI_lb = wrapToPi(deg2rad(cell2mat(cellfun(@(x) x(1),data.PO_CI,'UniformOutput', false))));
PO_CI_ub = wrapToPi(deg2rad(cell2mat(cellfun(@(x) x(2),data.PO_CI,'UniformOutput', false))));
PO_tuned = [data.Stat_tuned];
is_resp = cat(1,data.visually_responsive)';

clear PO_issig PO_dif
for tpa = 1:length(data.day)
    for tpb = 1:length(data.day)
        [tp1]=min([tpa,tpb]);
        [tp2]=max([tpa,tpb]);
        PO_dif_rad = wrapToPi(circ_dist(PO_all(tp1).*2,PO_all(tp2).*2)./2);
        
        if any(~PO_tuned([tp1 tp2])) || any(~is_resp([tp1 tp2]))
            PO_issig(tpa,tpb) = 0; %  if ether of the timepoints is not tuned or visually responseive it is considered not sig
        else   
            PO_issig(tpa,tpb) = ~all(circ_dist([PO_all(tp1),PO_all(tp1)].*2,[PO_CI_ub(tp2) PO_CI_lb(tp2)].*2)>0==[0 1]) && ...
                ~all(circ_dist([PO_all(tp2),PO_all(tp2)].*2,[PO_CI_ub(tp1) PO_CI_lb(tp1)].*2)>0==[0 1]);
        end
        
        PO_dif(tpa,tpb) = wrapTo180(rad2deg(PO_dif_rad));

        if showplots
            figure(105); hold off;
            polarplot([PO_all(tpa) PO_all(tpa)],[0 1],'b'); hold on;
            polarplot(PO_all(tpa)+[PO_CI_lb(tpa) PO_CI_lb(tpa)],rlim,'--b','LineWidth',2);
            polarplot(PO_all(tpa)+[PO_CI_ub(tpa) PO_CI_ub(tpa)],rlim,'--b');
            polarplot([PO_all(tpb) PO_all(tpb)],rlim,'r');
            polarplot(PO_all(tpb)+[PO_CI_lb(tpb) PO_CI_lb(tpb)],rlim,'--r','LineWidth',2);
            polarplot(PO_all(tpb)+[PO_CI_ub(tpb) PO_CI_ub(tpb)],rlim,'--r');
            title({['angle dif: ' num2str(wrapTo180(rad2deg(PO_dif_rad)),2), ', sig: ' num2str(PO_issig(tpa,tpb),2)]})
            thetaticks([0:60:300])
            thetaticklabels({'0' '30' '60' '+/-90' '-60' '-30'})
            set(gcf,'color','w');
            set(gca,'FontSize',15);
            pause(0.5)
        end
    end
end

PO_sig = PO_dif; PO_sig(~PO_issig) = nan;

if showplots
    figure; set(gcf,'color','w')
    subplot(1,2,1);
    histogram(PO_sig,[-90:5:90]); hold on
    histogram(PO_dif(logical(~PO_issig)),[-90:5:90]); hold on;
    title(['% sig changes: ' num2str(100*sum(PO_issig(:))/length(PO_issig(:)),2)])
    set(gca,'TickDir','out')
    xlabel('angle')
    ylabel('count')
    axis square
    
    subplot(1,2,2);
    polarhistogram(deg2rad((PO_dif(logical(PO_issig)))).*2,'BinWidth',pi/72,'Normalization','probability');
    thetaticks([0:60:300])
    thetaticklabels({'0' '30' '60' '+/-90' '-60' '-30'})
end

end
