function [PO_to_PrefO_count,...
    bin_centers,...
    PO_bin_frac_sig_ofTuned,...
    comp_type_label,...
    ratio_conTOdiv,...
    ratio_conTOdiv_CI,...
    dPOsig_perint,...
    dPOsig_perint_all,...
    output]...
    = StripeRearing_stability(data,timepoint_labels,varargin)
p = inputParser();
p.KeepUnmatched = false;
p.CaseSensitive = false;
p.StructExpand  = false;
p.PartialMatching = true;

addParameter(p, 'cell_filter' , 'resp_anytime') 
addParameter(p, 'comparison_type' , 'Sequential') 
addParameter(p, 'Binning_interval' , 30) 
addParameter(p, 'use_only_sig_changes', 0)
addParameter(p, 'permitted_ori_override',[])

parse(p, varargin{:})
names = fieldnames(p.Results);
for i=1:length(names)
    eval([names{i} '=p.Results.' names{i} ]);
end

data(cellfun(@isempty, data))= [];

PO_bins = (-(90+Binning_interval):Binning_interval:90) + (Binning_interval/2); % bin size 45° (rebinned)
bin_centers = -90:Binning_interval:90-Binning_interval;

% subselect cells based on filter
if ischar(cell_filter)
    if strcmp(cell_filter,'resp_anytime')
        for mouse_n = 1:length(data)
            idx_filter = cellfun(@(x) any(x),{data{mouse_n}.visually_responsive});
            data{mouse_n} = data{mouse_n}(find(idx_filter));
        end
    elseif contains(cell_filter,'resp_and_tuned_from_BL2toSR')
        idx_SR = str2num(cell_filter(end))+2;
        for mouse_n = 1:length(data)
            idx_filter1 = cellfun(@(x) all(x(2:idx_SR)),{data{mouse_n}.visually_responsive});
            idx_filter2 = cellfun(@(x) all(x(2:idx_SR)),{data{mouse_n}.Stat_tuned});
            data{mouse_n} = data{mouse_n}(find(idx_filter1&idx_filter2));
        end
    elseif contains(cell_filter,'resp_and_tuned_at_BL2xorSR')
        idx_SR = str2num(cell_filter(end))+2;
        for mouse_n = 1:length(data)
            idx_filtera = cellfun(@(x) logical(x(2)),{data{mouse_n}.Stat_tuned});
            idx_filterb = cellfun(@(x) x(2),{data{mouse_n}.visually_responsive});
            idx_filterc = cellfun(@(x) ~x(idx_SR),{data{mouse_n}.Stat_tuned});
            idx_filterd = cellfun(@(x) ~x(idx_SR),{data{mouse_n}.visually_responsive});
            idx_filter1 = (idx_filtera&idx_filterb) & (idx_filterc|idx_filterd);
            
            idx_filtera = cellfun(@(x) logical(x(idx_SR)),{data{mouse_n}.Stat_tuned});
            idx_filterb = cellfun(@(x) x(idx_SR),{data{mouse_n}.visually_responsive});
            idx_filterc = cellfun(@(x) ~x(2),{data{mouse_n}.Stat_tuned});
            idx_filterd = cellfun(@(x) ~x(2),{data{mouse_n}.visually_responsive});
            idx_filter2 = (idx_filtera&idx_filterb) & (idx_filterc|idx_filterd);
            
            data{mouse_n} = data{mouse_n}(find((idx_filter1)|(idx_filter2)));
        end
    else
        error('invalid cell filter')
    end
else
    for mouse_n = 1:length(data)
        idx_filter = find(ismember([data{mouse_n}.ROI_group],cell_filter));
        data{mouse_n} = data{mouse_n}(idx_filter);
    end
end

data = data(cellfun(@(x) ~isempty(x),data)); % remove mice with no cells left

% decide on interval type grouping
if strcmp(comparison_type,'Sequential')
    temp = diag(1:length(timepoint_labels)-1,-1)+diag(1:length(timepoint_labels)-1,1);
    temp(temp==0)=nan;
    comp_type = temp;
    
    for i = 1:max(comp_type(:))
        [x,y] = find(comp_type==i,1);
        comp_type_label{i} = [timepoint_labels{min(x,y)} '-' timepoint_labels{max(x,y)}];
    end
elseif strcmp(comparison_type,'BL2-All')
    BL2_idx = find(cellfun(@(x) strcmp(x,'BL2'),timepoint_labels));
    temp = zeros(length(timepoint_labels));
    counter = 1;
    for i = 1:length(timepoint_labels)
        for ii = 1:length(timepoint_labels)
            if i ~=ii && any([i,ii]==BL2_idx) && i>ii
                temp(i,ii)=counter;counter = counter+1;
            end
        end
    end
    temp = temp.'+temp; temp(temp==0)=nan;
    comp_type = temp;
    
    for i = 1:max(comp_type(:))
        [x,y] = find(comp_type==i,1);
        comp_type_label{i} = [timepoint_labels{min(x,y)} '-' timepoint_labels{max(x,y)}];
    end
end

max_tps = length(timepoint_labels);


issig_perint = cell(1,sum(~isnan(unique(comp_type))));  % variable to collect all index of is change is sig
dPOsig_perint = cell(1,sum(~isnan(unique(comp_type)))); % variable to collect all PO changes in based on time interval
start_PO_perint = cell(1,sum(~isnan(unique(comp_type))));  % variable to collect initial PO of sig changes
end_PO_perint = cell(1,sum(~isnan(unique(comp_type))));  % variable to collect initial PO of sig changes
dPOsig_perint_all = cell(1,sum(~isnan(unique(comp_type))));
dPOCI_dPOsig_perint_all = cell(1,sum(~isnan(unique(comp_type))));
start_PO_perint_all = cell(1,sum(~isnan(unique(comp_type))));  % variable to collect initial PO of sig changes
end_PO_perint_all = cell(1,sum(~isnan(unique(comp_type))));  % variable to collect initial PO of sig changes
dPOmin_perint = cell(1,sum(~isnan(unique(comp_type))));  % variable to collect all min detectable dPO
both_tuned_perint = cell(1,sum(~isnan(unique(comp_type))));  % variable to collect all cell tuned at both tps
both_visresp_perint = cell(1,sum(~isnan(unique(comp_type))));  % variable to collect all cell vis resp at both tps
dPOCI_dPOsig_perint = cell(1,sum(~isnan(unique(comp_type)))); 
start_PO_CI_perint = cell(1,sum(~isnan(unique(comp_type))));
end_PO_CI_perint = cell(1,sum(~isnan(unique(comp_type))));  % variable to collect initial PO of sig changes

cells_per_perint_all = cell(1,sum(~isnan(unique(comp_type)))); 
cells_per_perint_sig = cell(1,sum(~isnan(unique(comp_type)))); 

PO_tuned = [];
PO_CI_tuned = [];
PO_true_tuned = [];
OSI_tuned = [];
Max_Stim_resp_change = [];
ST_all = [];
VR_all = [];
for mouse_n = 1:length(data)
    
    % distribution of tuned PO
    if isfield(data{mouse_n},'permitted_ori') & isempty(permitted_ori_override)
        PO = wrapTo180((single(cat(2,data{mouse_n}.PO))-data{mouse_n}(1).permitted_ori).*2)./2;
    elseif ~isempty(permitted_ori_override)
        PO = wrapTo180((single(cat(2,data{mouse_n}.PO))-permitted_ori_override).*2)./2;
    else
        PO = wrapTo180((single(cat(2,data{mouse_n}.PO))).*2)./2;
        disp('not centering on permitted ori')
    end
    PO_CI = cellfun(@(x) rad2deg(abs(circ_dist(deg2rad(x(2)),deg2rad(x(1))))),cat(2,data{mouse_n}.PO_CI_norm)); % CI size
    PO_true = single(cat(2,data{mouse_n}.PO));
    OSI = single(cat(2,data{mouse_n}.OSnn));
    ST = single(cat(2,data{mouse_n}.Stat_tuned));
    VR = single(cat(1,data{mouse_n}.visually_responsive))';
        
    PO_tuned = cat(2,PO_tuned,PO(1:max_tps,:)); % obv this will not work if mice have different timepoints but for now its fine;
    PO_CI_tuned = cat(2,PO_CI_tuned,PO_CI(1:max_tps,:));
    PO_true_tuned = cat(2,PO_true_tuned,PO_true(1:max_tps,:)); % obv this will not work if mice have different timepoints but for now its fine;
    OSI_tuned = cat(2,OSI_tuned,OSI(1:max_tps,:));
    ST_all = cat(2,ST_all,ST(1:max_tps,:));
    VR_all = cat(2,VR_all,VR(1:max_tps,:));
    
    % sig dPO change matrix
    POdif_issig = single(cat(3,data{mouse_n}.POdif_issig));
    for i = 1:size(POdif_issig,1)
        POdif_issig(i,i,:)=nan;
    end
    
    % matrix of PO_CI changes
    dPOCI_POdif_sig = nan(max_tps,max_tps,length(PO_CI));
    dPOCI_POdif_all = nan(max_tps,max_tps,length(PO_CI));
    for tp1 = 1:max_tps
        for tp2 = 1:max_tps
            for celln = 1:size(PO_CI,2)
                if POdif_issig(tp1,tp2,celln)==1
                    if tp1<=tp2
                        dPOCI_POdif_sig(tp1,tp2,celln) = PO_CI(tp2,celln)-PO_CI(tp1,celln);
                    else
                        dPOCI_POdif_sig(tp1,tp2,celln) = PO_CI(tp1,celln)-PO_CI(tp2,celln);
                    end
                end
                if tp1<=tp2
                    dPOCI_POdif_all(tp1,tp2,celln) = PO_CI(tp2,celln)-PO_CI(tp1,celln);
                else
                    dPOCI_POdif_all(tp1,tp2,celln) = PO_CI(tp1,celln)-PO_CI(tp2,celln);
                end
            end
        end
    end
    
    % matrix of max stimulus response
    temp = cellfun(@(x) cellfun(@(y) max(y(:)),x(1:max_tps)),{data{mouse_n}.TC_mean},'UniformOutput',0);
    Max_Stim_resp = cat(2,temp{:});
    Max_Stim_resp_change_temp = real(log2(Max_Stim_resp./Max_Stim_resp(1,:)));
    Max_Stim_resp_change = cat(2,Max_Stim_resp_change,Max_Stim_resp_change_temp);
        
    % matrix of sig dPO magnitudes
    POdif_sig = cat(3,data{mouse_n}.POdif_sig).*-1;% the sign of the circ_dist function is oposite to as its caculates x-y not y-x
    POdif_all = cat(3,data{mouse_n}.POdif_all).*-1;% the sign of the circ_dist function is oposite to as its caculates x-y not y-x
        
    % min detectable dPO magnitude & cell is tuned at both tps
    
    clear start_PO_all end_PO_all dPO_min_min both_tunned both_visresp start_PO_CI_all end_PO_CI_all
    for celln = 1:length(data{mouse_n})
        if isfield(data{mouse_n},'permitted_ori') & isempty(permitted_ori_override)
            PO_all = wrapTo180((data{mouse_n}(celln).PO-data{mouse_n}(1).permitted_ori).*2)./2;
        elseif ~isempty(permitted_ori_override)
            PO_all = wrapTo180((data{mouse_n}(celln).PO-permitted_ori_override).*2)./2;
        else
            PO_all = wrapTo180((data{mouse_n}(celln).PO).*2)./2;
        end
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
    
    start_PO_sig = start_PO_all;
    end_PO_sig = end_PO_all;
    start_PO_CI_sig = start_PO_CI_all;
    end_PO_CI_sig = end_PO_CI_all;
    
    % only take significant changes
    start_PO_sig(find(POdif_issig~=1))=nan; % sig only
    end_PO_sig(find(POdif_issig~=1))=nan; % sig only
    start_PO_CI_sig(find(POdif_issig~=1))=nan; % sig only
    end_PO_CI_sig(find(POdif_issig~=1))=nan; % sig only
    
    for IT = 1:length(dPO_min_all)-1
        idx_IT = comp_type==IT;
        idx_IT = tril(idx_IT,-1);% lower half excluding diag
        [x,y] = find(idx_IT);
        
        % cells id mat
        cell_idx_mat_all = zeros(size(POdif_all));
        for i = 1:size(cell_idx_mat_all,3); cell_idx_mat_all(:,:,i) = mouse_n*10000+i;end
        cell_idx_mat_sig = cell_idx_mat_all;
        cell_idx_mat_all(isnan(POdif_all))=nan;
        cell_idx_mat_sig(isnan(POdif_sig))=nan;
        
        for i=1:length(x)
            % merge is sig change
            temp_issigs = squeeze(POdif_issig(x(i),y(i),:)); % collect indicese from interval length
            issig_perint{IT} = cat(1,issig_perint{IT},temp_issigs);
            
            % merge sig dPO magnitudes
            temp_sigdPOs = squeeze(POdif_sig(x(i),y(i),:));
            dPOsig_perint{IT} = cat(1,dPOsig_perint{IT},temp_sigdPOs);
            
            % merge PO_CI diff of sig dPOs
            temp = squeeze(dPOCI_POdif_sig(x(i),y(i),:));
            dPOCI_dPOsig_perint{IT} = cat(1,dPOCI_dPOsig_perint{IT},temp);
            
            % merge start PO of sig changes
            temp_start_PO = squeeze(start_PO_sig(x(i),y(i),:)); % collect indicese from interval length
            start_PO_perint{IT} = cat(1,start_PO_perint{IT},temp_start_PO);
            
            % merge end PO of sig changes
            temp_end_PO = squeeze(end_PO_sig(x(i),y(i),:)); % collect indicese from interval length
            end_PO_perint{IT} = cat(1,end_PO_perint{IT},temp_end_PO);
            
            % merge start PO_CI of sig changes
            temp_start_PO = squeeze(start_PO_CI_sig(x(i),y(i),:)); % collect indicese from interval length
            start_PO_CI_perint{IT} = cat(1,start_PO_CI_perint{IT},temp_start_PO);
            
            % merge end PO_CI of sig changes
            temp_end_PO = squeeze(end_PO_CI_sig(x(i),y(i),:)); % collect indicese from interval length
            end_PO_CI_perint{IT} = cat(1,end_PO_CI_perint{IT},temp_end_PO);
            
            % merge all dPO magnitudes
            temp_sigdPOs_all = squeeze(POdif_all(x(i),y(i),:));
            dPOsig_perint_all{IT} = cat(1,dPOsig_perint_all{IT},temp_sigdPOs_all);
            
            % merge all dPO_CI magnitudes
            temp_sigddPOCI_POs_all = squeeze(dPOCI_POdif_all(x(i),y(i),:));
            dPOCI_dPOsig_perint_all{IT} = cat(1,dPOCI_dPOsig_perint_all{IT},temp_sigddPOCI_POs_all);
            
            % merge start PO of all changes
            temp_start_PO_all = squeeze(start_PO_all(x(i),y(i),:)); % collect indicese from interval length
            start_PO_perint_all{IT} = cat(1,start_PO_perint_all{IT},temp_start_PO_all);
            
            % merge end PO of all changes
            temp_end_PO_all = squeeze(end_PO_all(x(i),y(i),:)); % collect indicese from interval length
            end_PO_perint_all{IT} = cat(1,end_PO_perint_all{IT},temp_end_PO_all);
            
            % merge min detectable dPO
            temp_dPOmin = squeeze(dPO_min_min(x(i),y(i),:));
            dPOmin_perint{IT} = cat(1,dPOmin_perint{IT},temp_dPOmin);
            
            % comparisons where cells are tunned at both tps
            temp_both_tuned = squeeze(both_tunned(x(i),y(i),:));
            both_tuned_perint{IT} = cat(1,both_tuned_perint{IT},temp_both_tuned);
            
            % comparisons where cells are visually responsive at both tps
            temp_both_visresp = squeeze(both_visresp(x(i),y(i),:));
            both_visresp_perint{IT} = cat(1,both_visresp_perint{IT},temp_both_visresp);
            
            % collect cell ids (not true cell id)
            cells_per_perint_all{IT} = cat(1,cells_per_perint_all{IT},squeeze(cell_idx_mat_all(x(i),y(i),:)));
            cells_per_perint_sig{IT} = cat(1,cells_per_perint_sig{IT},squeeze(cell_idx_mat_sig(x(i),y(i),:)));
        end
    end
end

%% binning is completely different in this version
% combine into 1D variables for boxplot/scatter
dPOsig_vect = [];
IT_bin_dPO_vect = [];
mindPO_vect = [];
for IT = 1:length(unique(comp_type(~isnan(comp_type))))
    if isempty(dPOsig_perint{IT})
        dPOsig_vect = cat(1,dPOsig_vect,nan);
        IT_bin_dPO_vect = cat(1,IT_bin_dPO_vect,IT);
    else
        dPOsig_vect = cat(1,dPOsig_vect,dPOsig_perint{IT});
        IT_bin_dPO_vect = cat(1,IT_bin_dPO_vect,ones(length(dPOsig_perint{IT}),1)*IT);
    end
    
    if isempty(dPOmin_perint{IT})
        mindPO_vect = cat(1,mindPO_vect,nan);
    else
        mindPO_vect = cat(1,mindPO_vect,dPOmin_perint{IT});
    end
end

% calculate fractions
numb_sig = cellfun(@(x) sum(~isnan(x)),dPOsig_perint,'UniformOutput', false); numb_sig = cat(2,numb_sig{:});
numb_tuned = cellfun(@(x) sum(~isnan(x)),dPOsig_perint_all,'UniformOutput', false); numb_tuned = cat(2,numb_tuned{:}); % only considering those that are conc tuned and vis resp (othewise they have no PO to compare)
numb_all = cellfun(@(x) length(x),dPOsig_perint_all,'UniformOutput', false); numb_all = cat(2,numb_all{:});
PO_bin_frac_sig_ofTuned = numb_sig./numb_tuned;

% calc number of cells in each orientation bin at each time point (output var)
for i = 1:size(PO_tuned,1)
    PO_to_PrefO_count(i,:) = histcounts(wrapTo180([PO_tuned(i,:)].*2)./2,PO_bins);
end
PO_to_PrefO_count(:,1) = PO_to_PrefO_count(:,1)+PO_to_PrefO_count(:,end);PO_to_PrefO_count(:,end)=[]; % combine the first and last bin (they should be equivalent)

% output for making all possible change plots
for IT_bin = 1:length(unique(comp_type(~isnan(comp_type))))
    output{IT_bin}(1,:)=start_PO_perint_all{IT_bin}; % this is normalized to permitted unless input specifies otherwise
    output{IT_bin}(2,:)=end_PO_perint_all{IT_bin}; 
    output{IT_bin}(3,:)=dPOsig_perint_all{IT_bin}; % change relative to true starting PO (not relative to permitted ori)
    output{IT_bin}(4,:)=issig_perint{IT_bin};
end

% calc bootstrapped CI for sig ratio
for IT_bin = 1:length(unique(comp_type(~isnan(comp_type))))
    theta_start_all = wrapTo180(([start_PO_perint_all{IT_bin}]')).*2;
    theta_end_all = wrapTo180(([end_PO_perint_all{IT_bin}]')).*2;
    theta_end_all(isnan(theta_start_all(:)))=[]; % get rid of all nans which are those with non vis or stat tuned timepoints
    theta_start_all(isnan(theta_start_all(:)))=[];
    theta_change_all = abs(theta_start_all)-abs(theta_end_all); % testing
    
    convergence_ = @(x) nanmedian(x);
    bin_edges = [-180:30:180];
    convergence_boot = bootstrp(1000,convergence_,theta_change_all./2);
    ratio_conTOdiv{IT_bin} = convergence_(theta_change_all./2);
    ratio_conTOdiv_CI{IT_bin} = [prctile(convergence_boot,5); prctile(convergence_boot,95)];
end


end