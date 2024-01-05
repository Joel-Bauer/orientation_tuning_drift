function [expanded_matrix, comparison_type_matrix_time, comparison_type_matrix_session, comparison_type_matrix_decay, comparison_type_matrix_time_from_onset_firstday, comparison_type_matrix_time_from_onset_secondday] = ...
    scale_time_course_matricies(input_matrix, all_imaging_days, included_days, latest_date, reptype, reptype_numb, Tau_decay, corr_matrix)

%%
% INPUTS:
%   input_matrix, an unexpanded matrix of comparison values (eg. spine
%       turnover, PSTH correlations...)
%   all_imaging_days, a vector  which idicates the number of days since 
%       the start of the imaging +1 (eg. 1 2 3 5 for mon, tus, wed, fri).
%       can be equal to or longer than input_matrix.
%   included_days, a vector of the same length as input_matrix which idicates
%       the number of days since the start of the imaging +1. if empty then
%       included_days = all_imaging_days
%   latest_date, the number of days from the first experiment to the last
%       +1. can be left empty (as  []) in which case  all_imaging_days(end) is
%       used.
%   reptype, a vector the same length as included_days indicating the
%       reptype of the experiment (within day comparison type). can be left as
%       [], in which case all comparisons will have the reptype 1.
%   reptype_numb, number of different with day rep types. can be left empty
%       in which case default is 1.
%   Tau_decay, decay constant used for the feature decay model (nonlinear
%       time depedent decay of measured feature)

%problem: exlcuded comparisons are not shown in the models but they are
%shown in the data. check

showplot = 0;
session_based_model = 'imaging day counter'
% session_based_model = 'delayed counter'
% fill in empty variables
if isempty(included_days)
    included_days = all_imaging_days;
end

if isempty(latest_date)
    max_numb_days = max(all_imaging_days); %What is the latest time points (in (days from session 1) + 1 )
else 
    max_numb_days = latest_date;
end

if isempty(reptype)
    reptype = ones(size(included_days,1),size(included_days,2));
    reptype(find(diff(all_imaging_days)==0)+1)=2;
end

% if strcmp(session_based_model,'imaging day counter')
%     reptype_numb = 1;
% elseif isempty(reptype_numb)
%     reptype_numb = 1;
% end

if isempty(reptype_numb)
    reptype_numb = 1;
end

if isempty(Tau_decay)
    Tau_decay = 2;
end

if ~exist('corr_matrix')
    corr_matrix = 0;
end
    
if reptype_numb == 1
    session_feature_vector_rep_type_all = ones(1,length(all_imaging_days));
    session_feature_vector_rep_type_included = ones(1,length(included_days));
else
    session_feature_vector_rep_type_included = reptype;
    session_feature_vector_rep_type_all = reptype ; %ones(1,length(all_imaging_days));
end


%% Make empty matrix of corrected size
max_numb_rep_types = reptype_numb;
empty_master_matrix = NaN(max_numb_days*max_numb_rep_types,max_numb_days*max_numb_rep_types);

%% Determine direction of change (so that the models can have ether increasing or decreasing amp)
% DONT USE DOES NOT WORK
%use the first 3 rows of the matrix to determine if they are increasing or
%decreasing on average.
% if 0<mean([mean(diff(input_matrix(1,2:end))),mean(diff(input_matrix(2,3:end))),mean(diff(input_matrix(3,4:end)))])
%     %increasing
%     corr_matrix = 0;
% elseif 0<mean([mean(diff(input_matrix(1,2:end))),mean(diff(input_matrix(2,3:end))),mean(diff(input_matrix(3,4:end)))])
%     %decreasing
%     corr_matrix = 1;
% end

%% Make time interval comparison type matrix 1 (to indicate which datapoints compare time points of different lengths of time)
axis_vector_rep_type = repmat(1:max_numb_rep_types,[1,max_numb_days]);
axis_vector_days = repmat(1:max_numb_days,[max_numb_rep_types,1]);
axis_vector_days = axis_vector_days(:);
comparison_type_matrix_time = zeros(max_numb_days*max_numb_rep_types);

for ii = 1:size(empty_master_matrix,1) % go row
    for i = 1:size(empty_master_matrix,2) % by row
        if axis_vector_days(i) == axis_vector_days(ii) %if it is a same day correlation
            % current within day comparison type 
            within_day_type = [axis_vector_rep_type(i) axis_vector_rep_type(ii)];
            within_day_comp_type1 = [1 2; 3 4; 2 1; 4 3]; % comparisons to repetitions with no eye shutter changes
            within_day_comp_type2 = [1 3; 1 4; 2 3; 2 4; 3 1; 4 1; 3 2; 4 2]; %comparisons to repetitions with eye shutter changes
            
            % within day comparisons are indicated by negative numbers 
            if ismember(within_day_type,within_day_comp_type1,'rows')
                comparison_type_matrix_time(i,ii) = -1; % comparisons types with no shutter change
            elseif ismember(within_day_type,within_day_comp_type2,'rows')
                comparison_type_matrix_time(i,ii) = -2; % comparisons types with a shutter change
            end
        else
            abs(axis_vector_days(i)-axis_vector_days(ii));
            if abs(axis_vector_days(i)-axis_vector_days(ii)) > 0 && abs(axis_vector_days(i)-axis_vector_days(ii)) 
                comparison_type_matrix_time(i,ii) = abs(axis_vector_days(i)-axis_vector_days(ii));
            end
        end
    end
end

%% Make session interval comparison type matrix 2 (to indicate which data points compare time points with different number of imaging sessions between them)

if strcmp(session_based_model,'include all')
    % number of imaging sessions (including all sessions)
    temp = [NaN,0:length(all_imaging_days)-2];
    comparison_type_matrix_session_nonexpanded = NaN(length(all_imaging_days),length(all_imaging_days));
    for i = 1:length(all_imaging_days)
        for ii = 1:length(all_imaging_days)
            comparison_type_matrix_session_nonexpanded(i,ii) = temp(abs(diff([i,ii]))+1);
        end
    end
    
elseif strcmp(session_based_model,'delayed counter')
    % delayed effect of imaging session (counter of session only increases
    % 'overnight'). if same day dS = 0. if imaged with no additional tp
    % inbetween then dS=1.
    comparison_type_matrix_session_nonexpanded = NaN(length(all_imaging_days),length(all_imaging_days));
    same_day_sessions = find(diff(all_imaging_days)==0)+1;
    delayed_effect_counter = 1:length(all_imaging_days);
    for j = 1:length(same_day_sessions)
        delayed_effect_counter(same_day_sessions(j)) = delayed_effect_counter(find(all_imaging_days == all_imaging_days(same_day_sessions(j)),1));
    end
    for i = 1:length(all_imaging_days)
        for ii = 1:length(all_imaging_days)
            if all_imaging_days(i) == all_imaging_days(ii) && i ~= ii
                comparison_type_matrix_session_nonexpanded(i,ii) =0;
            elseif i == ii
                continue
            else
                comparison_type_matrix_session_nonexpanded(i,ii) = abs(diff([delayed_effect_counter(i),delayed_effect_counter(ii)]));
            end
        end
    end
    
elseif strcmp(session_based_model,'imaging day counter')
    % number of imaging days (the number of days between the time points on which imaging was performed, regardless of how many reps per day)
    comparison_type_matrix_session_nonexpanded = NaN(length(all_imaging_days),length(all_imaging_days));
    n = 1;
    temp = 1;
    for i = 2:length(all_imaging_days)
        if all_imaging_days(i-1) < all_imaging_days(i)
            n = n+1;
        end
        temp(i) = n;
    end
    for i = 1:length(all_imaging_days)
        for ii = 1:length(all_imaging_days)
%             comparison_type_matrix_session_nonexpanded(i,ii) = abs(diff([temp(i),temp(ii)]))-1; % This makes all within day sessions nans and does not count the first imaging day as an interleaved session
            comparison_type_matrix_session_nonexpanded(i,ii) = abs(diff([temp(i),temp(ii)])); % count days of imaging between tps, within day sessions are nans and 1 day apart > 0;
            
%             % delayed counter, counting all imaging sessions not just
%             % imaging days
%             if temp(i) == temp(ii)
%                 comparison_type_matrix_session_nonexpanded(i,ii) = 0;
%             else
%                 comparison_type_matrix_session_nonexpanded(i,ii) = abs(diff([find(temp==temp(i),1),find(temp==temp(ii),1)]));
%             end
            
            if i == ii
                comparison_type_matrix_session_nonexpanded(i,ii) = NaN;
            end
        end
    end
    comparison_type_matrix_session_nonexpanded(comparison_type_matrix_session_nonexpanded<0) = NaN;
    
end

% expantion
comparison_type_matrix_2_expanded = empty_master_matrix;
for i = 1:size(comparison_type_matrix_2_expanded,1) % go row
    for ii = 1:size(comparison_type_matrix_2_expanded,2) %by row to fill in expanded matrix
        % Find coresponding pixels between expanded and nonexpanded
        % cell corr matrices
        current_day_i = ceil(i/max_numb_rep_types);
        current_rep_type_i = (i/max_numb_rep_types-current_day_i+1)*max_numb_rep_types;
        current_day_ii = ceil(ii/max_numb_rep_types);
        current_rep_type_ii = (ii/max_numb_rep_types-current_day_ii+1)*max_numb_rep_types;
        index_i = []; index_ii = [];
        same_days_i=all_imaging_days==current_day_i;
        same_reptype_i=session_feature_vector_rep_type_all==current_rep_type_i;
        index_i = find(same_days_i.*same_reptype_i);
        same_days_ii=all_imaging_days==current_day_ii;
        same_reptype_ii=session_feature_vector_rep_type_all==current_rep_type_ii ;
        index_ii = find(same_days_ii.*same_reptype_ii);
        
        % Coppy cor value from none expanded to expanded cell corr
        % matrix
        if ~isempty(index_i) && ~isempty(index_ii)
            comparison_type_matrix_2_expanded(i,ii) = comparison_type_matrix_session_nonexpanded(index_i,index_ii);
        end
    end
end
comparison_type_matrix_session = comparison_type_matrix_2_expanded;

%% Make time interval comparison type matrix 3 (to indicate which data points compare time points with different number of imaging sessions between them)
 tau1 = Tau_decay;
 t = 1:length(all_imaging_days);
 yinitial = exp(-t/(tau1+1));
 comparison_type_matrix_decay = NaN(length(all_imaging_days),length(all_imaging_days));
for i = 1:length(all_imaging_days)
    for ii = 1:length(all_imaging_days)
        if i ~= ii 
            comparison_type_matrix_decay(i,ii) = abs(exp(-i/(tau1+1))-exp(-ii/(tau1+1)));
        end
    end
end

comparison_type_matrix_3_expanded = empty_master_matrix;
for i = 1:size(comparison_type_matrix_3_expanded,1) % go row
    for ii = 1:size(comparison_type_matrix_3_expanded,2) %by row to fill in expanded matrix
        % Find coresponding pixels between expanded and nonexpanded
        % cell corr matrices
        current_day_i = ceil(i/max_numb_rep_types);
        current_rep_type_i = (i/max_numb_rep_types-current_day_i+1)*max_numb_rep_types;
        current_day_ii = ceil(ii/max_numb_rep_types);
        current_rep_type_ii = (ii/max_numb_rep_types-current_day_ii+1)*max_numb_rep_types;
        index_i = []; index_ii = [];
        same_days_i=all_imaging_days==current_day_i;
        same_reptype_i=session_feature_vector_rep_type_all==current_rep_type_i;
        index_i = find(same_days_i.*same_reptype_i);
        same_days_ii=all_imaging_days==current_day_ii;
        same_reptype_ii=session_feature_vector_rep_type_all==current_rep_type_ii ;
        index_ii = find(same_days_ii.*same_reptype_ii);
        
        % Coppy cor value from none expanded to expanded cell corr
        % matrix
        if ~isempty(index_i) && ~isempty(index_ii)
            comparison_type_matrix_3_expanded(i,ii) = comparison_type_matrix_decay(index_i,index_ii);
        end
    end
end
comparison_type_matrix_decay = comparison_type_matrix_3_expanded;
%% Make time from onset (of first session) comparison type matrix 4
% comparison_type_matrix_time_from_onset_firstday = NaN(length(all_imaging_days),length(all_imaging_days));
% for i = 1:length(all_imaging_days)
%     for ii = 1:length(all_imaging_days)
%         comparison_type_matrix_time_from_onset_firstday(i,ii) = min(i,ii);
%     end
% end

comparison_type_matrix_4_expanded = empty_master_matrix;
for i = 1:size(comparison_type_matrix_4_expanded,1) % go row
    for ii = 1:size(comparison_type_matrix_decay,2) %by row to fill in expanded matrix
        if ~isnan(comparison_type_matrix_session(i,ii))
            comparison_type_matrix_4_expanded(i,ii) = floor(min([i,ii])/max_numb_rep_types);
        end
    end
end
comparison_type_matrix_time_from_onset_firstday = comparison_type_matrix_4_expanded;

%% Make time from onset (of second session) comparison type matrix 5
% comparison_type_matrix_time_from_onset_secondday = NaN(length(all_imaging_days),length(all_imaging_days));
% for i = 1:length(all_imaging_days)
%     for ii = 1:length(all_imaging_days)
%         comparison_type_matrix_time_from_onset_secondday(i,ii) = max(i,ii);
%     end
% end

comparison_type_matrix_5_expanded = empty_master_matrix;
for i = 1:size(comparison_type_matrix_5_expanded,1) % go row
    for ii = 1:size(comparison_type_matrix_decay,2) %by row to fill in expanded matrix
        if ~isnan(comparison_type_matrix_session(i,ii))
            comparison_type_matrix_5_expanded(i,ii) = floor(max([i,ii])/max_numb_rep_types);
        end
    end
end
comparison_type_matrix_time_from_onset_secondday = comparison_type_matrix_5_expanded;

%%
comparison_type_matrix_time(isnan(comparison_type_matrix_time.*comparison_type_matrix_session)) = NaN;

if showplot
    figure; h = imagesc(comparison_type_matrix_time); hold on
    colorbar; set(h,'AlphaData',~isnan(comparison_type_matrix_time));
    
    figure; h = imagesc(comparison_type_matrix_session); hold on
    colorbar; set(h,'AlphaData',~isnan(comparison_type_matrix_session));
    
    figure; h = imagesc(comparison_type_matrix_decay); hold on
    colorbar; set(h,'AlphaData',~isnan(comparison_type_matrix_decay));
    
    figure; h = imagesc(comparison_type_matrix_time_from_onset_secondday); hold on
    colorbar; set(h,'AlphaData',~isnan(comparison_type_matrix_decay));
end

%% Make expanded correlation matrix
expanded_matrix = empty_master_matrix;
for i = 1:size(expanded_matrix,1) % go row
    for ii = 1:size(expanded_matrix,2) %by row to fill in expanded matrix
        % Find coresponding pixels between expanded and nonexpanded
        % cell corr matrices
        current_day_i = ceil(i/max_numb_rep_types);
        current_rep_type_i = (i/max_numb_rep_types-current_day_i+1)*max_numb_rep_types;
        current_day_ii = ceil(ii/max_numb_rep_types);
        current_rep_type_ii = (ii/max_numb_rep_types-current_day_ii+1)*max_numb_rep_types;
        index_i = []; index_ii = [];
        same_days_i=included_days==current_day_i;
        same_reptype_i=session_feature_vector_rep_type_included==current_rep_type_i;
        index_i = find(same_days_i.*same_reptype_i);
        same_days_ii=included_days==current_day_ii;
        same_reptype_ii=session_feature_vector_rep_type_included==current_rep_type_ii ;
        index_ii = find(same_days_ii.*same_reptype_ii);
        
        % Coppy cor value from none expanded to expanded cell corr
        % matrix
        if ~isempty(index_i) && ~isempty(index_ii);
            expanded_matrix(i,ii) = input_matrix(index_i,index_ii); %KEY VARIABLE
        end
        if i == ii 
            expanded_matrix(i,ii) = NaN;
        end
    end
end

expanded_matrix(isnan(comparison_type_matrix_session)) = nan; % exclude comparisons from the data matrix which are not included in the session comparison type used


if showplot
    figure; h = imagesc(expanded_matrix); color_test = cbrewer('seq','YlGnBu',100); colormap(color_test); colorbar; hold on
    caxis('auto'); set(h,'AlphaData',~isnan(expanded_matrix));
end


if strcmp(session_based_model,'imaging day counter')
    % under these conditions that matricies need to be compressed so that
    % all comparisons per day are averaged and the matrix resized
    
    comparison_type_matrix_time(comparison_type_matrix_time<0)=0;
    
    temp =repmat(1:max(all_imaging_days),reptype_numb,1);
    day_index = temp(:)';
    for i = 1:max(all_imaging_days) % go row
        for ii = 1:max(all_imaging_days) %by row to fill in expanded matrix
            i_idx = find(day_index==i);
            ii_idx = find(day_index==ii);
            
            all_values = expanded_matrix(i_idx,ii_idx);
            if max(expanded_matrix(:))<=1 & min(expanded_matrix(:))>=-1 & corr_matrix == 1
                M1(i,ii) = cocoe_mean(all_values);
            else
                M1(i,ii) = nanmean(all_values(:));
            end
            
            all_values = comparison_type_matrix_time(i_idx,ii_idx);
            M2(i,ii) = nanmean(all_values(:));
            
            all_values = comparison_type_matrix_session(i_idx,ii_idx);
            M3(i,ii) = nanmean(all_values(:));
            
            all_values = comparison_type_matrix_decay(i_idx,ii_idx);
            M4(i,ii) = nanmean(all_values(:));
            
            all_values = comparison_type_matrix_time_from_onset_firstday(i_idx,ii_idx);
            M5(i,ii) = nanmean(all_values(:));
            
            all_values = comparison_type_matrix_time_from_onset_secondday(i_idx,ii_idx);
            M6(i,ii) = nanmean(all_values(:));
        end
    end
    expanded_matrix = real(M1);
    comparison_type_matrix_time = real(M2);
    comparison_type_matrix_session = real(M3);
    comparison_type_matrix_decay = real(M4);
    comparison_type_matrix_time_from_onset_firstday = real(M5);
    comparison_type_matrix_time_from_onset_secondday = real(M6);
end




end

function [r_mean, numb_of_non_nan_elements] = cocoe_mean(r_vector) % correlation coefficient mean calculator, with graphic representation of transformation
r_vector(isnan(r_vector)) = []; % remove nan elements
numb_of_non_nan_elements = length(r_vector);

if numb_of_non_nan_elements == 0
    r_mean = nan;
elseif numb_of_non_nan_elements == 1
    r_mean = r_vector;
elseif numb_of_non_nan_elements > 1
    z_vector = arrayfun(@(r)atanh(r), r_vector);
    z_mean = mean(z_vector);
    r_mean = tanh(z_mean);
end

end