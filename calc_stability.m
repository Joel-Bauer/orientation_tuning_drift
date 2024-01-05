function [xvals_ofmeanfit_dT,yvals_ofmeanfit_dT,...
    xvals_ofmeanfit_dS,yvals_ofmeanfit_dS,...
    xvals_ofmeanfit_dT_allmice,yvals_ofmeanfit_dT_allmice,...
    xvals_ofmeanfit_dS_allmice,yvals_ofmeanfit_dS_allmice,...
    F_equation_dT] = calc_stability(data, metric, showplots, fig_numb)

if nargin<4 | isempty(fig_numb)
   fig_numb = nan; 
end

if showplots
    plot_peranimal_dTdS = 0;
    plot_overall = 1;
else
    plot_peranimal_dTdS = 0;
    plot_overall = 0;
end

for mouse_n = 1:length(data)
    clear temp_all
    for celln = 1:length(data{mouse_n})
        if any(strcmp(metric,{'PO_corr','r2er'}))
            temp_all(:,:,celln) = data{mouse_n}(celln).PO;
        else
            temp_all(:,:,celln) = eval(['data{mouse_n}(celln).' metric]);
        end
    end
    
    if strcmp(metric,'POdif_all') || strcmp(metric,'POdif_sig')  % if using the PO measure take abolute values
        temp_all = abs(temp_all);
    elseif strcmp(metric,'PO_corr')
        temp_all = squeeze(temp_all);
    end
        
    if contains(metric,'cor') || strcmp(metric,'PO_corr')% if metric is a correlation measure
        FT_mean = 1; % use fisher transform for mean calculations
    else 
        FT_mean = 0;
    end
    
    clear temp_mean
    if strcmp(metric,'PO_corr') % needs to be calculated first
        for tp_a=1:size(temp_all,1)
            for tp_b=1:size(temp_all,1)
                temp_mean(tp_a,tp_b) = circ_corrcc_withUniCorrection(wrapToPi(deg2rad(temp_all(tp_a,:)*2)),wrapToPi(deg2rad(temp_all(tp_b,:)*2)),1);
            end
        end
    elseif strcmp(metric,'r2er') % needs to be calculated first
        for tp_a=1:size(temp_all,1)
            for tp_b=1:size(temp_all,1)
                if tp_a==tp_b
                    temp_mean(tp_a,tp_b)=nan;
                else
                    temp1 = arrayfun(@(x) permute(x.TC_reps{tp_a},[3,2,1]),data{mouse_n}(:),'UniformOutput', false);
                    temp2 = arrayfun(@(x) permute(x.TC_reps{tp_b},[3,2,1]),data{mouse_n}(:),'UniformOutput', false);
                    temp1 = cat(1,temp1{:});
                    temp2 = cat(1,temp2{:});
                    temp_mean(tp_a,tp_b) = nanmedian(r2er_n2n_JB(temp1,temp2));
                end
            end
        end
     elseif contains(metric,'cor')
        for tp_a=1:size(temp_all,1)
            for tp_b=1:size(temp_all,2)
                temp_mean(tp_a,tp_b) = cocoe_mean(squeeze(temp_all(tp_a,tp_b,:)));
            end
        end
    else
        for tp_a=1:size(temp_all,1)
            for tp_b=1:size(temp_all,2)
                temp_mean(tp_a,tp_b) = nanmedian(squeeze(temp_all(tp_a,tp_b,:)));
            end
        end
    end
    
    % comparison type matrices
    [expanded_measure{mouse_n}, comparison_type_matrix_time{mouse_n}, comparison_type_matrix_session{mouse_n}] = ...
        scale_time_course_matricies(temp_mean, data{mouse_n}(1).day', [], [], [], 4, [], FT_mean);
    [dTdS_mat] = display_feature_dTvsdS_mean(expanded_measure(mouse_n),comparison_type_matrix_time(mouse_n),comparison_type_matrix_session(mouse_n),0,0,0);
    
    if plot_peranimal_dTdS
        if mouse_n==1
            figure
        end
        subplot(3,length(data),mouse_n)
        display_feature_TvsT(expanded_measure(mouse_n),0,1);
        colormap('winter')
        if contains(metric,'cor'); caxis([0 1]); end
        if isfield(data{mouse_n}(1),'model')
            title(data{mouse_n}(1).model,'Interpreter','none')
        elseif isfield(data{mouse_n}(1),'mouse')
            title(data{mouse_n}(1).mouse,'Interpreter','none')
        end
        
        subplot(3,length(data),length(data)+mouse_n)
        display_feature_dTvsdS_mean(expanded_measure(mouse_n),comparison_type_matrix_time(mouse_n),comparison_type_matrix_session(mouse_n),0,1,1);
        colormap('winter')
        if contains(metric,'cor'); caxis([0 1]); end
        title([''])
        
        subplot(3,length(data),length(data)*2+mouse_n)
        
        clear dT_mat
        if contains(metric,'cor')
            for dT_n = 1:size(dTdS_mat{1},2)
                dT_mat(dT_n) = cocoe_mean(dTdS_mat{1}(:,dT_n));
            end
        else
            for dT_n = 1:size(dTdS_mat{1},2)
                dT_mat(dT_n) = nanmedian(dTdS_mat{1}(:,dT_n));
            end
        end
        
        plot([1:length(dT_mat)]-1,dT_mat); 
        box off;
        xlabel('Delta Time'); ylabel(metric,'Interpreter','none')
        set(gca,'TickDir','out')
        pbaspect([1 1 1])
        if contains(metric,'cor'); ylim([0 1]); end
        title([' ']);
    end
end

if plot_peranimal_dTdS
    a = suptitle(metric);
    a.Interpreter = 'none';
end

% overall plots
if plot_overall
    if isnan(fig_numb)
        [~, ~,xvals_ofmeanfit_dT,yvals_ofmeanfit_dT,xvals_ofmeanfit_dT_allmice,yvals_ofmeanfit_dT_allmice,...
            F_equation_dT] = ...
            display_feature_decay_curve2(...
            comparison_type_matrix_time,expanded_measure,[],[],1);
        [~, ~,xvals_ofmeanfit_dS,yvals_ofmeanfit_dS,xvals_ofmeanfit_dS_allmice,yvals_ofmeanfit_dS_allmice] = ...
            display_feature_decay_curve2(...
            comparison_type_matrix_session,expanded_measure,[],[],0);
    else
        [~, ~,xvals_ofmeanfit_dT,yvals_ofmeanfit_dT,xvals_ofmeanfit_dT_allmice,yvals_ofmeanfit_dT_allmice,...
            F_equation_dT] = ...
            display_feature_decay_curve2(...
            comparison_type_matrix_time,expanded_measure,[],fig_numb,1);
        [~, ~,xvals_ofmeanfit_dS,yvals_ofmeanfit_dS,xvals_ofmeanfit_dS_allmice,yvals_ofmeanfit_dS_allmice] = ...
            display_feature_decay_curve2(...
            comparison_type_matrix_session,expanded_measure,[],fig_numb,0);
    end
    
    xlabel('dT');ylabel(metric,'Interpreter','none');
else
        [~, ~,xvals_ofmeanfit_dT,yvals_ofmeanfit_dT,xvals_ofmeanfit_dT_allmice,yvals_ofmeanfit_dT_allmice,...
            F_equation_dT] = display_feature_decay_curve2(...
            comparison_type_matrix_time,expanded_measure,[],[],0);
        [~, ~,xvals_ofmeanfit_dS,yvals_ofmeanfit_dS,xvals_ofmeanfit_dS_allmice,yvals_ofmeanfit_dS_allmice] = display_feature_decay_curve2(...
            comparison_type_matrix_session,expanded_measure,[],[],0);        
end

end
%%

