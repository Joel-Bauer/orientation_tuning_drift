function data = roi_correlations(data,channels,dataset)
clear all_roi_corrs all_randroi_corrs collect_rois
for mouse_n = 1:length(data)
    for cell_n = 1:length(data{mouse_n})
        days = 1:length(data{mouse_n}(cell_n).day);
        roi_corr_matrix = nan(days(end),days(end));
        for day_1 = days
            for day_2 = days
                if day_1<day_2
                    if isfield(data{mouse_n},'suite2P_template')
                        day_1_roi = double(data{mouse_n}(cell_n).suite2P_template{day_1});
                        day_2_roi = double(data{mouse_n}(cell_n).suite2P_template{day_2});
                    else
                        day_1_roi = data{mouse_n}(cell_n).roi_template{day_1};
                        day_2_roi = data{mouse_n}(cell_n).roi_template{day_2};
                    end
                    if ~isempty(day_1_roi) && ~isempty(day_2_roi)
                        roi_crop_size = size(day_1_roi,1);
                                               
                        % take average
                        day_1_roia = day_1_roi(1:roi_crop_size,1:roi_crop_size,1:channels);
                        day_2_roia = day_2_roi(1:roi_crop_size,1:roi_crop_size,1:channels);
                        try 
                            day_1_roib = day_1_roi(1:roi_crop_size,roi_crop_size+1:roi_crop_size*2,1:channels);
                        catch
                            day_1_roib = nan(size(day_1_roia));
                        end
                        try
                            day_2_roib = day_2_roi(1:roi_crop_size,roi_crop_size+1:roi_crop_size*2,1:channels);
                        catch
                            day_2_roib = nan(size(day_2_roia));
                            
                        end
                        day_1_roi = nanmean(cat(4,day_1_roia,day_1_roib),4);
                        day_2_roi = nanmean(cat(4,day_2_roia,day_2_roib),4);
                        
                        % remove black pixels (roi outlines)
                        day_1_roi(all(day_1_roi==0,3)) = nan;
                        day_2_roi(all(day_1_roi==0,3)) = nan;
                        
                        % rescale channels (not really needed for cor)
                        for ch = 1:channels
                            day_1_roi(:,:,ch) = rescale(day_1_roi(:,:,ch));
                            day_2_roi(:,:,ch) = rescale(day_2_roi(:,:,ch));
                        end                     
                        roi_corr_matrix(day_1,day_2) = corr(day_1_roi(:),day_2_roi(:));
                    end
                elseif day_1 == day_2 
                    roi_corr_matrix(day_1,day_2) = 0.5;
                end
                
            end
        end
        roi_corr_matrix = nansum(cat(3,roi_corr_matrix,triu(roi_corr_matrix)'),3);
        data{mouse_n}(cell_n).ROI_correlations = roi_corr_matrix;
    end
end
end
%%