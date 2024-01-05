function [tau, mean_values_yandx,xvals_ofmeanfit,yvals_ofmeanfit,xvals_ofallfits,yvals_ofallfits,F_equation] = display_feature_decay_curve2(varx,vary,varc,fig,showplots)

%%
if nargin<5
    showplots = 1;
end

if showplots
    if ~isempty(fig)
        figure(fig); set(gcf,'color','w')
    else
        figure; set(gcf,'color','w')
    end
end

x_val_all = [];
y_val_all  =[];
for mouse = 1:size(varx,2)
    clear var1_current var2_current var3_current vary_current varc_current
    varx_current = varx{mouse};
    temp = triu(ones(size(varx_current,1),size(varx_current,1))); temp(temp==0)=nan;
    varx_current = varx_current.*temp;
    clear popvec_corr_at_dt day1_corr_at_dt
    for i = 1:size(varx_current,1)
        [idx_row,idx_col] = find(varx_current==i); %
        for ii = 1:length(idx_row)
            vary_current{i,ii} = vary{mouse}(idx_row(ii),idx_col(ii));
            if ~isempty(varc)
                varc_current{i,ii} = varc{mouse}(idx_row(ii),idx_col(ii));
            end
        end
    end
    
    clear y_val x_val c_val
    y_val = []; x_val = []; c_val = [];
    for i = 1:size(vary_current,1)
        y_val = cat(2,y_val,[vary_current{i,:}]);
        x_val = cat(2,x_val,[i*ones(1,size([vary_current{i,:}],2))]);
        if ~isempty(varc)
            c_val= cat(2,c_val,[varc_current{i,:}]);
        end
    end
    y_val_all = cat(2,y_val_all,y_val);
    x_val_all = cat(2,x_val_all,x_val);
    
    clear ymean
    for i = unique(x_val)
        ymean(i)=mean(y_val(find(x_val==i)));
    end
    
    mean_values_yandx{mouse}(1,:) = ymean;
    mean_values_yandx{mouse}(2,:) = unique(x_val);
    max_dt_mouse(mouse) = max(x_val);
    
    try
        f = @(a,x) a(1)*exp(-a(2)*x)+a(3);
        A{mouse} = fmincon(@(a) sum((y_val - f(a,x_val)).^2), [0.5 0.5 0.5],[],[],[],[],[-Inf -Inf 0],[Inf Inf Inf])';
    catch
        A{mouse} = [];
    end
    if showplots & size(varx,2)<10
        if size(varx,2)>1
            ax(mouse) = subplot(size(varx,2),2,mouse*2-1);
        end
        if ~isempty(varc)
            scatter(x_val,y_val,[],c_val,'filled'); set(gca,'TickDir','out'); hold on;
            colorbar
        else
            scatter(x_val,y_val,[],'k','filled','MarkerFaceAlpha',0.5); set(gca,'TickDir','out'); hold on
        end
        plot(unique(x_val),ymean,'k-','Linewidth',2)
        try
            plot(0:0.2:size(vary_current,1),f(A{mouse},0:0.2:size(vary_current,1)),'r','Linewidth',2);
        end
        box off; pbaspect([2 1 1])
        
        if min(y_val)>0 & max(y_val)<0
            ylim([0 1]);
        end
    end
    
end

if showplots & size(varx,2)<10
    if size(varx,2)>1
        linkaxes(ax)
        linkprop(ax,{'XTick'})
    end
    ylim([0 1]); pbaspect([2 1 1])
    xticks([0:5:max(xlim)])
    if ~isempty(varc)
        linkprop(ax,{'CLim'})
        colormap([0:0.01:1;zeros(1,101);1:-0.01:0]')
    end
end

%     if size(varx,2)<10
%         for mouse = 1:size(varx,2)
%             subplot(size(varx,2),1,mouse)
%             xlim([0 max(max_dt_mouse)+1])
%         end
%     end
try
    A_mean = mean([A{:}],2);
    A_all= fmincon(@(a) sum((y_val_all - f(a,x_val_all)).^2), [0.5 0.5 0.5],[],[],[],[],[-Inf -Inf 0],[Inf Inf Inf])';
    xvals_ofmeanfit = 1:0.2:max(max_dt_mouse);
    yvals_ofmeanfit = f(A_all,1:0.2:max(max_dt_mouse));
    
    for  mouse = 1:size(varx,2)
        xvals_ofallfits{mouse} = 1:0.2:max_dt_mouse(mouse);
        yvals_ofallfits{mouse} = f(A{mouse},1:0.2:max_dt_mouse(mouse));
    end
    
    if size(varx,2)>1
        if showplots
            subplot(size(varx,2),2,2:2:size(varx,2)*2)
            for  mouse = 1:size(varx,2)
                plot(1:0.2:max_dt_mouse(mouse),f(A{mouse},1:0.2:max_dt_mouse(mouse)),'Color',[0.5 0.5 0.5]); hold on
            end
            plot(xvals_ofmeanfit,yvals_ofmeanfit,'r','LineWidth',2);
            xlim([0 max(max_dt_mouse)+1 ]); set(gca,'TickDir','out');
            if min(y_val)>0 & max(y_val)<0
                ylim([0 1]);
            end
            ylim([0 1]);pbaspect([2 1 1]); box off
        end
    end
    tau = [A{:}];
    tau = tau(2,:);
    F_equation = [num2str(A_all(1),3) '*e^(' num2str(-1*A_all(2),3) '*x) + ' num2str(A_all(3),3)];
    
catch
    tau = [];
    yvals_ofmeanfit = [];
    xvals_ofmeanfit = [];
end
end