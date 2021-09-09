function plot_multivar(sigma,accept,xStep,x_keep,x0,numIt,p,...
                       paramsVaried,axisLabels,maxNumBins,L_type,path,f_max,t_max)
                   
% get useful info
numParams = length(paramsVaried);
numPanelSide = numParams + 1;

% make gridded plots of all independent parameters
for i = 1:numParams
    for j = i:numParams
        
        if i ~= j
            
            % get number of bins
            numBins = length(unique(x_keep(paramsVaried(i),:)));
            if numBins > maxNumBins
                numBins = maxNumBins;
            end
            
            % get point with max density for current variable pair- use log
            % version for t0 (parameter 4)
            if paramsVaried(j) == 4
                xFit = getFitLog(x_keep,[paramsVaried(j) paramsVaried(i)],numBins,x0,'x');
            else
                xFit = getFit(x_keep,[paramsVaried(j) paramsVaried(i)],numBins,x0);
            end
            
            % get correct subplot indices
            figInd = sub2ind([numPanelSide,numPanelSide],j,i);
            subplot(numPanelSide,numPanelSide,figInd);           
            
            % make dscatter density plot of results- use log
            % version for t0 (parameter 4)
            dscatter(x_keep(paramsVaried(j),:)',x_keep(paramsVaried(i),:)','BINS',[numBins,numBins]);
            if paramsVaried(j) == 4
                dscatter(x_keep(paramsVaried(j),:)',x_keep(paramsVaried(i),:)','BINS',[numBins,numBins],'LOGX',true);
                xticks([1e-2 0.05 1e-1 0.5 1e0 5 1e1 50 1e2]);
                xticklabels({'0.01','0.05','0.1','0.5','1','5','10','50','100'})
            end
            
            % set axes positions
            ax = gca;
            ax.YAxisLocation = "right";
            ax.XAxisLocation = "top";
            
            % plot dashed red lines intersecting at max density point
            xline(xFit(paramsVaried(j)),"r--");
            yline(xFit(paramsVaried(i)),"r--");
            
            % set axes limits based on range of values 
            xlim([min(x_keep(paramsVaried(j),:)),max(x_keep(paramsVaried(j),:))]);
            ylim([min(x_keep(paramsVaried(i),:)),max(x_keep(paramsVaried(i),:))]);
            
            % force full outlines
            box on;
           
            % remove x-axis labels if current subplot does not fall along
            % edge of plot grid
            if figInd > length(paramsVaried)
                xticklabels(gca,{})
            else
                xlabel(axisLabels(paramsVaried(j)),'Color',[0 0.4470 0.7410])
            end
            
            % only plots on the right edge will have y-axis labels (and
            % these will all be M_0 plots, so remove y labels here) 
            yticklabels(gca,{})
            
        end
        
    end

    % get correct subplot indices for histogram
    histInd = sub2ind([numPanelSide,numPanelSide],i,i);
    subplot(numPanelSide,numPanelSide,histInd);
   
    % make histogram, using log handling if t0
    if paramsVaried(i) == 4
        [~,edges] = histcounts(log10(x_keep(paramsVaried(i),:)),numBins);
        histogram(x_keep(paramsVaried(i),:),10.^edges,'FaceColor',[0 0.4470 0.7410],...
                 'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',1);
        set(gca,'xscale','log');
        xticks([1e-2 0.05 1e-1 0.5 1e0 5 1e1 50 1e2]);
        xticklabels({'0.01','0.05','0.1','0.5','1','5','10','50','100'})
    else
        histogram(x_keep(paramsVaried(i),:),numBins,'FaceColor',[0 0.4470 0.7410],...
                 'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',1);
    end
    
    % set axes limits based on range of values and add labels
    xlim([min(x_keep(paramsVaried(i),:)),max(x_keep(paramsVaried(i),:))]);
    xlabel(axisLabels(paramsVaried(i)),'Color',[0 0.4470 0.7410])
end

% report parameters and settings for MCMC
subplot(numPanelSide,numPanelSide,numPanelSide*numPanelSide-(numPanelSide-1))
yticklabels(gca,{})
xticklabels(gca,{})
set(gca, 'visible', 'off')
text(0,1,string("MCMC parameters" + newline + "----------------------------" + newline + ...                            
                       "h_i step: " + xStep(1) + " m    h_w step: " + xStep(2) + " m" + newline + ...
                       "X_{stat} step: " + xStep(3) + " m    t_0 step: log10(" + round(10^xStep(4)) + ") s" + newline + ...
                       "Number of iterations: " + numIt + newline + "Sigma: " + sigma + newline + ...
                       "Liklihood function: " + L_type + newline + "Accepted " + ...
                       round(100*sum(accept)/length(accept)) + "% of proposals" + newline + newline + ...
                       "Model parameters" + newline + "----------------------------" + newline + ... 
                       "Sampling Frequency: " + 2*f_max + " Hz" + newline + "Duration: " + t_max + " s"))                   
            
% set figure size and title
set(gcf,'Position',[10 10 1200 1000])
sgtitle("Result of MCMC inversion after " + numIt + " iterations" + newline)

saveas(gcf,path + "run" + p + "_multivar_dscatter.png")
close(gcf)

end