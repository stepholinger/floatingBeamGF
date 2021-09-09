function plot_multivar(mcmc_object,out_path)
            
% set internal parameters
max_num_bins = 100;

% get useful info
x_step = [mcmc_object.steps.h_i;mcmc_object.steps.h_w;mcmc_object.steps.stat_dist;...
          mcmc_object.steps.t0;mcmc_object.steps.h_c_initial;mcmc_object.steps.h_c_final];
x_keep = mcmc_object.results.x_keep;
axis_labels = ["Ice thickness (m)", "Water depth (m)", "X_{stat} (km)",...
    "t_0 (s)","Initial h_c (% ice thickness)","Final h_c (% ice thickness)"];
param_ind = [1,2,3,4,5,6];
params_varied = param_ind(x_step ~= 0);
num_params = length(params_varied);
num_panel_per_side = num_params + 1;

% convert X_stat to km and crevasse heights to percentage
x_keep(3,:) = x_keep(3,:)/1000;
x_keep(5:6,:) = 100*(0.5+x_keep(5:6,:));

% make gridded plots of all independent parameters
for i = 1:num_params
    
    % get number of bins
    num_bins = length(unique(x_keep(params_varied(i),:)));
    if num_bins > max_num_bins
        num_bins = max_num_bins;
    end
    
    for j = i:num_params
        
        if i ~= j
            
            % get point with max density for current variable pair
            x_fit = get_fit(mcmc_object,[params_varied(j) params_varied(i)],num_bins);
            
            % get correct subplot indices
            fig_ind = sub2ind([num_panel_per_side,num_panel_per_side],j,i);
            subplot(num_panel_per_side,num_panel_per_side,fig_ind);           
            
            % make dscatter density plot of results
            dscatter(x_keep(params_varied(j),:)',x_keep(params_varied(i),:)','BINS',[num_bins,num_bins]);
            
            % set axes positions
            ax = gca;
            ax.YAxisLocation = "right";
            ax.XAxisLocation = "top";
            
            % plot dashed red lines intersecting at max density point
            xline(x_fit(params_varied(j)),"r--");
            yline(x_fit(params_varied(i)),"r--");
            
            % set axes limits based on range of values 
            xlim([min(x_keep(params_varied(j),:)),max(x_keep(params_varied(j),:))]);
            ylim([min(x_keep(params_varied(i),:)),max(x_keep(params_varied(i),:))]);
            
            % force full outlines
            box on;
           
            % remove x-axis labels if current subplot does not fall along
            % edge of plot grid
            if fig_ind > length(params_varied)
                xticklabels(gca,{})
            else
                xlabel(axis_labels(params_varied(j)),'Color',[0 0.4470 0.7410])
            end
            
            % only plots on the right edge will have y-axis labels (and
            % these will all be M_0 plots, so remove y labels here) 
            yticklabels(gca,{})
            
        end
        
    end

    % get correct subplot indices for histogram
    hist_ind = sub2ind([num_panel_per_side,num_panel_per_side],i,i);
    subplot(num_panel_per_side,num_panel_per_side,hist_ind);
   
    % make histogram, using log handling if t0
    histogram(x_keep(params_varied(i),:),num_bins,'FaceColor',[0 0.4470 0.7410],...
                 'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',1);
    
    % set axes limits based on range of values and add labels
    xlim([min(x_keep(params_varied(i),:)),max(x_keep(params_varied(i),:))]);
    xlabel(axis_labels(params_varied(i)),'Color',[0 0.4470 0.7410])
end

% report parameters and settings for MCMC
subplot(num_panel_per_side,num_panel_per_side,num_panel_per_side*num_panel_per_side-(num_panel_per_side-1))
yticklabels(gca,{})
xticklabels(gca,{})
set(gca, 'visible', 'off')

if sum(contains(["basal","surface","hydrostatic"],mcmc_object.mode)) > 0 
    text(0,1,string("MCMC parameters" + newline + "----------------------------" + newline + ...                            
                    "h_i step: " + x_step(1) + " m    h_w step: " + x_step(2) + " m" + newline + ...
                    "X_{stat} step: " + x_step(3) + " m    t_0 step: " + x_step(4) + " s" + newline + ...
                    "h_c step: " + x_step(5) + newline + "Number of iterations: " + mcmc_object.num_iterations + newline + ... 
                    "Sigma: " + mcmc_object.sigma + newline + "Liklihood function: " + mcmc_object.L_type + newline + ...
                    "Accepted " + round(100*sum(mcmc_object.results.accept)/length(mcmc_object.results.accept)) + ...
                    "% of proposals" + newline + newline + "Model parameters" + newline + ...
                    "----------------------------" + newline + "Sampling Frequency: " + 1/mcmc_object.model.dt + " Hz" + newline + ...
                    "Duration: " + mcmc_object.model.t_max + " s" + newline + "Crevasse type: " + mcmc_object.mode))           
else
    text(0,1,string("MCMC parameters" + newline + "----------------------------" + newline + ...                            
                    "h_i step: " + x_step(1) + " m    h_w step: " + x_step(2) + " m" + newline + ...
                    "X_{stat} step: " + x_step(3) + " m    t_0 step: " + x_step(4) + " s" + newline + ...
                    "h_c_initial step: " + x_step(5) + "    h_c_final step: " + x_step(6) + newline + ...
                    "Number of iterations: " + mcmc_object.num_iterations + newline + "Sigma: " + mcmc_object.sigma + newline + ...
                    "Liklihood function: " + mcmc_object.L_type + newline + "Accepted " + ...
                     round(100*sum(mcmc_object.results.accept)/length(mcmc_object.results.accept)) + "% of proposals" + newline + newline + ...
                     "Model parameters" + newline + "----------------------------" + newline + ... 
                     "Sampling Frequency: " + 1/mcmc_object.model.dt + " Hz" + newline + "Duration: " + mcmc_object.model.t_max + " s" + newline + ...
                     "Crevasse type: " + mcmc_object.mode))         
end
% set figure size and title
set(gcf,'Position',[10 10 1200 1000])
sgtitle("Result of MCMC inversion after " + mcmc_object.num_iterations + " iterations" + newline)

saveas(gcf,out_path + ".png")
close(gcf)

end