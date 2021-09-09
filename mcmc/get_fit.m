function x_fit = get_fit(mcmc_object,params_varied,num_bins)
    
    % set xFit
    x_fit = [mcmc_object.x0.h_i,mcmc_object.x0.h_w,mcmc_object.x0.stat_dist,mcmc_object.x0.t0,mcmc_object.x0.h_c_initial,mcmc_object.x0.h_c_final];

    % get matrix back from dscatter
    figure(2) 
    [~,~,col] = dscatter(mcmc_object.results.x_keep(params_varied(1),:)',mcmc_object.results.x_keep(params_varied(2),:)','BINS',[num_bins num_bins]);    
    close(gcf);
    
    % get index of highest density point
    [~,fit_idx] = max(col);
    
    % get values at that index
    x_fit(params_varied(1)) = mcmc_object.results.x_keep(params_varied(1),fit_idx);
    x_fit(params_varied(2)) = mcmc_object.results.x_keep(params_varied(2),fit_idx);
    
end