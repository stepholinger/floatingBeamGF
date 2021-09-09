function mcmc_object = run_GF_mcmc(param_file,out_path)

% this code runs a version of P. Segall's mcmc code adapted for fitting the beam deflection
% Green's function with observed seismic data

% get parameter function
param_func = str2func(param_file);

% run setup to initialize mcmc and model parameters and objects
try
    mcmc_object = centroid_mcmc_setup(param_func);
catch
    mcmc_object = mseed_mcmc_setup(param_func);
end

% if h_i and h_w parameters are fixed, calculate matrix of green's function
% solutions once before mcmc
if sum([mcmc_object.steps.h_i,mcmc_object.steps.h_w]) == 0
    G = semiAnalyticGreenFunction(mcmc_object.model);
    
    % take spatial derivative if moment source
    if mcmc_object.source == "moment"
        [~,G] = gradient(G,mcmc_object.model.dx);
    end
    
    mcmc_object.G_matrix = G;
end

% generate intial Green's function
[G_0,xcorr_coef] = GF_func_mcmc(mcmc_object);
mcmc_object.G_0 = G_0;

% calculate initial liklihood
L_0 = liklihood(G_0,mcmc_object,xcorr_coef);
mcmc_object.L_0 = L_0;

% run mcmc
tic;
mcmc_object = mcmc(mcmc_object);
runtime = toc;

if contains("errors",fieldnames(mcmc_object))

      % give output upon completion
    fprintf("\n--------------------------------\n")
    fprintf("\nRan into error on iteration " + mcmc_object.iteration + "!\n")
    
    % save results if output filename provided
    if exist("out_path","var")
        save(out_path + ".mat","mcmc_object")
    end
    
else
    
    % give output upon completion
    fprintf("\n--------------------------------\n")
    fprintf("\nFinished MCMC in " + runtime + " seconds\n")
    fprintf("\nAccepted " + round((sum(mcmc_object.results.accept)/mcmc_object.num_iterations)*100) + " %% of proposals\n\n");

    % save results if output filename provided
    if exist("out_path","var")
        save(out_path + ".mat","mcmc_object")
    end

    % make plot
    plot_multivar(mcmc_object,out_path)
end

end