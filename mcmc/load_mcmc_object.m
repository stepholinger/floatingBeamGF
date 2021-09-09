function mcmc_object = load_mcmc_object(model,data,stat_dist,t0,source,mode,h_c_initial,h_c_final)

% set model parameters
mcmc_object.model = model;

% set data to fit with mcmc
mcmc_object.data = data;

% set inital parameters for mcmc
mcmc_object.x0.h_i = model.h_i;
mcmc_object.x0.h_w = model.h_w;
mcmc_object.x0.stat_dist = stat_dist;
mcmc_object.x0.t0 = t0;

% set source and crevasse mode parameters- crevasse heights will do nothing 
% if using erf or gaussian stf
mcmc_object.source = source;
mcmc_object.mode = mode;
mcmc_object.x0.h_c_initial = h_c_initial;
mcmc_object.x0.h_c_final = h_c_final;

end
