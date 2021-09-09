function mcmc_object = set_mcmc_bounds(mcmc_object,h_i_bounds,h_w_bounds,stat_dist_bounds,t0_bounds,h_c_initial_bounds,h_c_final_bounds)

% set bounds for mcmc paramters
mcmc_object.bounds.h_i = h_i_bounds;
mcmc_object.bounds.h_w = h_w_bounds;

% if station distance or stf duration are less than spatial / temporal
% resolution of model, bound at those values
if stat_dist_bounds(1) < mcmc_object.model.dx
    stat_dist_bounds(1) = mcmc_object.model.dx;
end
mcmc_object.bounds.stat_dist = stat_dist_bounds;
if t0_bounds(1) < mcmc_object.model.dt
    t0_bounds(1) = mcmc_object.model.dt;
end
mcmc_object.bounds.t0 = t0_bounds;

% set crevasse height bounds- will do nothing if using erf or gaussian stf
mcmc_object.bounds.h_c_initial = h_c_initial_bounds;
mcmc_object.bounds.h_c_final = h_c_final_bounds;

end
