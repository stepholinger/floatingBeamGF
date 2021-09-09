function mcmc_object = set_mcmc_steps(mcmc_object,h_i_step,h_w_step,stat_dist_step,t0_step,h_c_initial_step,h_c_final_step)
% set bounds for mcmc
mcmc_object.steps.h_i = h_i_step;
mcmc_object.steps.h_w = h_w_step;
mcmc_object.steps.stat_dist = stat_dist_step;
mcmc_object.steps.t0 = t0_step;

% set crevasse height bounds- will do nothing if using erf or gaussian stf
mcmc_object.steps.h_c_initial = h_c_initial_step;
mcmc_object.steps.h_c_final = h_c_final_step;

end
