function mcmc_object = centroid_mcmc_setup(param_func)

% read parameter file
[L, f_max, t_max, h_i, h_w, stat_dist, t0, source, mode, h_c_initial, h_c_final,...
h_i_bounds, h_w_bounds, stat_dist_bounds, t0_bounds, h_c_initial_bounds, h_c_final_bounds,...
h_i_step, h_w_step, stat_dist_step, t0_step, h_c_initial_step, h_c_final_step,...
sigma, num_iterations, L_type, data_path, data_set] = param_func();  

% make model object
model = loadParameters(L,f_max,t_max,h_i,h_w);

% read in k-shape mean event
data = read_centroid_mcmc(data_path,data_set);

% make MCMC object, set MCMC initial parameters, bounds, and step sizes
mcmc_object = load_mcmc_object(model,data,stat_dist,t0,source,mode,h_c_initial,h_c_final);
mcmc_object = set_mcmc_bounds(mcmc_object,h_i_bounds,h_w_bounds,stat_dist_bounds,t0_bounds,h_c_initial_bounds,h_c_final_bounds);
mcmc_object = set_mcmc_steps(mcmc_object,h_i_step,h_w_step,stat_dist_step,t0_step,h_c_initial_step,h_c_final_step);

% set MCMC run parameters
mcmc_object.sigma = sigma;
mcmc_object.num_iterations = num_iterations;
mcmc_object.L_type = L_type;
mcmc_object.function = "GF_func_mcmc";
mcmc_object.iteration = 1;

% make fields for MCMC results
mcmc_object.results.x_keep = zeros(length(fieldnames(mcmc_object.x0)),num_iterations);
mcmc_object.results.alpha_keep = zeros(1,num_iterations);
mcmc_object.results.accept = zeros(1,num_iterations);
mcmc_object.results.L_keep = zeros(1,num_iterations);

% field for scaling value if non-crevasse stf
if sum(contains(["erf","gaussian"],mcmc_object.mode)) > 0
    mcmc_object.results.M_frac = zeros(1,num_iterations);
end

% fill first column of x_keep with initial parameters
mcmc_object.results.x_keep(1,1) = h_i;
mcmc_object.results.x_keep(2,1) = h_w;
mcmc_object.results.x_keep(3,1) = stat_dist;
mcmc_object.results.x_keep(4,1) = t0;
mcmc_object.results.x_keep(5,1) = h_c_initial;
mcmc_object.results.x_keep(6,1) = h_c_final;
