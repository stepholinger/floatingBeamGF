% script to launch multiple runs of run_GF_mcmc with different parameter files in parallel 

% make a list of parameter files
param_file_list = ["hydrostatic_moment_big_shear_t0_only"];

% make a list of output files
out_path = "/home/setholinger/Documents/Projects/PIG/modeling/mcmc/final_mcmc/";
out_file_list = param_file_list;

errors = strings(length(param_file_list),2);
parfor p = 1:length(param_file_list)
    try
        run_GF_mcmc(param_file_list(p),out_path + out_file_list(p))
    catch ME
        errors(p,:) = [param_file_list(p),getReport(ME)]; 
    end
end