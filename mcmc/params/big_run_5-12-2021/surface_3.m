function [L, f_max, t_max, h_i, h_w, stat_dist, t0, source, mode, h_c_initial, h_c_final,...
          h_i_bounds, h_w_bounds, stat_dist_bounds, t0_bounds, h_c_initial_bounds, h_c_final_bounds,...
          h_i_step, h_w_step, stat_dist_step, t0_step, h_c_initial_step, h_c_final_step,...
          sigma, num_iterations, L_type, data_path, data_set] = surface_3()

% set fixed model parameters
L = 1e7;
f_max = 2;
t_max = 1000;
source = "moment";
mode = "surface";

% set mcmc inital parameters
h_i = 400;
h_w = 590;
stat_dist = 20000;
t0 = 3;
h_c_initial = 0.01;
h_c_final = -0.01;

% set mcmc bounds for each parameter
h_i_bounds = [0,1000];
h_w_bounds = [0,1000];
stat_dist_bounds = [0,40000];
t0_bounds = [0,20];
h_c_initial_bounds = [-0.5,0.5];
h_c_final_bounds = [-0.5,0.5];

% set mcmc step size for each parameter
h_i_step = 0;
h_w_step = 0;
stat_dist_step = 1000;
t0_step = 0.5;
h_c_initial_step = 0.025;
h_c_final_step = 0.025;

% set mcmc run parameters
sigma = 100;
num_iterations = 100000;
L_type = "xcorr";

% set paths to clustering results to get average waveform
num_cluster = 2;
centroid_number = 1;
data_set = "/mean_wave";
data_path = "/home/setholinger/Documents/Projects/PIG/detections/templateMatch/multiTemplate/run3/short_3D_clustering/modified_k_shape/";
data_path = data_path + num_cluster + "/cluster_" + string(centroid_number) + "_mean_wave_0.05-1Hz.h5";

end