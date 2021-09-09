function [L, f_max, t_max, h_i, h_w, stat_dist, t0, source, mode, h_c_initial, h_c_final,...
          h_i_bounds, h_w_bounds, stat_dist_bounds, t0_bounds, h_c_initial_bounds, h_c_final_bounds,...
          h_i_step, h_w_step, stat_dist_step, t0_step, h_c_initial_step, h_c_final_step,...
          sigma, num_iterations, L_type, data_path, fs, event_bounds] = erf_moment_big_rift()

% set fixed model parameters
L = 1e7;
f_max = 2;
t_max = 1000;
source = "moment";
mode = "erf";

% set mcmc inital parameters
h_i = 400;
h_w = 590;
stat_dist = 29000;
t0 = 2;
h_c_initial = 0;
h_c_final = 0;

% set mcmc bounds for each parameter
h_i_bounds = [0,1000];
h_w_bounds = [0,1000];
stat_dist_bounds = [10000,40000];
t0_bounds = [0,20];
h_c_initial_bounds = [-0.5,0.5];
h_c_final_bounds = [-0.5,0.5];

% set mcmc step size for each parameter
h_i_step = 0;
h_w_step = 0;
stat_dist_step = 1000;
t0_step = 0.66;
h_c_initial_step = 0;
h_c_final_step = 0;

% set mcmc run parameters
sigma = 100;
num_iterations = 10000;
L_type = "xcorr";

% set path to data and event time and duration in seconds
fs = 100;
data_path = "/media/Data/Data/PIG/MSEED/noIR/PIG2/HHZ/2013-04-15.PIG2.HHZ.noIR.MSEED";
hour = 11;
minute = 30;
second = 53;
start_time = ((hour*60)+minute)*60+second;
duration = 1000;
event_bounds = [start_time,start_time + duration];

end