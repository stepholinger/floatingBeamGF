function [G_stf,G,stf] = conv_GF(mcmc_object,G)

% get current iteration number and parameters
it = mcmc_object.iteration;
h_i = mcmc_object.results.x_keep(1,it);
h_w = mcmc_object.results.x_keep(2,it);
stat_dist = mcmc_object.results.x_keep(3,it);
t0 = mcmc_object.results.x_keep(4,it);
if sum(contains(["basal","surface","hydrostatic"],mcmc_object.mode)) > 0 
    h_c_initial = mcmc_object.results.x_keep(5,it);
    h_c_final = mcmc_object.results.x_keep(6,it);
end

% caclulate max pressure and moment
P_max = mcmc_object.model.rho_i * mcmc_object.model.g * mcmc_object.model.h_i;
M_max =  mcmc_object.model.rho_i * mcmc_object.model.g * mcmc_object.model.h_i^3 / 12 * 0.072;

% get source time function and convolve with Green's function
if sum(contains(["basal","surface","hydrostatic"],mcmc_object.mode)) > 0 
    
    % call stf function with crevasse mode
    stf = source_time_function(mcmc_object.model,t0,mcmc_object.mode,h_c_initial,h_c_final);
    
    % convolve stf with Green's function
    G_pad = zeros(1,length(stf));
    G_pad(ceil(end/2):ceil(end/2)+length(G)-1) = G;
    G_stf_pad = ifft(fft(G_pad).*fft(stf));
    
    % get rid of padding
    G_stf = G_stf_pad;
    
elseif mcmc_object.mode == "erf"
    
    % scale Green's function output
    if mcmc_object.source == "moment"
        G = G * M_max;
    elseif mcmc_object.source == "force"
        G = G * P_max;
    end    
    
    % call stf function with erf mode
    stf = source_time_function(mcmc_object.model,t0,mcmc_object.mode);
    
    % convolve stf with Green's function
    G_pad = zeros(1,length(stf));
    G_pad(ceil(end/2):ceil(end/2)+length(G)-1) = G;
    G_stf_pad = ifft(fft(G_pad).*fft(stf));
    
    % remove padding
    G_stf = G_stf_pad(ceil(end/2):end);
    
elseif mcmc_object.mode == "gaussian"
    
    % scale Green's function output
    if mcmc_object.source == "moment"
        G = G * M_max;
    elseif mcmc_object.source == "force"
        G = G * P_max;
    end    
    
    % call stf function with gaussian mode
    stf = source_time_function(mcmc_object.model,t0,mcmc_object.mode);
    
    % convolve stf with Green's function
    G_pad = zeros(1,length(stf));
    G_pad(ceil(end/2):ceil(end/2)+length(G)-1) = G;
    G_stf_pad = ifft(fft(G_pad).*fft(stf));
    
    % get rid of padding
    G_stf = G_stf_pad(ceil(end/2):end);
    
end
    
end