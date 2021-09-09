function [dGdt,xcorr_coef,M_frac] = GF_func_mcmc(mcmc_object)

% this function is for use with P. Segall's MCMC code- instead of actually
% calling the Green's function model, it pulls already created Green's
% functions (spatial derivative has already been taken). Therefore, L, f_max, and t_max are fixed, and already set by
% this point. In general, they will be L = 1e7, f_max = 1, t_max = 1000

% get current iteration number and parameters
it = mcmc_object.iteration;
stat_dist = mcmc_object.results.x_keep(3,it);

% get index of desired position
[~,locIdx] = min(abs(mcmc_object.model.x - stat_dist));

% get trace closest seismometer location
G = mcmc_object.G_matrix(locIdx,:);

% pass matrix of Green's function solutions to convGF for stf convolution
[G_stf,~,~] = conv_GF(mcmc_object,G);

% take time derivative to get velocity seismogram
dGdt = gradient(G_stf,mcmc_object.model.dt);

% cross correlate with data
dGdt = dGdt(1:10*length(mcmc_object.data));
[xcorr_trace,lag] = xcorr([mcmc_object.data,zeros(1,9*length(mcmc_object.data))],dGdt,"coeff");
[xcorr_coef,lag_idx] = max(xcorr_trace);  

slide = lag(lag_idx);

if slide > 0
    dGdt = [zeros(1,abs(slide)),dGdt];
    dGdt = dGdt(1:length(mcmc_object.data));
end
if slide < 0
    dGdt = dGdt(abs(slide):end);
    if length(dGdt) < length(mcmc_object.data)
        dGdt = [dGdt,zeros(1,length(mcmc_object.data)-length(dGdt))];
    elseif length(dGdt) > length(mcmc_object.data)
        dGdt = dGdt(1:length(mcmc_object.data));
    end
end

% get scaling value if non-crevasse stf
if sum(contains(["erf","gaussian"],mcmc_object.mode)) > 0
    M_frac = max(abs(mcmc_object.data))/max(abs(dGdt));
    dGdt = dGdt * M_frac;
else
    M_frac = NaN;
end

end