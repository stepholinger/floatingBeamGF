function L = liklihood(G,mcmc_object,xcorr_coef)
if mcmc_object.L_type == "standard" || mcmc_object.L_type == "xcorr"
    normArg = (mcmc_object.data-G)./mcmc_object.data;

elseif mcmc_object.L_type == "modified"
    normArg = (mcmc_object.data-G);    
end

normArg(isnan(normArg)) = 0;
L = -0.5/mcmc_object.sigma^2 * norm(normArg)^2;

if mcmc_object.L_type == "xcorr" && nargin > 2
    L = L/xcorr_coef;
end

end