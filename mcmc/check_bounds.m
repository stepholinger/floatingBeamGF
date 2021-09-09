function check = check_bounds(mcmc_object,x_prop,x_bounds)

% if using a crevasse source time function, verify that inital and final
% crevasse heights are sensical
if mcmc_object.mode == "basal" || mcmc_object.mode == "hydrostatic"
    stf_check = logical(x_prop(5) < x_prop(6));
elseif mcmc_object.mode == "surface"
    stf_check = logical(x_prop(5) > x_prop(6));
else
    stf_check = 1;
end

% make sure all proposed parameters are greater than lower bound and less
% than upper bound
prop_check = logical(min(x_prop > x_bounds(:,1)) & min(x_prop < x_bounds(:,2)));

check = (stf_check + prop_check)/2;

end