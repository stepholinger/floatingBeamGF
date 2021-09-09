function mcmc_object = mcmc(mcmc_object)

% check to make sure input is a function
fun = fcnchk(mcmc_object.function);

%number of elements in x
num_params = length(mcmc_object.results.x_keep(:,1));

% set a few parameters
x = mcmc_object.results.x_keep(:,1);
x_step = [mcmc_object.steps.h_i;mcmc_object.steps.h_w;mcmc_object.steps.stat_dist;mcmc_object.steps.t0;mcmc_object.steps.h_c_initial;mcmc_object.steps.h_c_final];
x_bounds = [mcmc_object.bounds.h_i;mcmc_object.bounds.h_w;mcmc_object.bounds.stat_dist;mcmc_object.bounds.t0;mcmc_object.bounds.h_c_initial;mcmc_object.bounds.h_c_final];
L = mcmc_object.L_0;
H = NaN;

% start iteration
try
    for it = 2:mcmc_object.num_iterations
        
        disp(' ');
        disp(['Starting iteration ' num2str(it) ' of ' num2str(mcmc_object.num_iterations)]);

        % generate proposal
        x_prop = x + x_step.* 2 .* (rand(num_params,1)-0.5);
        
        % update mcmc object
        mcmc_object.iteration = it;
        mcmc_object.results.x_keep(:,it) = x_prop;

        % check bounds for all parameters
        if check_bounds(mcmc_object,x_prop,x_bounds) == 1

            % generate prediction for proposal
            [d_Prop,xcorr_coef,M_frac] = fun(mcmc_object);

            % calculate likelihood      
            L_prop = liklihood(d_Prop,mcmc_object,xcorr_coef);
            %fprintf("Lprop: " + Lprop + ", " + "L0: " + L0 + "\n")

            % make random number to compare to ratio of likelihoods
            u = rand;

            % compute hastings ratio
            H = exp(L_prop-L);

            % accept proposal 
            if (L == 0 || u <= min([H,1]))
                 fprintf("Accepted proposal ( " + u + " < " + min([H,1]) +")\n")
                 x = x_prop;
                 L = L_prop;
                 mcmc_object.results.accept(it) = 1;
            end 
            
        end

        % save results
        mcmc_object.results.x_keep(:,it) = x;
        mcmc_object.results.L_keep(it) = L;
        mcmc_object.results.alpha_keep(it) = H;
        
        % save scaling value if non-crevasse stf
        if sum(contains(["erf","gaussian"],mcmc_object.mode)) > 0
            mcmc_object.results.M_frac(it) = M_frac;
        end      
        
    end
    
catch ME
    mcmc_object.errors.x = x_prop;
    mcmc_object.errors.message = getReport(ME);
    mcmc_object.errors.iteration = it;

end