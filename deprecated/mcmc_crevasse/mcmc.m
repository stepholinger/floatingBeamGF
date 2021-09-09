function  [x_keep,L_keep,count,alpha_keep,accept] = mcmc(func,data,x0,xStep,xBounds,...
                 sigma,numIt,L0,liklihoodType)

% [x_keep, L_keep, count] = mcmc(func,data,x0,xstep,sigma,Niter,varargin)
%
% subroutine for MCMC sampling using Metropolis-Hasting w/ normal
% distribution.
%
% Inputs:
%       func  = function name eg:   mcmc('travel_time',....)
%               parameter required for function in varargin
%       data  = vector of observations
%       x0    = initial estimate of parameter vector
%       xStep = step size in all parameter directions
%       xBnds = bounds (
%       sigma = sigma of normal distribution 
%       numIt = number of iterations
%
% Outputs:
%       x_keep = array of samples
%       L_keep = likelihood of samples
%       count  = number of accepted. Acceptance ratio is count/Niter
%
%  P Segall: 2012
% TODO? Should I add something to seed rand? Better if this is done in
% calling program
% rng('shuffle');

% check to make sure input is a function
fun  = fcnchk(func);

%number of elements in x
numParams = length(x0);

% check dimensions of bounds
if( size(xBounds,1) ~= numParams || size(xBounds,2) ~= 2)
    disp('Dimension of xBounds is not valid')
    return
end

% make empty arrays for storing results
x_keep = zeros(numParams,numIt,'single'); 
L_keep = zeros(1,numIt,'single'); 
alpha_keep = zeros(1,numIt,'single'); 
accept = zeros(1,numIt,'single');

% set a few parameters
x = x0;
L = L0;
H = NaN;

% make counter and start iteration
count = 0;
for k = 1:numIt
disp(' ');
disp(['Starting iteration ' num2str(k) ' of ' num2str(numIt)]);
    
    % generate proposal
    xProp = x + xStep.* 2 .* (rand(1,numParams)-0.5);
    
    % make sure h_c_final is larger than h_c_initial
    if xProp(6) > xProp(5)

        % check bounds 
        if (min(xProp' > xBounds(:,1)) & min(xProp' < xBounds(:,2)))

            % deal with log t0
            %xProp(4) = 10^(xProp(4));

            % generate prediction for proposal
            [dProp,corrCoef] = fun(xProp,data);

            % deal with log t0
            %xProp(4) = log10(xProp(4));

            % calculate likelihood      
            Lprop = liklihood(dProp,data,sigma,liklihoodType,corrCoef);
            %fprintf("Lprop: " + Lprop + ", " + "L0: " + L0 + "\n")

            % make random number to compare to ratio of likelihoods
            u = rand;

            % compute hastings ratio
            H = exp(Lprop-L);

            % accept proposal 
            if (L == 0 || u <= min([H,1]))
                 fprintf("Accepted proposal ( " + u + " < " + min([H,1]) +")\n")
                 count = count+1;
                 x = xProp;
                 L = Lprop;
                 accept(k) = 1;
            end 
        end
    end
    
    % save results
    x_keep(:,k) = x;
    L_keep(k) = L;
    alpha_keep(k) = H;
    
end

end