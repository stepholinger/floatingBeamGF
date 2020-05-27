% set output path
path = "/home/setholinger/Documents/Projects/PIG/modeling/mcmc/RC/";

% set data parameters
statDist = 10000;
t0 = 5.5;
f_max = 10;
t_max = 250;
t = 1/(2*f_max):1/(2*f_max):t_max;
nt = t_max*(2*f_max);

% set flag for toy problem
test = 0;

% get real data
fname = "/media/Data/Data/PIG/MSEED/noIR/PIG2/HHZ/2012-04-02.PIG2.HHZ.noIR.MSEED";
dataStruct = rdmseed(fname);

% extract trace
trace = extractfield(dataStruct,'d');
fs = 100;

% resample data to 1 Hz
fsNew = f_max*2;
trace = resample(trace,fsNew,100);

% set event bounds
startTime = ((15*60+18)*60+50)*fsNew;
endTime = startTime + nt;

% trim data to event bounds
eventTrace = trace(startTime:endTime-1);

% remove scalar offset using first value
eventTrace = eventTrace - eventTrace(1);

% if toy problem mode, set waveform to recover
if test
    testParams = [450,600,statDist,t0,f_max,t_max];
    [eventTrace,~] = GF_func_mcmc(testParams,eventTrace);
end

% find index of max value
[~,dataMaxIdx] = max(eventTrace);

% set mcmc parameters
% x0 goes like this: [h_i,h_w,statDist,t0,f_max,t_max,dataMaxIdx]
% f_max, t_max, and dataMaxIdx MUST have 0 step size in xStep
% t0 is log now! so the value here is like 10^x
x0_vect = {[350,800,10000,1,f_max,t_max]};
xStepVect = {[100,100,1000,log10(50),0,0]};
xBounds = [0,1000;
           0,1000;
           0,100000;
           -1,2;
           0,f_max+1;
           0,t_max+1;];
sigmaVect = [16];
t_max_vect = [t_max];
numIt = 1000;
L_type_vect = ["modified"];
axisLabels = ["Ice thickness (m)", "Water depth (m)", "X_{stat} (km)","t_0 (s)"];
paramLabels = ["h_i","h_w","Xstat","t0"];
maxNumBins =  100;

try
    parpool;
    poolobj = gcp;
catch
    fprintf("Using existing parpool...\n")
end

tic;

%parfor p = 1:length(sigmaVect)
for p = 1:length(sigmaVect)      
    
    % get parameters for run
    xStep = xStepVect{p};
    sigma = sigmaVect(p);
    L_type = L_type_vect(p);
    t_max = t_max_vect(p);
    x0 = x0_vect{p};

    % trim data if needed
    t = 1/(2*f_max):1/(2*f_max):t_max;
    nt = t_max*(2*f_max);
    eventTraceTrim = eventTrace(1:nt);
    
    % record which two parameters will be varied this run
    paramInd = [1,2,3,4,5,6];
    paramsVaried = paramInd(xStep ~= 0);
    
    % deal with log t0
    x0(4) = 10^(x0(4));
    
    % generate intial Green's function
    [G_0,eventAlign,M_frac_0] = GF_func_mcmc(x0,eventTraceTrim);

    % deal with log t0
    x0(4) = log10(x0(4));
    
    % calculate initial liklihood
    L0 = liklihood(G_0,eventAlign,sigma,L_type);
    
    % run mcmc
    [x_keep,L_keep,count,alpha_keep,accept,M_frac] = mcmc('GF_func_mcmc',eventTraceTrim,...
                                              x0,xStep,xBounds,sigma,numIt,M_frac_0,L0,L_type);
                                          
    % give output
    fprintf("Accepted " + round((sum(accept)/numIt)*100) + " %% of proposals\n");
    
    % make plots for bivariate runs
    if length(paramsVaried) == 2

        % find best-fit parameters using 2D histogram
        numBins = length(unique(x_keep(paramsVaried(1),:)));
        if numBins > maxNumBins
            numBins = maxNumBins;
        end
        
        % get best fit parameters
        xFit = getFit(x_keep,paramsVaried,numBins,x0);

        % run model for best fit parameters
        [G_fit,eventAlign,M_fit] = GF_func_mcmc(xFit,eventTraceTrim);

        % calculate solution liklihood
        L_fit = liklihood(G_fit,eventAlign,sigma,L_type);

        % convert X_stat to km
        x_keep(3,:) = x_keep(3,:)/1000;
        
        % 'unlog' t0
        x_keep(4,:) = 10.^x_keep(4,:);
        
        % save results
        resultStruct = struct('xFit',xFit,'L_fit',L_fit,'G_fit',G_fit,'G_0',G_0,'L_keep',L_keep,...
                              'x_keep',x_keep,'x0',x0,'xStep',xStep,'M_frac',M_frac,'L_type',L_type,...
                              'xBounds',xBounds,'sigma',sigma,'numIt',numIt,'labels',paramLabels);
        parsave("run" + p + "_results.mat",resultStruct)   
        
        % call plotting functions
        plot_bivar(x_keep,xFit,numIt,p,paramsVaried,axisLabels,paramLabels,numBins,path)
        plot_M_frac(x_keep,M_frac,M_fit,xFit,numIt,p,paramsVaried,axisLabels,paramLabels,numBins,path)
        plot_start_wave(t,eventAlign,sigma,L0,M_frac_0,G_0,x0,numIt,xStep,p,path)
        plot_fit_wave(t,eventAlign,sigma,L_fit,M_fit,G_fit,xFit,numIt,xStep,p,accept,L_type,path)

    else

        % convert X_stat to km
        x_keep(3,:) = x_keep(3,:)/1000;
        
        % 'unlog' t0
        x_keep(4,:) = 10.^x_keep(4,:);
        
        % save results
        resultStruct = struct('G_0',G_0,'L_keep',L_keep,'x_keep',x_keep,'x0',x0,...
                              'xStep',xStep,'M_frac',M_frac,'xBounds',xBounds,'L_type',L_type,...
                              'sigma',sigma,'numIt',numIt,'labels',paramLabels,'accept',accept,...
                              'f_max',f_max,'t_max',t_max);
        parsave(path + "run" + p + "_results.mat",resultStruct)
    
        % call plotting functions
        plot_multivar(sigma,accept,xStep,x_keep,M_frac,x0,numIt,....
                      p,paramsVaried,axisLabels,maxNumBins,L_type,path,f_max,t_max)
    end
    
end

runtime = toc;