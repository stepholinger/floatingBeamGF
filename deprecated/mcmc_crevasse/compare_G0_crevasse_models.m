% this code runs a version of P. Segall's mcmc code adapted for fitting the beam deflection
% Green's function with observed seismic data

% set output path
path = "/home/setholinger/Documents/Projects/PIG/modeling/mcmc/";

% set some parameters
statDist = 21000;
t0 = 8.75;
h_i_avg = 400;
h_w_avg = 590;
f_max = 1;
t_max = 250;

% set crevasse height stuff- these values range from -0.5 to 0.5
h_c_initial = -0.05;
h_c_final = 0.05;

% construct some variables
t = 1/(2*f_max):1/(2*f_max):t_max;
nt = t_max*(2*f_max);

% get real data
centroid = 0;
numCluster = 2;
centroid_fs = 2;
if centroid
    centroids = h5read(numCluster + "_cluster_predictions_0.05-1Hz.h5","/centroids");
    centroids = squeeze(centroids)';
    centroids = centroids(:,1:1001);
    %medAmps = h5read(numCluster + "_cluster_median_amplitudes.h5","/median_amplitudes");
    %centroid_win = [400,900;1,500;100,600;100,300;200,700;200,700;400,900;0,500;0,500;400,900];
    centroid_win = [100,600;100,600];
else
    
    fname = "/media/Data/Data/PIG/MSEED/noIR/PIG2/HHZ/2012-06-20.PIG2.HHZ.noIR.MSEED";
    dataStruct = rdmseed(fname);

    % extract trace
    trace = extractfield(dataStruct,'d');
    fs = 100;

    % resample data to 1 Hz
    fsNew = f_max*2;
    trace = resample(trace,fsNew,100);

    % set event bounds
    hr = 20;
    min = 42;
    sec = 0;
    startTime = ((hr*60+min)*60+sec)*fsNew;
    nt = t_max*(2*f_max);
    endTime = startTime + nt;

    % trim data to event bounds
    eventTrace = trace(startTime:endTime-1);

    % remove scalar offset using first value
    eventTrace = eventTrace - eventTrace(1);
    
    % find index of max value
    [~,dataMaxIdx] = max(eventTrace);
end

% get current centroid and scale by median cluster amplitude
if centroid

    eventTrace = centroids(p,centroid_win(p,1):centroid_win(p,2));
    eventTrace = eventTrace/max(abs(eventTrace));
    eventTrace = eventTrace*medAmps(p);

    % resample centroid
    if f_max < 1
        centroid_fs = centroid_fs/f_max;
        eventTrace = resample(eventTrace,2,centroid_fs);
    else
        eventTrace = resample(eventTrace,f_max*2,centroid_fs);
    end
    % find index of max value
    [~,dataMaxIdx] = max(eventTrace);
end

% trim data if needed
t = 1/(2*f_max):1/(2*f_max):t_max;
nt = t_max*(2*f_max);
eventTraceTrim = eventTrace(1:nt);

% generate intial Green's function for each crevasse mode
x0 = [h_i_avg,h_w_avg,statDist,t0,h_c_initial,h_c_final,f_max,t_max];
[G_0_basal,corrCoef] = GF_func_mcmc(x0,"basal",eventTraceTrim);

x0 = [h_i_avg,h_w_avg,statDist,t0,h_c_final,h_c_initial,f_max,t_max];
[G_0_surface,corrCoef] = GF_func_mcmc(x0,"surface",eventTraceTrim);

x0 = [h_i_avg,h_w_avg,statDist,t0,h_c_initial,h_c_final,f_max,t_max];
[G_0_hydrostatic,corrCoef] = GF_func_mcmc(x0,"hydrostatic",eventTraceTrim);

