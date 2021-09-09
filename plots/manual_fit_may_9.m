% strategy is to match wave shape with amplitudes normalized to 1, then use
% scaling factor to match amplitude


% set model parameters
L = 1e7;
f_max = 2;
t_max = 1000;
h_i = 400;
h_w = 590;
statDist = 20000;
t0 = 500;
source = "moment";
mode = "hydrostatic";
h_c_initial = -0.5;
h_c_final = -0.5+0.8945;

% set length of viewing window in seconds 
viewLen = 7200;

% run model
[model,dGdt,G,stf] = calcGF(L,f_max,t_max,h_i,h_w,statDist,t0,source,mode,h_c_initial,h_c_final);

% find index of max value
[modelMax,modelMaxIdx] = max(dGdt);

% get average waveform
centroid = 0;
numCluster = 2;
centroid_number = 1;
fsNew = 2;
centroidPath = "/home/setholinger/Documents/Projects/PIG/detections/templateMatch/multiTemplate/run3/short_3D_clustering/modified_k_shape/";

if centroid
    centroids = h5read(centroidPath + numCluster + "/cluster_" + string(centroid_number) + "_mean_wave_0.05-1Hz.h5","/shear_mean_wave");
    centroids = squeeze(centroids)';
    centroids = centroids(:,200:600);
    eventTrace = centroids;

else
    
    % get real data
    fname = "/media/Data/Data/PIG/MSEED/noIR/PIG2/HHZ/2012-05-09.PIG2.HHZ.noIR.MSEED";
    dataStruct = rdmseed(fname);
    
    % extract trace
    trace = extractfield(dataStruct,'d');
    fs = 100;

    % resample data to 1 Hz
    fsNew = f_max*2;
    if f_max < 0.5
        fs = fs/(f_max*2);
        trace = resample(trace,1,fs);
    else
        trace = resample(trace,fsNew,fs);
    end

    % filter out frequencies above the gaussian smoothing frequency (t0/2)
    % choose filter band and design filter
    freq = 1/(t0/2);
    [b,a] = butter(4,freq/(fsNew/2),'low');
    trace = filtfilt(b,a,trace);
    
    % set event bounds
    startTime = ((18*60)*60)*fsNew;
    endTime = startTime + viewLen*1/model.dt;
    if endTime > length(trace)
        eventTrace = trace(startTime:end);
        eventTrace = [eventTrace,zeros(1,endTime-startTime-length(eventTrace))];
    else
        % trim data to event bounds
        eventTrace = trace(startTime:endTime-1);
    end
end

% remove scalar offset using first value
eventTrace = eventTrace - eventTrace(1);

% find index of max value
[dataMax,dataMaxIdx] = max(eventTrace);

% align maximum values by padding with zeros
if modelMaxIdx > dataMaxIdx
    slide = modelMaxIdx-dataMaxIdx;
    eventTrace = [zeros(1,slide),eventTrace(1:end-slide)];
else
    slide = dataMaxIdx-modelMaxIdx;
    dGdt = [zeros(1,slide),dGdt(1:end-slide)];
end    

% make plot
%figure(1)
%plot(model.t(1:viewLen*fsNew),dGdt(1:viewLen*fsNew))
%hold on;
%plot(model.t(1:viewLen*fsNew),eventTrace(1:viewLen*fsNew))

figure(1)
plot(dGdt(1:length(eventTrace)))
hold on;
plot(eventTrace)

% make normalized plot
%figure(2)
%plot(model.t(1:viewLen*fsNew),dGdt(1:viewLen*fsNew)/max(abs(dGdt(1:viewLen*fsNew))))
%hold on;
%plot(model.t(1:viewLen*fsNew),eventTrace(1:viewLen*fsNew)/max(abs(eventTrace(1:viewLen*fsNew))))