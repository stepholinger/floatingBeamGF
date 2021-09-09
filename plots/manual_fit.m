% strategy is to match wave shape with amplitudes normalized to 1, then use
% scaling factor to match amplitude


% set model parameters
L = 1e6;
f_max = 2;
t_max = 1000;
h_i = 400;
h_w = 590;
statDist = 29000;
t0 = 2;
source = "moment";
mode = "erf";
h_c_initial = -0.12;
h_c_final = 0.05;

% set length of viewing window in seconds 
view_len = 1000;

% run model
[model,dGdt,G,stf] = calcGF(L,f_max,t_max,h_i,h_w,statDist,t0,source,mode,h_c_initial,h_c_final);

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
    
    fname = "/media/Data/Data/PIG/MSEED/noIR/PIG2/HHZ/2013-04-15.PIG2.HHZ.noIR.MSEED";
    dataStruct = rdmseed(fname);

    % extract trace
    trace = extractfield(dataStruct,'d');
    fs = 100;

    % resample data to 1 Hz
    fsNew = f_max*2;
    trace = resample(trace,fsNew,100);

    % set event bounds
    event_len = 500;
    hour = 11;
    minute = 30;
    second = 53;
    startTime = ((hour*60+minute)*60+second)*fsNew;
    nt = event_len*(2*f_max);
    endTime = startTime + nt;

    eventTrace = trace(startTime:endTime-1);

end

% remove scalar offset using first value
eventTrace = eventTrace - eventTrace(1);

% normalize before cross correlation if not using crevasse stf
if sum(contains(["erf","gaussian"],mode)) > 0
    if source == "force" 
        P_frac = max(eventTrace)/max(dGdt);
        dGdt = dGdt * P_frac;
    elseif source == "moment"
        M_frac = max(eventTrace)/max(dGdt);
        dGdt = dGdt * M_frac;
    end
end

% cross correlate with data
eventTrace = [zeros(1,event_len),eventTrace,zeros(1,length(dGdt)-length(eventTrace)-event_len)];
[xcorr_trace,lag] = xcorr(eventTrace,dGdt,"coeff");
[xcorr_coef,lag_idx] = max(xcorr_trace);  

slide = lag(lag_idx);

if slide > 0
    dGdt = [zeros(1,abs(slide)),dGdt];
    dGdt = dGdt(1:length(eventTrace));
end
if slide < 0
    dGdt = dGdt(abs(slide):end);
    if length(dGdt) < length(eventTrace)
        dGdt = [dGdt,zeros(1,length(eventTrace)-length(dGdt))];
    end
end

% make plot
%figure(1)
%plot(model.t(1:viewLen*fsNew),dGdt(1:viewLen*fsNew))
%hold on;
%plot(model.t(1:viewLen*fsNew),eventTrace(1:viewLen*fsNew))

%figure(1)
%plot(dGdt(1:length(eventTrace)))
%hold on;
%plot(eventTrace)

% make normalized plot
figure(2)
plot(model.t(1:view_len*fsNew),dGdt(1:view_len*fsNew))
hold on;
plot(model.t(1:view_len*fsNew),eventTrace(1:view_len*fsNew))
xlabel('Time (seconds)');
ylabel('Velocity (m/s)');
%title("P = " + num2str(P_frac) + "*rho_i*g*h")