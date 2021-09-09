% code to make a multi-panel figure comparing source time functions and
% moment vs single force for supplement.

% set model parameters
L = 1e7;
f_max = 2;
t_max = 1000;
h_i = 400;
h_w = 590;
statDist = [18000,21000,16000,19000,21500];
t0_vect = [2,3,2,3,5];
sources = ["force","force","moment","moment","moment"];
modes = ["erf","gaussian","erf","gaussian","hydrostatic"];
h_c_initial = -0.40;
h_c_final = -0.35;

% set length of viewing window in seconds 
view_len = 1000;

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
    
    fname = "/media/Data/Data/PIG/MSEED/noIR/PIG2/HHZ/2012-03-26.PIG2.HHZ.noIR.MSEED";
    dataStruct = rdmseed(fname);

    % extract trace
    trace = extractfield(dataStruct,'d');
    fs = 100;

    % resample data to 1 Hz
    fsNew = f_max*2;
    trace = resample(trace,fsNew,100);

    % set event bounds
    event_len = 1000;
    hour = 8;
    minute = 54;
    second = 26;
    startTime = ((hour*60+minute)*60+second)*fsNew;
    nt = event_len*(2*f_max);
    endTime = startTime + nt;

    data = trace(startTime:endTime-1);

end

% make normalized plot
for n = 1:length(sources)*3
    disp(n)
    if n == 1
        mode = modes(1);
        source = sources(1);
        t0 = t0_vect(1);
        [model,dGdt,G,stf] = calcGF(L,f_max,t_max,h_i,h_w,statDist(1),t0,source,mode);
    
    elseif n == 4
        mode = modes(2);
        source = sources(2);
        t0 = t0_vect(2);
        [model,dGdt,G,stf] = calcGF(L,f_max,t_max,h_i,h_w,statDist(2),t0,source,mode);

    elseif n == 7
        mode = modes(3);
        source = sources(3);
        t0 = t0_vect(3);
        [model,dGdt,G,stf] = calcGF(L,f_max,t_max,h_i,h_w,statDist(3),t0,source,mode);
        
    elseif n == 10
         mode = modes(4);
         source = sources(4);
         t0 = t0_vect(4);
         [model,dGdt,G,stf] = calcGF(L,f_max,t_max,h_i,h_w,statDist(4),t0,source,mode);

    elseif n == 13
         mode = modes(5);
         source = sources(5);
         t0 = t0_vect(5);
         [model,dGdt,G,stf] = calcGF(L,f_max,t_max,h_i,h_w,statDist(5),t0,source,mode,h_c_initial,h_c_final);

    end
    
    % remove scalar offset using first value
    eventTrace = data;
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

    subplot(5,3,n)
    
    
    if sum([1,4,7,10,13]==n) > 0 
        plot(model.t(1:view_len*fsNew),G(1:view_len*fsNew))
        title("G (" + source + ")")
    end
    
    if sum([2,5,8,11,14]==n) > 0 
        stf_plot_len = t0*f_max*2*2;
        if mode == "hydrostatic"
            plot(model.t(1:1+stf_plot_len+100),stf(floor(end/2)-100:floor(end/2)+stf_plot_len)) 
        elseif mode == "erf"
            plot(model.dt:model.dt:stf_plot_len*10*model.dt,stf(1:stf_plot_len*10)) 
        elseif mode == "gaussian"
            [~,offset_index] = find(stf > 0.01);
            plot(stf(offset_index)) 
        end
        title("stf (" + mode + ")")
    end
    
    if mod(n,3) == 0
        plot(model.t(1:view_len*fsNew),dGdt(1:view_len*fsNew))
        hold on;
        plot(model.t(1:view_len*fsNew),eventTrace(1:view_len*fsNew))
        title("G*stf")
        xlim([0,view_len])
    end
    
    % each row corresponds to one Green's function / stf combo

end

%xlabel('Time (seconds)');
%ylabel('Velocity (m/s)');
%title("P = " + num2str(P_frac) + "*rho_i*g*h")