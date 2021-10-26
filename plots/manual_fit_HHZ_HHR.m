% set some parameters
f_max = 1;
source = "moment";
mode = "gaussian";
fs_data = 100;

% read vertical data
fname = "/media/Data/Data/PIG/MSEED/noIR/2013-03-26/2013-03-26T08:58:00.PIG3.HHZ.MSEED";
dataStruct = rdmseed(fname);

% extract trace
z_data = extractfield(dataStruct,'d');

% resample data to 1 Hz
fs_new = f_max*2;
z_data = resample(z_data,fs_new,fs_data);
z_data = z_data(1:end-1);

% read radial data
fname = "/media/Data/Data/PIG/MSEED/noIR/2013-03-26/2013-03-26T08:58:00.PIG3.HHR.MSEED";
dataStruct = rdmseed(fname);

% extract trace
r_data = extractfield(dataStruct,'d');

% resample data to 1 Hz
r_data = resample(r_data,fs_new,fs_data);
r_data = r_data(1:end-1);

% run model
[z_model,r_model,stf,model] = calcGF(1e7,1,840,400,590,28000,5.5,source,mode);
max_z_model = max(abs(z_model));
max_r_model = max(abs(r_model));

% normalize both model outputs to the max of both data (must preserve
% relative amplitudes) and align max amplitude phases
max_val = max(abs([z_data,r_data]));
if max(abs(z_model)) > max(abs(r_model))
   
    if contains(mode,["erf","gaussian"]) 
        % normalize both to max of data
        z_model = z_model./max_z_model.*max_val;
        r_model = r_model./max_z_model.*max_val;
    end
    
    % cross correlate appropriate traces
    [xcorr_trace,lag] = xcorr(z_data,z_model,"coeff");
    [xcorr_coef,lag_idx] = max(xcorr_trace);  
    slide = lag(lag_idx);
    
else
    
    if contains(mode,["erf","gaussian"]) 
        % normalize both to max of data
        z_model = z_model./max_r_model.*max_val;
        r_model = r_model./max_r_model.*max_val;
    end
    
    % cross correlate appropriate traces
    [xcorr_trace,lag] = xcorr(r_data,r_model,"coeff");
    [xcorr_coef,lag_idx] = max(xcorr_trace);  
    slide = lag(lag_idx);
end

% align both traces
if slide > 0
    z_model = [zeros(1,abs(slide)),z_model];
    z_model = z_model(1:length(z_data));
    r_model = [zeros(1,abs(slide)),r_model];
    r_model = r_model(1:length(r_data));

elseif slide < 0
    z_model = z_model(abs(slide):end);
    z_model = [z_model,zeros(1,length(z_data)-length(z_model))];
    r_model = r_model(abs(slide):end);
    r_model = [r_model,zeros(1,length(r_data)-length(r_model))];
end

% plot vertical data and model output together
subplot(2,1,1)
plot(model.t,z_data,'Color',[0.85,0.325,0.098],"LineWidth",1)
hold on;
plot(model.t,z_model,'Color',[0,0.447,0.741],"LineWidth",1)
xlabel("Time (s)");
ylabel("Velocity (m/s)")
ylim([-1*max_val,max_val]);
title("Vertical Component")
legend("Observed seismogram","Modeled seismogram")

% plot radial data and model output together
subplot(2,1,2)
plot(model.t,r_data,'Color',[0.85,0.325,0.098],"LineWidth",1)
hold on;
plot(model.t,r_model,'Color',[0,0.447,0.741],"LineWidth",1)
xlabel("Time (s)");
ylabel("Velocity (m/s)")
ylim([-1*max_val,max_val]);
title("Horizontal Component")