% set model parameters
L = 1e7;
f_max = 2;
t_max = 1000;
h_i = 400;
h_w = 590;
source_dist_initial = 20000;
source_dist_final = 15000;
rupture_vel = 20;
t0 = 200;
source = "moment";
mode = "basal";
h_c_initial = 0;
h_c_final = 0.5;

% set length of viewing window in seconds 
viewLen = 7200;

% run model to get dx
[model,dGdt,G,stf] = calcGF(L,f_max,t_max,h_i,h_w,source_dist_initial,t0,source,mode,h_c_initial,h_c_final);

% set up vector for superposed waves
max_t_offset = round((abs(source_dist_final-source_dist_initial))/rupture_vel*(1/model.dt));
dGdt_sup = zeros(1,max_t_offset+length(dGdt));  
    
% superpose solutions along rupture length
for d = source_dist_initial:-model.dx:source_dist_final

    % run model
    [model,dGdt,G,stf] = calcGF(L,f_max,t_max,h_i,h_w,d,t0,source,mode,h_c_initial,h_c_final);
    
    % get time offset
    t_offset = round((abs(d-source_dist_initial))/rupture_vel*(1/model.dt));

    % offset dGdt to account for time between radiation from each point
    %if t_offset > 0
    dGdt_offset = [zeros(1,t_offset),dGdt,zeros(1,max_t_offset-t_offset)];
    %else
    %    dGdt_offset = [dGdt,zeros(1max_t_offset-t_offset)];
    %end
    
    dGdt_sup = dGdt_sup + dGdt_offset;
end

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
%freq = 1/(t0/2);
%[b,a] = butter(4,freq/(fsNew/2),'low');
%trace = filtfilt(b,a,trace);

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

% remove scalar offset using first value
eventTrace = eventTrace - eventTrace(1);

% cross correlate with data
dGdt_sup = dGdt_sup(1:10*length(eventTrace));
[xcorr_trace,lag] = xcorr([eventTrace,zeros(1,9*length(eventTrace))],dGdt_sup,"coeff");
[xcorr_coef,lag_idx] = max(xcorr_trace);  

slide = lag(lag_idx);

if slide > 0
    dGdt_sup = [zeros(1,abs(slide)),dGdt_sup];
    dGdt_sup = dGdt_sup(1:length(eventTrace));
end
%if slide < 0
%    dGdt_sup = dGdt_sup(abs(slide):end);
%    if length(dGdt_sup) < length(eventTrace)
%        dGdt_sup = [dGdt_sup,zeros(1,length(eventTrace)-length(dGdt_sup))];
%    elseif length(dGdt_sup) > length(eventTrace)
%        dGdt_sup = dGdt_sup(1:length(eventTrace));
%    end
%end

figure(1)
t_vect = model.dt:model.dt:length(eventTrace)*model.dt;
plot(t_vect,dGdt_sup(1:length(eventTrace)))
hold on;
plot(t_vect,eventTrace)
