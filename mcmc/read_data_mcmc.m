function data = read_data_mcmc(data_file,fs,f_max,event_bounds)

% read mseed file
dataStruct = rdmseed(data_file);

% extract trace
trace = extractfield(dataStruct,'d');

% resample data to 1 Hz
fs_new = f_max*2;
if fs_new < 1
    fs = fs/(f_max*2);
    trace = resample(trace,1,fs);
else
    trace = resample(trace,fs_new,fs);
end

% get event bound indices from start and end times
start_time = event_bounds(1)*fs_new;
end_time = event_bounds(2)*fs_new;

% trim data to event bounds
event_trace = trace(start_time:end_time-1);

% remove scalar offset using first value
data = event_trace - event_trace(1);

end
