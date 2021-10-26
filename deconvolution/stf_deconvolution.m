% script to get a source time function by deconvolving theoretical Green's
% function from observations

% set some parameters
f_max = 1;
fs_data = 100;

% read centroid
centroidPath = "/home/setholinger/Documents/Projects/PIG/detections/templateMatch/multiTemplate/run3/short_3D_clustering/modified_k_shape/";
centroids = h5read("/home/setholinger/Documents/Projects/PIG/detections/templateMatch/multiTemplate/run3/short_3D_clustering/modified_k_shape/2/cluster_1_mean_wave_0.05-1Hz.h5","/rift_mean_wave");
z_data = squeeze(centroids)';
z_data = z_data(1:1000);

% run model for moment case
[z_model,r_model,stf,model] = calcGF(1e7,1,500,400,590,20000,5,"moment","erf");

% get the green's function that has not been convolved with an stf
G = model.G_vert;

% pad with a bunch of zeros for deconvolution and to make model and data
% lengths the same
G_pad = zeros(1,10000);
G_pad(ceil(end/2):ceil(end/2)+length(G)-1) = G;
z_data_pad = zeros(1,10000);
z_data_pad(ceil(end/2):ceil(end/2)+length(z_data)-1) = z_data;

% deconvolve the Green's function from the observation
stf = ifft(fft(z_data_pad)./fft(G_pad));

% test the stf
z_data_test = ifft(fft(stf).*fft(G_pad));

% test a lowpass filtered verstion of the stf
% choose filter band and design filter
freq = 0.1;
[b,a] = butter(4,freq/(2/2),'low');
stf_filt = filtfilt(b,a,stf);
z_data_test = ifft(fft(stf_filt).*fft(G_pad));
