function data = read_centroid_mcmc(centroid_path,data_set)

% read h5 file
centroids = h5read(centroid_path,data_set);
centroids = squeeze(centroids)';

% get HHZ only (first 3rd of trace)
data = centroids(1:floor(end/3));

end