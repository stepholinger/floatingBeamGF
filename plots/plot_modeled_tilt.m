% script to get tilt timeseries approximation from model by taking extra
% spatial derivative (dG/dx) where G is the displacement field generated

% set some variables
stat_dist = 20000;

% solve Green's function
[model,dGdt,G,stf] = calcGF(1e7,1,1000,400,590,stat_dist,10,"moment","basal",-0.5,0.5);

% take spatial derivative of displacements to get tilt
[~,tilt_mat] = gradient(model.G_mat,model.dx);

% get index of desired position
[~,loc_idx] = min(abs(model.x - stat_dist));

% get tilt timeseries at desired position
tilt = tilt_mat(loc_idx,:);

% make plot
plot(model.t,tilt)
ylabel("Tilt (m/m)")
xlabel("Time (s)")

% get acceleration using tilt relation
acceleration = model.g * tilt;

% integrate to get velocity
horizontal_velocity = cumtrapz(model.t,acceleration);

% make plot of velocity
plot(model.t,horizontal_velocity)
ylabel("Horizontal velocity (m/s)")
xlabel("Time (s)")