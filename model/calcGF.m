function [dGdt_vert,dGdt_horz,stf,model] = calcGF(L,f_max,t_max,h_i,h_w,stat_dist,t0,source,mode,h_c_initial,h_c_final)

% function that provides modeled vertical seismogram (dGdt_vert) and radial velocity
% seismogram (dGdt_horz) for a floating beam with a single-force or bending moment
% source
% input parameters:
% L: length of model spatial domain
% f_max: highest desired frequency to resolve
% h_i: beam thickness
% h_w: water thickness
% stat_dist: distnace from source to receiver
% t0: pulse-width of source time function
% source: source type, either "moment" or "force"
% mode: source time function mode ("erf","gaussian","basal","surface","hydrostatic")
% h_c_initial: initial crevasse height for crevasse source time function options
% h_c_final: final crevasse height for crevasse source time function options
% h_c parameters can be left blank for non-crevasse stf options

% make model object
model = loadParameters(L,f_max,t_max,h_i,h_w,stat_dist,t0,source,mode);
if nargin > 9
    model.h_c_initial = h_c_initial;
    model.h_c_final = h_c_final;    
end

% caclulate max pressure and moment
model.P_max = model.rho_i * model.g * h_i;
model.M_max =  model.rho_i * model.g * h_i^3 / 12 * 0.072;

% get index of desired position
[~,loc_idx] = min(abs(model.x - stat_dist));

% run model
G_mat = semiAnalyticGreenFunction(model);

% take spatial derivative if moment source
if source == "moment"
    [~,G_mat] = gradient(G_mat,model.dx);
end

% attach matrix of Green's function solutions to the model
model.G_mat = G_mat;

% get vertical Green's function closest seismometer location
G_vert = G_mat(loc_idx,:);
model.G_vert = G_vert;

% get tilt Green's function
[~,tilt_mat] = gradient(model.G_mat,model.dx);
G_tilt = tilt_mat(loc_idx,:);
model.G_tilt = G_tilt;

% convolve vertical Green's function with stf
[G_vert_stf,~] = convolve_stf(G_vert,model);

% take time derivative to get velocity seismogram
dGdt_vert = gradient(G_vert_stf,model.dt);

% convolve tilt Green's function with stf
[G_tilt_stf,stf] = convolve_stf(G_tilt,model);

% get acceleration timeseries from tilt
acc_horz = G_tilt_stf * model.g;

% integrate to get horizontal (radial) velocity timeseries
dGdt_horz = cumtrapz(model.t,acc_horz);

% taper to remove static offset if using asymmetric stf (not observed by seismometers)
if mode == "erf"
    % estimate event duration to select taper window length
    cumsum_amp = cumsum(abs(G_tilt_stf));
    cumsum_amp_der = diff(cumsum_amp);
    norm_cumsum_amp_der = cumsum_amp_der/max(cumsum_amp_der);
    win_indices = find(norm_cumsum_amp_der>0.01);
    win_len = win_indices(end);
    
    % apply taper window
    taper_win = tukeywin(win_len,0.3);
    taper_win = [taper_win',zeros(1,length(dGdt_horz)-win_len)];
    dGdt_horz = taper_win.*dGdt_horz;
end

end