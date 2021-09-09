% set parameters
L = 1e7;
f_max = 2; 
t_max = 5000;
h_i = 100;
h_w = 590;
t0 = 8.75;
h_c_initial = 0;
h_c_final = 0.25;
mode = 'basal';
source = 'moment';

% setup below shows artifacts from stf corner
%t0 = 100;
%h_c_initial = -0.15;
%h_c_final = 0.25;

% make model object
model = loadParameters(L,f_max,t_max,h_i,h_w);

% caclulate max pressure and moment
P_max = 916 * 9.8 * h_i;
M_max =  916 * 9.8 * h_i^3 / 12 * 0.072;

% set length of animated portion of beam in km
view_length = 100;
view_idx = floor(view_length*1000/model.dx);

% run model
G = semiAnalyticGreenFunction(model);

% take spatial derivative if moment source
if source == "moment"
    [~,G] = gradient(G,model.dx);
end

% make array for waves after convolution with source time function
G_stf = zeros(view_idx,model.nt);

% get source time function and convolve with Green's function
if sum(contains(["basal","surface","hydrostatic"],mode)) > 0 
     
    % call stf function with crevasse mode
    stf = source_time_function(model,t0,mode,h_c_initial,h_c_final);
    
    % convolve stf with Green's function
    for x = 1:view_idx
        G_pad = zeros(size(model.t)*2-1);
        G_pad(ceil(end/2):end) = G(end/2+x,:);
        G_stf_pad = ifft(fft(G_pad).*fft(stf));
        G_stf(x,:) = G_stf_pad(1:ceil(end/2));
    end

elseif mode == "erf"
    
    % scale Green's function output
    if source == "moment"
        G = G * M_max;
    elseif source == "force"
        G = G * P_max;
    end   
    
    % call stf function with erf mode
    stf = source_time_function(model,t0,mode);
    
    % convolve stf with Green's function
    for x = 1:view_idx
        G_pad = zeros(size(model.t)*2-1);
        G_pad(ceil(end/2):end) = G(end/2+x,:);
        G_stf_pad = ifft(fft(G_pad).*fft(stf));
        G_stf(x,:) = G_stf_pad(1:ceil(end/2));
        %G_stf(x,:) = G_stf_pad(ceil(end/2):end);
    end
    
elseif mode == "gaussian"
    
    % scale Green's function output
    if source == "moment"
        G = G * M_max;
    elseif source == "force"
        G = G * P_max;
    end   
    
    % call stf function with gaussian mode
    stf = source_time_function(model,t0,mode);
    
    for x = 1:view_idx
        G_pad = zeros(size(model.t)*2-1);
        G_pad(ceil(end/2):end) = G(end/2+x,:);
        G_stf_pad = ifft(fft(G_pad).*fft(stf));
        G_stf(x,:) = G_stf_pad(1:ceil(end/2));
    end
end

% do plot setup
ymax = max(G_stf,[],'all');
ymin = min(G_stf,[],'all');

ax = axes();
for t = 1:model.nt
    plot(ax,[1:view_idx]*model.dx/1000,G_stf(:,t))
    ylim(ax,[ymin,ymax])
    xlim(ax,[0,view_idx*model.dx/1000])
    xlabel("Position (km)")
    ylabel("Deflection amplitude (m)")
    title("Ice shelf bending forced with " + mode + " " + source + " source")
    drawnow    
end