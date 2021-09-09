function [G_stf,stf] = convolve_stf(G,model)

% convolve Green's function with stf
if sum(contains(["basal","surface","hydrostatic"],model.mode)) > 0 
    
    % call stf function with crevasse mode
    stf = source_time_function(model,model.t0,model.mode,model.h_c_initial,model.h_c_final);
    
    % convolve stf Green's function
    G_pad = zeros(1,length(stf));
    G_pad(ceil(end/2):ceil(end/2)+length(G)-1) = G;
    G_stf_pad = ifft(fft(G_pad).*fft(stf));
    
    % get rid of padding
    G_stf = G_stf_pad;
    G_stf = G_stf(1:model.nt);
    
elseif model.mode == "erf"
    
    % scale Green's function output
    if model.source == "moment"
        G = G * model.M_max;
    elseif model.source == "force"
        G = G * model.P_max;
    end    
    
    % call stf function with erf mode
    stf = source_time_function(model,model.t0,model.mode);
    
    % convolve stf with Green's function
    G_pad = zeros(1,length(stf));
    G_pad(ceil(end/2):ceil(end/2)+length(G)-1) = G;
    G_stf_pad = ifft(fft(G_pad).*fft(stf));
    
    % remove padding
    G_stf = G_stf_pad(ceil(end/2):end);
    G_stf = G_stf(1:model.nt);

elseif model.mode == "gaussian"
    
    % scale Green's function output
    if model.source == "moment"
        G = G * model.M_max;
    elseif model.source == "force"
        G = G * model.P_max;
    end    
    
    % call stf function with gaussian mode
    stf = source_time_function(model,model.t0,model.mode);
    
    % convolve stf with Green's function
    G_pad = zeros(1,length(stf));
    G_pad(ceil(end/2):ceil(end/2)+length(G)-1) = G;
    G_stf_pad = ifft(fft(G_pad).*fft(stf));
    
    % get rid of padding
    G_stf = G_stf_pad(ceil(end/2):end);
    G_stf = G_stf(1:model.nt);
    
end

end