function stf = source_time_function(model,t0,mode,h_c_initial,h_c_final)

% if enough arguments were supplied, call the crevasse moment function
if nargin > 3
    
    % call moment_curve function and get portion of moment curve corresponding to h_c_initial and h_c_final
    [m_curve,c_ratio] = moment_curve(model,t0,mode,h_c_initial,h_c_final);
    [~,h_c_initial_idx] = min(abs(c_ratio - h_c_initial));
    [~,h_c_final_idx] = min(abs(c_ratio - h_c_final));
    m_curve = m_curve(h_c_initial_idx:h_c_final_idx);
    m_curve_len = length(m_curve);
    delta_m_curve = m_curve - m_curve(1);
    c_ratio = c_ratio(h_c_initial_idx:h_c_final_idx);
    
    % pad m curve on either side with 0s and 1s
    pad_len = 10000;
    delta_m_curve_pad = [zeros(1,pad_len),delta_m_curve,ones(1,pad_len)*delta_m_curve(end)];
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   
    % interpolation approaches 
    
    % interpolate with cubic spline (endpoints conditions specified)
    %interp_points = [1:200,1020:1300,2200:2400];
    %interp_vals = delta_m_curve_pad(interp_points);
    %stf_func = csape(interp_points,interp_vals,'clamped',[0,0]);
    %stf = fnval(stf_func,1:length(delta_m_curve_pad));
    
    % interpolate with smoothing cubic spline with even weighting over
    % entire delta_m_curve_pad
    stf_func = csaps(1:length(delta_m_curve_pad),delta_m_curve_pad,1e-8);
    stf = fnval(stf_func,1:length(delta_m_curve_pad));
    
    % interpolate with smoothing cubic spline with even weighting over
    % specific regions of delta_m_curve_pad
    %interp_points = [1:pad_len-1000,pad_len+100:pad_len+length(delta_m_curve)-100,pad_len+length(delta_m_curve)+1000:pad_len*2+length(delta_m_curve)];
    interp_points = sort([1:100:length(delta_m_curve_pad),find(islocalmax(abs(delta_m_curve_pad)))]);
    interp_vals = delta_m_curve_pad(interp_points);
    stf_func = csaps(interp_points,interp_vals,1e-4);
    stf = fnval(stf_func,1:length(delta_m_curve_pad));
    
    % interpolate with smoothing cubic spline with weighting vector over
    % entire delta_m_curve_pad
    %weight_vect = [100000000*ones(1,100),ones(1,pad_len-100),10000*ones(1,length(delta_m_curve)),ones(1,pad_len-100),100000000*ones(1,100)];
    %stf_func = csaps(1:length(delta_m_curve_pad),delta_m_curve_pad,1e-8,[],weight_vect);
    %stf = fnval(stf_func,1:length(delta_m_curve_pad));
    
    % interpolate with smoothing spline
    %stf_func = spaps(1:100:length(delta_m_curve_pad),delta_m_curve_pad(1:100:length(delta_m_curve_pad)),1e9);
    %stf = fnval(stf_func,1:length(delta_m_curve_pad));
   
    
    % get interploation points- two endpoints and 1 midpoint of m_curve region
    %curve_len = length(stf);
    %interp_points = [1,2,floor(curve_len/2),curve_len-1,curve_len];
    %interp_vals = [stf(1),stf(1),stf(floor(curve_len/2)),stf(end),stf(end)];
 
    % get interpolation points by sampling input
    %[interp_vals,interp_points] = datasample(stf,ramp_t_len);
    %interp_points = 1:10000:length(stf);
    %interp_points = [1,floor(m_curve_len/4),floor(m_curve_len/2),floor(3*m_curve_len/4),m_curve_len,floor(curve_len/2),curve_len];
    %interp_vals = stf(interp_points);

    % interpolate entire stf with spline to avoid corners
    %stf = csape(interp_points,interp_vals,'periodic',[0,0]);
    %stf = fnval(stf,1:1:curve_len);
    
    % pchip interpolation
    interp_points = unique(sort([1:100:length(delta_m_curve_pad),find(islocalmax(abs(delta_m_curve_pad)))]));
    interp_vals = delta_m_curve_pad(interp_points);
    stf = pchip(interp_points,interp_vals,1:length(delta_m_curve_pad));
    
    % add erf ramp to make stf periodic
    ramp_len = 1000000;
    ramp_t = [-ramp_len/2:model.dt:ramp_len/2];
    ramp = (erf(ramp_t/10000)+1)/2;
    [~,start_index] = max(find(ramp == 0));
    [end_index,~] = min(find(ramp == 1));
    ramp = ramp*stf(end);
    stf = [stf,fliplr(ramp(start_index:end_index))];
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %smoothing approaches
    ramp_len = 1000000;
    ramp_t = [-ramp_len/2:model.dt:ramp_len/2];
    ramp = (erf(ramp_t/10000)+1)/2;
    [~,start_index] = max(find(ramp == 0));
    [end_index,~] = min(find(ramp == 1));
    ramp = ramp*delta_m_curve(end);
    pad_len = 2*model.nt-length(delta_m_curve);
    stf = [delta_m_curve_pad,ones(1,pad_len-1)*delta_m_curve(end),fliplr(ramp(start_index:end_index))];
    
    % smooth stf with moving average
    %stf = smooth(stf,10000)';
    
    % smooth stf with savitzky-Golay
    %stf = sgolayfilt(stf,10,10001);
   
    % manual gaussian filter
    new_t = [0:model.dt:length(stf)/4];
    %new_t = [0:1:length(stf)];
    gauss = exp((((new_t-(model.dt*length(new_t))/2)/t0).^2)/-2)/t0/sqrt(2*pi);
    %gauss = gauss./max(gauss);
    gauss_stf = ifft(fft(gauss(1:end-1)).*fft(stf));

    % matlab gaussian filter (probably does better job of preserving
    % amplitude
    stf = smoothdata(stf,'gaussian',t0*1/model.dt/2);

    
    
    
    % weird erf concatenation idea
    %new_t = [-100000/2:model.dt:100000/2];
    %erf_stf = (erf(new_t/(t0/4))+1)/2;
    %[~,start_index] = max(find(erf_stf == 0));
    %[end_index,~] = min(find(erf_stf == 1));
    %erf_stf = erf_stf(start_index:end_index);

    %erf_stf = [min(delta_m_curve)*erf_stf(1:800),(min(delta_m_curve)-delta_m_curve(end))*fliplr(erf_stf(1:800))+delta_m_curve(end)];
    %erf_stf = [min(delta_m_curve)*erf_stf,(min(delta_m_curve)-delta_m_curve(end))*fliplr(erf_stf)+delta_m_curve(end)];

    % pad stf to correct length for convolution and add smooth ramp to make stf periodic
    %ramp_len = (length(erf_stf)-(end_index-start_index));
    %ramp = interp1q([0,ramp_len]',[1,0]',[1:ramp_len]');
    %stf = [erf_stf(start_index:end_index),ramp'];

    % add erf ramp to make stf periodic
    %ramp_len = 1000000;
    %ramp_t = [-ramp_len/2:model.dt:ramp_len/2];
    %ramp = (erf(ramp_t/10000)+1)/2;
    %[~,start_index] = max(find(ramp == 0));
    %[end_index,~] = min(find(ramp == 1));
    %ramp = ramp*erf_stf(end);
    %stf = [erf_stf,fliplr(ramp(start_index:end_index))];


    
    
    
% if 3 or less arguments were supplied, use the simple stfs
else
      
    % internal parameter sets length of vector used to construct special
    % functions erf and gaussian pulse
    new_t_len = 100000;
    
    % if half pulse source type, use error function to make stf
    if mode == "erf"
        % make error function pulse
        try
            new_t = [-new_t_len/2:model.dt:new_t_len/2];
            erf_stf = (erf(new_t/t0)+1)/2;
            [~,start_index] = max(find(erf_stf == 0));
            [end_index,~] = min(find(erf_stf == 1));
            erf_stf = erf_stf(start_index:end_index);
            
            % pad stf to correct length for convolution and add smooth ramp to make stf periodic
            %ramp_len = (length(erf_stf)-(end_index-start_index));
            %ramp = interp1q([0,ramp_len]',[1,0]',[1:ramp_len]');
            %stf = [erf_stf(start_index:end_index),ramp'];
            
            % add erf ramp to make stf periodic
            ramp_len = 1000000;
            ramp_t = [-ramp_len/2:model.dt:ramp_len/2];
            ramp = (erf(ramp_t/10000)+1)/2;
            [~,start_index] = max(find(ramp == 0));
            [end_index,~] = min(find(ramp == 1));
            stf = [erf_stf,fliplr(ramp(start_index:end_index))];
            
            
        if isempty(start_index) || isempty(end_index)
            error("Use a larger t_max for a pulse with t0 = " + string(t0))
        end
        catch
            error("Use a larger t_max for a pulse with t0 = " + string(t0))
        end
        
    elseif mode == "gaussian"   
        % make gaussian pulse      
        try
            new_t = [0:model.dt:new_t_len];
            gauss_stf = exp((((new_t-(model.dt*length(new_t))/2)/t0).^2)/-2)/t0/sqrt(2*pi);
            gauss_stf = gauss_stf./max(gauss_stf);
            [~,offset_index] = max(find(gauss_stf(1:ceil(end/2)) == 0));
            stf = [gauss_stf(offset_index:end),zeros(1,offset_index)];
            
        if isempty(offset_index)
            error("Use a longer new_t vector for gaussian pulse with t0 = " + string(t0))
        end
        catch
            error("Use a longer new_t vector for gaussian pulse with t0 = " + string(t0))
        end
    end
end