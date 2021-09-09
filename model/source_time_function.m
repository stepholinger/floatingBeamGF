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

    % make erf ramp to ensure periodicity of stf 
    ramp_len = 1000000;
    ramp_t = [-ramp_len/2:model.dt:ramp_len/2];
    ramp = (erf(ramp_t/10000)+1)/2;
    [~,start_index] = max(find(ramp == 0));
    [end_index,~] = min(find(ramp == 1));
    ramp = ramp*delta_m_curve(end);
    pad_len = 2*model.nt-length(delta_m_curve);
    stf = [zeros(1,length(ramp(start_index:end_index))+length(delta_m_curve)),delta_m_curve,fliplr(ramp(start_index:end_index))];

    % apply gaussian smoothing filter to smooth corners in stf
    stf = smoothdata(stf,'gaussian',t0*1/model.dt/2);
    
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