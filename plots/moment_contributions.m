% set model input parameters and make model
L = 1e7;
f_max = 2;
t_max = 250;
h_i = 400;
h_w = 590;
mode = "basal";
h_c_inital = -0.5;
h_c_final = 0.5;
model = loadParameters(L,f_max,t_max,h_i,h_w);

% set parameters
rho_i = model.rho_i;
rho_w = model.rho_w;
h_i = model.h_i;
h_w = -h_i/2 + h_i * rho_i/rho_w;
g = 9.8;
fs = 1/model.dt;

% calculate c_ratio spacing to correspond to desired t0 and h_c range:
% we want the portion of the curve between h_c_initial and h_c_final to
% be t0*fs long
m_overburden = -(1/12)*g*rho_i*h_i^3;
m_waterline_buoyancy = (1/6)*(1/(rho_w^2))*g*(rho_i^3)*(h_i^3) - (1/4)*(1/(rho_w))*g*(rho_i^2)*(h_i^3);
m0 = m_overburden - m_waterline_buoyancy;

if mode == "basal"
    
    % make vector of crevasse height ratios from -0.5 to 0.5
    c_ratio_spacing = (h_c_final-h_c_initial)/(t0*fs);
    c_ratio = -0.5:c_ratio_spacing:0.5;

    % make storage vector
    moments = zeros(length(c_ratio),4);

    % calculate moment for each crevasse fraction
    for c = 1:length(c_ratio)

        h_bc = c_ratio(c)*h_i;
        
        % calculate first part of moment expression
        m_overburden = -(1/12)*g*rho_i*h_i^3;
        m_crevasse_overburden = (1/3)*g*(rho_i)*(h_bc^3) - (1/4)*g*(rho_i)*h_i*(h_bc^2) + (1/48)*g*(rho_i)*(h_i^3);
        m_waterline_buoyancy = (1/6)*(1/(rho_w^2))*g*(rho_i^3)*(h_i^3) - (1/4)*(1/(rho_w))*g*(rho_i^2)*(h_i^3);
        m_crevasse_buoyancy = (1/2)*g*(h_bc^2)*(h_i)*(rho_i) - (1/8)*g*(h_i^3)*(rho_i) - (1/3)*g*(h_bc^3)*(rho_w) - ...
                       (1/4)*g*(h_bc^2)*(h_i)*(rho_w) + (1/48)*g*(h_i^3)*(rho_w);
        if h_bc < h_w
            m = m_overburden - m_crevasse_buoyancy - m_crevasse_overburden;
            moments(c,3) = m_crevasse_buoyancy;
        else
            m = m_overburden - m_waterline_buoyancy - m_crevasse_overburden;
            moments(c,3) = m_waterline_buoyancy;
        end
        moments(c,1) = m;
        moments(c,2) = m_overburden;
        moments(c,4) = m_crevasse_overburden;
    end
    
    m_curve = moments';
end
