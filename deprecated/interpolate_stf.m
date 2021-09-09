function [spline] = interpolate_stf(model,m_curve,c_ratio,mode,h_c_initial,h_c_final,interp_len)

% RIGHT NOW THIS IS JUST THE FIRST SPLINE

% set parameters
rho_i = model.rho_i;
rho_w = model.rho_w;
fs = 1/model.dt;
g = model.g;
h_i = model.h_i;
h_w = -h_i/2 + h_i * rho_i/rho_w;


if mode == "basal"
    
    % get points for data vector
    m_initial = m_curve(1);
    m_final = m_curve(end);
    
    h_c = (h_c_initial+interp_len)*h_i;

    [~,m1_idx] = min(abs(c_ratio - (h_c_initial+interp_len)));
    m1 = m_curve(m1_idx);
    
    % write out derivative expressions
    if h_c < h_w
        m1_1st_der = -(h_c^2)*g*rho_i - (1/2)*h_c*g*h_i*rho_i + (h_c^2)*g*rho_w + (1/2)*h_c*g*h_i*rho_w;
        m1_2nd_der = -2*h_c*g*rho_i - (1/2)*g*h_i*rho_i + 2*h_c*g*rho_w + (1/2)*g*h_i*rho_w;
    else
        m1_1st_der = -(h_c^2)*g*rho_i + (1/2)*h_c*g*h_i*rho_i;
        m1_2nd_der = -2*h_c*g*rho_i + (1/2)*g*h_i*rho_i;
    end

    h_c = (h_c_final-interp_len)*h_i;
    
    [~,m2_idx] = min(abs(c_ratio - (h_c_final-interp_len)));
    m2 = m_curve(m2_idx);
    
    % write out derivative expressions
    if h_c < h_w
        m2_1st_der = -(h_c^2)*g*rho_i - (1/2)*h_c*g*h_i*rho_i + (h_c^2)*g*rho_w + (1/2)*h_c*g*h_i*rho_w;
        m2_2nd_der = -2*h_c*g*rho_i - (1/2)*g*h_i*rho_i + 2*h_c*g*rho_w + (1/2)*g*h_i*rho_w;
    else
        m2_1st_der = -(h_c^2)*g*rho_i + (1/2)*h_c*g*h_i*rho_i;
        m2_2nd_der = -2*h_c*g*rho_i + (1/2)*g*h_i*rho_i;
    end
    
    % ENDPOINTS. 1ST, AND 2ND DERIVATIVES
    % fill data vector for inverting for coefficients
    %d = [m1,0,0,m2,0,0];
    
    % make model matrix of cubic spline basis functions
    %m = [[cubic_basis(0,1,0),cubic_basis(0,2,0),cubic_basis(0,3,0),cubic_basis(0,4,0)];
    %         [cubic_basis(0,1,1),cubic_basis(0,2,1),cubic_basis(0,3,1),cubic_basis(0,4,1)];
    %         [cubic_basis(0,1,2),cubic_basis(0,2,2),cubic_basis(0,3,2),cubic_basis(0,4,2)];
    %         [cubic_basis(1,1,0),cubic_basis(1,2,0),cubic_basis(1,3,0),cubic_basis(1,4,0)];
    %         [cubic_basis(1,1,1),cubic_basis(1,2,1),cubic_basis(1,3,1),cubic_basis(1,4,1)];
    %         [cubic_basis(1,1,2),cubic_basis(1,2,2),cubic_basis(1,3,2),cubic_basis(1,4,2)]];
     
    % ENDPOINTS AND 2ND DERIVATIVES
    % fill data vector for inverting for coefficients
    d = [m_initial, 0, m1, m1_2nd_der, m2, m2_2nd_der, m_final, 0];
   
    % make model matrix of cubic spline basis functions
    m = [[cubic_basis(0,1,0),cubic_basis(0,2,0),cubic_basis(0,3,0),cubic_basis(0,4,0)];
             [cubic_basis(0,1,2),cubic_basis(0,2,2),cubic_basis(0,3,2),cubic_basis(0,4,2)];
             [cubic_basis(interp_len,1,0),cubic_basis(interp_len,2,0),cubic_basis(interp_len,3,0),cubic_basis(interp_len,4,0)];
             [cubic_basis(interp_len,1,2),cubic_basis(interp_len,2,2),cubic_basis(interp_len,3,2),cubic_basis(interp_len,4,2)];
             [cubic_basis(1-interp_len,1,0),cubic_basis(1-interp_len,2,0),cubic_basis(1-interp_len,3,0),cubic_basis(1-interp_len,4,0)];
             [cubic_basis(1-interp_len,1,2),cubic_basis(1-interp_len,2,2),cubic_basis(1-interp_len,3,2),cubic_basis(1-interp_len,4,2)];
             [cubic_basis(1,1,0),cubic_basis(1,2,0),cubic_basis(1,3,0),cubic_basis(1,4,0)];
             [cubic_basis(1,1,2),cubic_basis(1,2,2),cubic_basis(1,3,2),cubic_basis(1,4,2)]];
     
    % ENDPOINTS AND 1ST DERIVATIVES
    % fill data vector for inverting for coefficients
    %d = [m1,0,m2,0];
    
    % make model matrix of cubic spline basis functions
    %m = [[cubic_basis(0,1,0),cubic_basis(0,2,0),cubic_basis(0,3,0),cubic_basis(0,4,0)];
    %         [cubic_basis(0,1,1),cubic_basis(0,2,1),cubic_basis(0,3,1),cubic_basis(0,4,1)];
    %         [cubic_basis(1,1,0),cubic_basis(1,2,0),cubic_basis(1,3,0),cubic_basis(1,4,0)];
    %         [cubic_basis(1,1,1),cubic_basis(1,2,1),cubic_basis(1,3,1),cubic_basis(1,4,1)]];     
          
    % ENDPOINTS, 1ST DERIVATIVE FOR POINT 1, 2ND DERIVATIVE FOR POINT 2
    % fill data vector for inverting for coefficients
    %d = [m1,0,m2,d2m_dx2];
        
    % make model matrix of cubic spline basis functions
    %m = [[cubic_basis(0,1,0),cubic_basis(0,2,0),cubic_basis(0,3,0),cubic_basis(0,4,0)];
    %         [cubic_basis(0,1,1),cubic_basis(0,2,1),cubic_basis(0,3,1),cubic_basis(0,4,1)];
    %         [cubic_basis(1,1,0),cubic_basis(1,2,0),cubic_basis(1,3,0),cubic_basis(1,4,0)];
    %         [cubic_basis(1,1,2),cubic_basis(1,2,2),cubic_basis(1,3,2),cubic_basis(1,4,2)]];     
             
    % do the inversion
    c = m\d';
         
    % plot the resulting spline
    n_pts = length(m_curve);
    spline = zeros(n_pts-1,1);
    for s = 1:n_pts-1
        spline(s) = c(1)*cubic_basis(s/n_pts,1,0)+c(2)*cubic_basis(s/n_pts,2,0)+c(3)*cubic_basis(s/n_pts,3,0)+c(4)*cubic_basis(s/n_pts,4,0);
    end
    % make plot to check that slope is correct
    %plot(m_curve); hold on
    %plot(m_curve(1) - dm_dx*h_c + dm_dx*(h_c_initial*h_i:(h_c_final-h_c_initial)*h_i/length(m_curve):h_c_final*h_i))

elseif mode == "surface"


elseif mode == "hydrostatic"


end


end