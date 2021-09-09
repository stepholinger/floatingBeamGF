% set model parameters
L = 1e6;
f_max = 2;
t_max = 1000;
h_i = 400;
h_w = 590;
statDist = 21500;
t0 = 5;
source = "moment";
mode = "hydrostatic";

figure(2)
hold on;
h_c_limit = -0.5 + 0.8945;
cols = ["red","blue","green","black"];
ind = 1;
for h_c = -0.5:0.25:(h_c_limit-0.1)
    
    h_c_initial = h_c;
    h_c_final = h_c + 0.1;
    fprintf("Initial: " + num2str(h_c_initial) + "     Final: " + num2str(h_c_final));
    
    % run model
    [model,dGdt,G,stf] = calcGF(L,f_max,t_max,h_i,h_w,statDist,t0,source,mode,h_c_initial,h_c_final);
    
    plot(dGdt,'color',cols(ind))
    ind = ind + 1;
end