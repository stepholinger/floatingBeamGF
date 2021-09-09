L = 1e7;
f_max = 10;
t0 = 10;
t_max = 250;
h_i_avg = 400;
h_w_avg = 590;

% make model object
model = loadParameters(L,f_max,t_max,h_i_avg,h_w_avg);

[m_curve,c_ratio,m0] = moment_curve(model,t0,"basal",-0.5,0.5);
plot(m_curve/m0)
hold on
[m_curve,c_ratio,m0] = moment_curve(model,t0,"surface",0.5,-0.5);
plot(m_curve/m0,'Color',[0.85,0.325,0.098])
hydrostatic_final_ratio = (-h_i_avg/2+h_i_avg*model.rho_i/model.rho_w)/h_i_avg;
[m_curve,c_ratio,m0] = moment_curve(model,t0,"hydrostatic",-0.5,hydrostatic_final_ratio);
plot(m_curve/m0,'Color',[0.929,0.694,0.125])
yline(1,"--")
legend("Hydrostatic","m_0")
%legend("Basal","Surface","Hydrostatic","m_0")
title("Moments for Crevasse Growth Models")
ylabel("Bending Moment (Fraction m_0)")
xticks([0,length(m_curve)/4,length(m_curve)/2,3*length(m_curve)/4,length(m_curve)])
xlim([0,length(m_curve)])
xticklabels({"0%","25%","50%","75%","100%"})
xlabel("Percent Ice Thickness Fractured")
