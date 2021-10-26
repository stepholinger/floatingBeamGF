% compares vertical and radial modeled traces and displays stf
t0 = 10;
f_max = 1;

% run model for stf 1
[z,h,stf,model] = calcGF(1e7,f_max,1000,400,590,10000,10,"moment","erf");
%[~,~,m0] = moment_curve(model,t0,"basal",-0.3,0.25);

% plot stf
subplot(2,3,1)
%scaled_stf = stf/m0;
%plot(model.t(1:t0*f_max*3),scaled_stf(end/2:end/2+t0*f_max*3-1))
plot(stf(1:model.nt))
ylabel("Bending moment (m/m0)")
title("STF")

% plot vertical
subplot(2,3,2)
plot(model.t,z)
ylabel("Velocity (m/s)")
title("Vertical")

% plot horizontal
subplot(2,3,3)
plot(model.t,h)
ylabel("Velocity (m/s)")
title("Horizontal")

% scale to max of horizontal or vertical
max_val = max(max(abs(z)),max(abs(h)));
subplot(2,3,2)
ylim([-1*max_val,max_val])
subplot(2,3,3)
ylim([-1*max_val,max_val])

% run model for stf 2
[z,h,stf,model] = calcGF(1e7,f_max,500,400,590,10000,10,"force","erf");
%[~,~,m0] = moment_curve(model,t0,"basal",0,0.25);

% plot stf
subplot(2,3,4)
%scaled_stf = stf/m0;
%plot(model.t(1:t0*f_max*3),scaled_stf(end/2:end/2+t0*f_max*3-1))
plot(stf(1:model.nt))
xlabel("Time (s)");
ylabel("Bending moment (m/m0)")

% plot vertical
subplot(2,3,5)
plot(model.t,z)
xlabel("Time (s)");
ylabel("Velocity (m/s)")

% plot horizontal
subplot(2,3,6)
plot(model.t,h)
xlabel("Time (s)");
ylabel("Velocity (m/s)")

% scale to max of horizontal or vertical
max_val = max(max(abs(z)),max(abs(h)));
subplot(2,3,5)
ylim([-1*max_val,max_val])
subplot(2,3,6)
ylim([-1*max_val,max_val])


