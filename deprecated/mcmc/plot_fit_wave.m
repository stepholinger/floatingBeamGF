function plot_fit_wave(t,eventAlign,sigma,L_fit,M_frac_fit,G_fit,xFit,numIt,xStep,p,accept,L_type,path)

% plot resulting waveform and starting waveform
figure('Position', [0 1 1000 1000])
hold on;
plot(t,G_fit,'LineWidth',2.5,'Color',[0.8500 0.3250 0.0980]);
plot(t,eventAlign,'LineWidth',2.5,'Color',[0 0.4470 0.7410]);
plot(t,eventAlign+std(eventAlign),'Color',[0 0.4470 0.7410],'LineStyle','--')
plot(t,eventAlign-std(eventAlign),'Color',[0 0.4470 0.7410],'LineStyle','--')
ylim([-1.5,1.5])
title("Cluster " + string(p-1) + " Centroid and Best-Fit Waveform")
xlabel("Time (s)")
ylabel("Normalized Velocity")
[k,mk] = max(G_fit);            
% text(mk*4/3,k/1.5,string("Best-fit parameters" + newline + "----------------------------" + ...
%                         newline + "h_i: " + xFit(1) + " m     h_w: " ...
%                         + xFit(2) + " m" + newline + "X_{stat}: " + ...
%                         xFit(3)/1000 + " km" + "     t_0: " + xFit(4) + ...
%                         newline + "M/M_0: " + M_frac_fit + newline))
% text(mk*4/3,k/3,string("MCMC parameters" + newline + "----------------------------" + newline + ...                            
%                        "h_i step: " + xStep(1) + " m    h_w step: " + xStep(2) + " m" + newline + ...
%                        "X_{stat} step: " + xStep(3) + " m    t_0 step: " + xStep(4) + " s" + newline + ...
%                        "Number of iterations: " + numIt + newline + "Sigma: " + sigma + newline + ...
%                        "Liklihood function: " + L_type + newline + "L: " + L_fit + newline + "Accepted " + ...
%                        round(100*sum(accept)/length(accept)) + "% of proposals"))                   
l = legend("MCMC best-fit model","Centroid","1\sigma");
set(l,'Location','southwest');
set(gcf,'Position',[10 10 1000 800])
set(findall(gcf,'-property','FontSize'),'FontSize',16)
hold off
saveas(gcf,path + "centroid" + string(p-1) + "_fit_wave.png")
%close(gcf)

end