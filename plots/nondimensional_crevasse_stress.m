% this is a function to plot the nondimensional stress field up the 
% crevassed ice column as a function of nondimensional ice thickness 
% and crevasse height

% set parameters
resolution = 0.001;
mode = "surface";

% make model object for moment curve calculation
L = 1e7;
f_max = 2;
t0 = 10;
t_max = 1000;
h_i_avg = 400;
h_w_avg = 590;

% make model object
model = loadParameters(L,f_max,t_max,h_i_avg,h_w_avg);
rho_i = model.rho_i;
rho_w = model.rho_w;

% calculate water height h_w
h_w = -1/2 + rho_i/rho_w;

% define the two stress functions that we will call based on vertical
% position and crevasse model
opposed_overburden = @(y) (-1/2+rho_w/(2*rho_i)+y*(rho_w-rho_i)/rho_i)*(1 - rho_i/rho_w)^(-1);
unopposed_overburden = @(y) (1/2 - y)*(1 - rho_i/rho_w)^(-1);

% iterate through nondimensional position up ice column (y/h_i)
y_vect = -0.5:resolution:0.5;
if mode == "hydrostatic"
    h_c_vect = -0.5:resolution:h_w;
else
    h_c_vect = -0.5:resolution:0.5;
end
sigma_net = zeros(length(y_vect),length(h_c_vect));
for y = 1:length(y_vect)
    
    % iterate through nondimensional crevasse height (h_bc/h_i)
    if mode == "basal"

        for h_bc = 1:length(h_c_vect)
            if h_c_vect(h_bc) < h_w 
                if y_vect(y) < h_c_vect(h_bc)
                    sigma_net(y,h_bc) = opposed_overburden(y_vect(y));
                else
                    sigma_net(y,h_bc) = 0;
                end
            else
                if y_vect(y) < h_w 
                    sigma_net(y,h_bc) = opposed_overburden(y_vect(y));
                elseif y_vect(y) > h_w && y_vect(y) < h_c_vect(h_bc)
                    sigma_net(y,h_bc) = unopposed_overburden(y_vect(y));
                elseif y_vect(y) > h_c_vect(h_bc)
                    sigma_net(y,h_bc) = 0;
                end
            end    
        end
        
     elseif mode == "surface"
             
         for h_sc = 1:length(h_c_vect)   
             if h_c_vect(h_sc) > h_w
                if y_vect(y) > h_c_vect(h_sc)
                    sigma_net(y,h_sc) = unopposed_overburden(y_vect(y));
                else
                    sigma_net(y,h_sc) = 0;
                end
             else
                if y_vect(y) > h_w
                    sigma_net(y,h_sc) = unopposed_overburden(y_vect(y));
                elseif y_vect(y) < h_w && y_vect(y) > h_c_vect(h_sc)
                    sigma_net(y,h_sc) = opposed_overburden(y_vect(y));
                elseif y_vect(y) < h_c_vect(h_sc)
                    sigma_net(y,h_sc) = 0;
                end 
             end
         end
         
    elseif mode == "hydrostatic"

         for h_bc = 1:length(h_c_vect)   
            h_sc =  (1/2 - h_w)*(h_w - h_c_vect(h_bc))/(h_w + 1/2) + h_w;
            if y_vect(y) > h_sc
                sigma_net(y,h_bc) = unopposed_overburden(y_vect(y));
            elseif y_vect(y) > h_c_vect(h_bc) && y_vect(y) < h_sc
                sigma_net(y,h_bc) = 0;
            elseif y_vect(y) < h_c_vect(h_bc)
                sigma_net(y,h_bc) = opposed_overburden(y_vect(y));
            end
             
         end 
    end
end

% surface/contour plot
fig1 = figure(1);
imagesc(h_c_vect,y_vect,sigma_net)
set(gca,'YDir','normal')
xticks([-0.5,0,h_w,0.5])
xticklabels(["-h_i/2","0","h_w","h_i/2"])
xlabel("Crevasse height (" + mode + ")")
yticks([-0.5,0,h_w,0.5])
yticklabels(["-h_i/2","0","h_w","h_i/2"])
ylabel("Position along ice column")
title("Nondimensional net stress along ice column")
h = colorbar;
ylabel(h,"$\sigma_{\Sigma}^\prime$","Interpreter","latex","fontsize",20)
ylim(h,[0,1])

% get moment curve
if mode == "surface"
    [m_curve,c_ratio,m0] = moment_curve(model,1/resolution/4,mode,0.5,-0.5);
elseif mode == "basal"
     [m_curve,c_ratio,m0] = moment_curve(model,1/resolution/4,mode,-0.5,0.5);
elseif mode == "hydrostatic"
    [m_curve,c_ratio,m0] = moment_curve(model,1/resolution/4,mode,-0.5,h_w);
end
% set up animated plots of stress as crevasse height changes
figure(2)
t = tiledlayout(1,1);
ax2 = axes(t);
ax2_top = axes(t);

for h_c = 1:length(h_c_vect)
    
    if mode == "surface"
        plot(ax2,sigma_net(:,length(h_c_vect)-h_c+1),y_vect)
        yline(ax2,h_c_vect(length(h_c_vect)-h_c+1),'--','h_{sc}','Color','r')
        ax2.Box = 'off';
        xlabel(ax2,"$\sigma_{\Sigma}^\prime$","Interpreter","latex","fontsize",20)
        xlim(ax2,[0,1])
        ylabel(ax2,"Position along ice column")
        yticks(ax2,[-0.5,0,h_w,0.5])
        yticklabels(ax2,["-h_i/2","0","h_w","h_i/2"])
        ylim(ax2,[-0.5,0.5])
        m_curve_flipped = fliplr(m_curve);
        plot(ax2_top,m_curve_flipped(length(h_c_vect)-h_c+1:end)/m0,y_vect(length(h_c_vect)-h_c+1:end),'Color',[0.85,0.325,0.098])
        yticks(ax2_top,[])
        yticklabels(ax2_top,[])
        
    elseif mode == "basal"
        plot(ax2,sigma_net(:,h_c),y_vect)
        yline(ax2,h_c_vect(h_c),'--','h_{bc}','Color','r')
        ax2.Box = 'off';
        xlabel(ax2,"$\sigma_{\Sigma}^\prime$","Interpreter","latex","fontsize",20)
        xlim(ax2,[0,1])
        ylabel(ax2,"Position along ice column")
        yticks(ax2,[-0.5,0,h_w,0.5])
        yticklabels(ax2,["-h_i/2","0","h_w","h_i/2"])
        ylim(ax2,[-0.5,0.5])
        plot(ax2_top,m_curve(1:h_c)/m0,y_vect(1:h_c),'Color',[0.85,0.325,0.098])
        yticks(ax2_top,[])
        yticklabels(ax2_top,[])
        
    elseif mode == "hydrostatic"
        plot(ax2,sigma_net(:,h_c),y_vect)
        yline(ax2,h_c_vect(h_c),'--','h_{bc}','Color','r','LabelVerticalAlignment','bottom')
        yline(ax2,0.5-(0.5-h_w)/length(h_c_vect)*h_c,'--','h_{sc}','Color','r')
        ax2.Box = 'off';
        xlabel(ax2,"$\sigma_{\Sigma}^\prime$","Interpreter","latex","fontsize",20)
        xlim(ax2,[0,1])
        ylabel(ax2,"Position along ice column")
        yticks(ax2,[-0.5,0,h_w,0.5])
        yticklabels(ax2,["-h_i/2","0","h_w","h_i/2"])
        ylim(ax2,[-0.5,0.5])
        plot(ax2_top,m_curve(1:h_c/(0.5+h_w))/m0,y_vect(1:h_c/(0.5+h_w)),'Color',[0.85,0.325,0.098])
        yticks(ax2_top,[-0.5,-0.25,0,0.25,0.5])
        yticklabels(ax2_top,["0%","25%","50%","75%","100%"])
        ylabel(ax2_top,"Percentage of ice thickness fractured")
    end
    
    ax2_top.XAxisLocation = 'top';
    ax2_top.YAxisLocation = 'right';
    ax2_top.Color = 'none';
    ax2_top.Box = 'off';
    title(ax2_top,"Net stress and bending moment as " + mode + " crevasse grows")
    xlabel(ax2_top,"$m/m_0$","Interpreter","latex","fontsize",20)
    xlim(ax2_top,[-1,2])
    xticks(ax2_top,[-0.5,0,0.5,1,1.5])
    ylim(ax2_top,[-0.5,0.5])
    drawnow    
end
