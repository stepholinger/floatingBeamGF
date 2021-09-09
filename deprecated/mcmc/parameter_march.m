% load mcmc results
load('/home/setholinger/Documents/Projects/PIG/modeling/mcmc/RC/run4_results.mat')

% construct some variables
statDist = 10000;
t0 = 5.5;
f_max = 0.5;
t_max = 250;
t = 1/(2*f_max):1/(2*f_max):t_max;
nt = t_max*(2*f_max);

% set flag for toy problem
test = 0;

% get real data
fname = "/media/Data/Data/PIG/MSEED/noIR/PIG2/HHZ/2012-04-02.PIG2.HHZ.noIR.MSEED";
dataStruct = rdmseed(fname);

% extract trace
trace = extractfield(dataStruct,'d');
fs = 100;

% resample data
fsNew = f_max*2;
trace = resample(trace,fsNew,100);

% set event bounds
startTime = ((15*60+18)*60+50)*fsNew;
endTime = startTime + nt;

% trim data to event bounds
eventTrace = trace(startTime:endTime-1);

% remove scalar offset using first value
eventTrace = eventTrace - eventTrace(1);

%x0_vect = [80,600,6500,6,0.5,250;...
%         85,450,6000,6,0.5,250;...
%         100,450,6500,7.5,0.5,250;...
%         80,600,5750,6,0.5,250;...
%         80,600,6500,6.75,0.5,250;...
%         80,450,7000,6,0.5,250];

x0_vect = [250,600,6500,6,0.5,250];%...
         %85,450,6000,6,0.5,250;...
         %100,450,6500,7.5,25,250;...
         %80,600,5750,6,0.5,250;...
         %80,600,6500,6.75,0.5,250;...
         %80,450,7000,6,0.5,250];
     
axisLabels = ["Ice thickness (m)", "Water depth (m)", "X_{stat} (km)","t_0 (s)"];
paramLabels = ["h_i","h_w","Xstat","t0"];
%paramsVariedVect = [1,2;1,3;1,4;2,3;2,4;3,4];
paramsVariedVect = [1,2];
%stepVect = [10,100;10,400;10,0.5;100,400;-100,0.075;-700,0.5];
stepVect = [25,0;0,400;0,10;100,0;0,0;0,0];
numSteps = 5;
maxNumBins = 100;

% set up colored lines
c = hot(floor(numSteps*1.75));

for n = 1:size(paramsVariedVect,1)

    paramsVaried = paramsVariedVect(n,:);
    step = [0,0,0,0,0,0];
    step(paramsVaried(1)) = stepVect(n,1);
    step(paramsVaried(2)) = stepVect(n,2);

    % get starting point from x0
    x0 = x0_vect(n,:);
    x_start = x0 - numSteps/2*step;
    
    t = 1/(2*f_max):1/(2*f_max):t_max;    
    nt = t_max*(2*f_max);
    eventTraceTrim = eventTrace(1:nt);
    
    % get number of bins
    numBins = length(unique(x_keep(paramsVaried(1),:)));
    if numBins > maxNumBins
        numBins = maxNumBins;
    end
    
    % make figure
    figure(n)
    subplot(1,4,1:2)
    hold on; 
    
    % make dscatter
    dscatter(x_keep(paramsVaried(1),:)',x_keep(paramsVaried(2),:)','BINS',[numBins,numBins]);

    % force full outlines
    box on;

    % remove x-axis labels if current subplot does not fall along
    % edge of plot grid
    xlabel(axisLabels(paramsVaried(1)));
    ylabel(axisLabels(paramsVaried(2)));
    
    % make array for storing model output
    G_mat = zeros(numSteps,nt);
    
    for i = 1:numSteps

        % get point
        x = x_start + (i-1) * step;
        
        % generate Green's function
        [G,eventAlign,M_fit] = GF_func_mcmc(x,eventTraceTrim);
        
        % fill matrix of model output
        G_mat(i,:) = G;
        
        % convert x_stat to km
        x(3) = x(3)/1000;
        
        % plot points on dscatter plot
        subplot(1,4,1:2)
        scatter(x(paramsVaried(1)),x(paramsVaried(2)),50,c(i,:),'filled')
        hold on;
        
        %subplot(1,4,3:4)
        %plot(G_mat(i,:),'Color',c(i,:))
        %hold on;
    end
    
    G_align = zeros(numSteps,nt);

    % make plot of 1st trace
    subplot(1,4,3:4)
    plot(G_mat(1,:),'Color',c(1,:))
    legend;
    hold on;
    % align traces with cross correlation
    for g = 2:numSteps

        % cross correlate
        [G_xcorr,lags] = xcorr(G_mat(1,:),G_mat(g,:));
        [coef,lagIdx] = max(abs(G_xcorr));
        lag = lags(lagIdx);
        if lag <= 0
            G_align = [zeros(1,abs(lag)),G_mat(g,:)];
        end

        if lag > 0
            G_align = G_mat(g,abs(lag):end);
        end
            
         % make plot of all traces
         subplot(1,4,3:4)
         plot(G_align,'Color',c(g,:))
         legend;
         hold on;
        
    end
    
end