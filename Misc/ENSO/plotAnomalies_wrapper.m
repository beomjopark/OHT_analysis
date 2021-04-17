%% Wrapper function for localMLE and Anomaly Prediction
%% Run after selection / Merge and MeanField
addpath(genpath('../gsw_matlab'));
addpath(genpath('./cbrewer'));

%% Parameters
Params_LatFlux_Step2

typeTag = 'int'
targetVar = 'lon'
isAbsolute = true
isAdjusted = false
iterEM = 3

%% Run Selection
if is2step && strcmp(typeTag, 'int') % This case is for intlatflux/intlonflux
    intStartList = [10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900] % Full upper ocean - 5 % 50:50:500;

    fprintf('Target pressure: %d to %d\n', min(intStartList), max(intStartList));
    typeTag = strcat(typeTag, targetVar); %'intlatlon'
    verticalSelection = strcat(num2str(min(intStartList)),'_',num2str(max(intStartList)));

%    plotAnomaliesMovie(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute);
%   plotAnomaliesHov(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, ...
%    	[190, 240]);

    
    plotIndices;
    close all;
    
    noMaskIdx = 76:(180-76+1);
    [latGrid,longGrid] = meshgrid(linspace(-89.5,89.5,180),linspace(20.5,379.5,360));

    intStartCell = {[10, 900], [10, 100], [100, 300], [300, 900]};
    intStartCell = {[10, 100], [10, 100], [100, 300], [300, 900]};
    pmCell = {};
    nCell = numel(intStartCell);

    figure;
    for ii = 1:nCell
        verticalSelection = strcat(num2str(min(intStartCell{ii})),'_',num2str(max(intStartCell{ii})));
        [dateRange, pmCell{ii}] = plotAnomaliesHov_core(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, ...
            [190, 240]);
        pmCell{ii} = pmCell{ii} / 10^4

        switch targetVar
            case 'lat'
                ylimit = [-155, 155];
                cbLabel = 'W m^{-2}';
            otherwise
                pmCell{ii} = pmCell{ii} ./ 1e+06;
                cbLabel = 'MW m^{-1}';
                ylimit = quantile(pmCell{ii}(:), [0.01 0.99]);        
                ylimit = [-150, 150];
        end        
        
        fig1 = subplot(nCell+1,1,ii);
        imagesc(dateRange,(latGrid(1,noMaskIdx)), pmCell{ii},...
            'AlphaData', double(~isnan(pmCell{ii})));

        caxis(ylimit);
        cLims = caxis;
        colormap(darkb2r(cLims(1), cLims(2)));
        cb = colorbar;
        cb.Label.String = cbLabel;
%        cb.Title.String = 'W m^{-2}';
    
        ax1 = gca;
        ax1.XTick = [dateRange(ONI.MON ==1); dateRange(end) + 20];
        ax1.XGrid = 'on';
        datetick('x','yyyy','keepticks');

        ax1.XAxis.MinorTick = 'on';
        ax1.XAxis.MinorTickValues = dateRange;
    
        ax1.YDir = 'normal';
        ylabel('Latitude');
    
        set(ax1,'fontsize', 14);
        title(['Pressure: ', num2str(min(intStartCell{ii})), ' - ', num2str(max(intStartCell{ii})), ' dbar']);
    end
    
    fig1 = subplot(nCell+1,1,nCell+1);
    plot(ONI.DATE, ONI.ANOM, 'color', 'black', 'LineWidth',2);
    %bar(ONI.DATE, ONI.ANOM, 'color', 'black', 'LineWidth',2);
    hold on;
    plot(ONI.DATE, ONI.ANOM_EL, 'color', 'red', 'LineWidth',2);
    plot(ONI.DATE, ONI.ANOM_LA, 'color', 'blue', 'LineWidth',2);
    hold off;
    cb = colorbar;

    ax1 = gca;
    ax1.XTick = [ONI.DATE(ONI.MON ==1); ONI.DATE(end)+20];
    datetick('x','yyyy','keepticks');
    ax1.XGrid = 'on';
    
    ax1.XAxis.MinorTick = 'on';
    ax1.XAxis.MinorTickValues = ONI.DATE;

    yline(0, 'k', 'LineStyle', '-');
    yline(ulev, 'r', 'LineStyle',':', 'LineWidth',1.5);
    yline(llev, 'b', 'LineStyle',':', 'LineWidth',1.5);

    set(ax1,'fontsize', 14);
    title('Oceanic Nino Index');
else
    for intStartIdx = 1:numel(intStartList)
        intStart = intStartList(intStartIdx); % Parfor requires increasing 1 index

        fprintf('Target pressure: %d\n', intStart);
        verticalSelection = strcat('Relative', num2str(intStart)); %'MidMeso';%'MidMeso';%'UpperOcean';%'Mikael';

        if is2step  % This case is for each latflux / lonflux
%            plotAnomaliesMovie(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute);
            plotAnomaliesDivMovie(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute);
            %            plotAnomaliesHov(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowSize, minNumberOfObs, is2step, isDeriv, targetVar,...
%                159.5);
        else
            plotAnomaliesMovie(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, [], [], isAdjusted, isAbsolute);            
%            computeMeanAnomalies(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, [], []);            
        end
        close all
    end
end