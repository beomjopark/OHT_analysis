function computeMeanAnomalies(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM)

    isPlot = false

    if isempty(isAdjusted)
      isAdjusted = false;
    end
    if isempty(isAbsolute)
      isAbsolute = false;
    end

    switch numel(windowSize)
      case 3
        windowSizeMean = windowSize(1);
        windowSizeCov = windowSize(2);
        windowSizeKrig = windowSize(3);
      case 2
        windowSizeMean = windowSize(1);
        windowSizeCov = windowSize(2);
        windowSizeKrig = windowSize(1);
       otherwise
        windowSizeMean = windowSize;
        windowSizeCov = windowSize;
        windowSizeKrig = windowSize;
    end

    if windowSizeMean == windowSizeCov
        windowSizeTag = num2str(windowSizeMean)
    else
        windowSizeTag = [num2str(windowSizeMean),'_',num2str(windowSizeCov)]
    end
    windowSizeFullTag = [windowSizeTag, '_', num2str(windowSizeKrig)]

    if isempty(windowType)
        windowType = 'box'; % expected 'spherical'
    end

    switch windowType
        case 'spherical'
          windowTypeTag = 'spherical'
          windowSizeMargined = windowSize * 2;
        otherwise
          windowTypeTag = []
          windowSizeMargined = windowSize;
    end

    if isAdjusted
        adjustTag = 'Adjusted';    
        adjustNumTag = ['Adjusted', num2str(nAdjust)];
        windowSizeTag = windowSizeFullTag;
    else
      adjustTag = [];
      adjustNumTag = [];
    end

    if isAbsolute
      absoluteTag = 'Absolute';
    else
      absoluteTag = [];
    end   

    if isempty(iterEM) || (iterEM == 0)
        EMTag = []
    else
        EMTag = ['EM',num2str(iterEM)]
    end

    isFullMonth = (numel(month) == 12)
    if isFullMonth
      EMOutTag = ['Full', EMTag] %% For anomaly, add Full even for 0 iterEM
    else
      EMOutTag = EMTag;
    end
    if isFullMonth && ~isempty(EMTag)
        EMTag = ['Full', EMTag]
    end

    yearRange = strsplit(dataYear, '_');
    startYear = str2num(yearRange{2});
    endYear = str2num(yearRange{3});

    if isDeriv
        switch targetVar
            case 'lon'
                targetLabel = 'Meridional';
            case 'lat'
                targetLabel = 'Zonal';
        end
    else
        targetLabel = [];
    end

    vs = strsplit(verticalSelection, '_');
    isIntFlux = (numel(vs) == 2); %intlat/intlon has intStart_intEnd
    if isIntFlux
      intStart = str2double(vs{1});
      intEnd = str2double(vs{2});
      intMid = mean([intStart, intEnd]);
    else % Relative'intStart'
      intStart = str2double(verticalSelection(9:end));
    end
    % Legacy support
    switch verticalSelection
        case 'MidMeso'
            intStart = 600;
        case 'FullMeso'
            intStart = 200;
    end

    %tag = 'TrendSeason_t2_SpaceExp';
    %tag = 'TrendSeasonSpaceTimeExp';
    %tag = 'TrendSeasonSpaceExp';

    % Load Data mask
    if is2step
        if ~isIntFlux
            data = load(['./Data/','int','TempDens','Prof','PchipPotTemp',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);
            intEnd = data.intEnd;
            clear data;
        end
        load(['./Data/dataMask',typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);
%        load(['./Data/dataMask','target','Temp',verticalSelection,dataYear,'_',num2str(minNumberOfObs),'.mat']);
    else
        % Data here is used to grab intEnd
        data = load(['./Data/',typeTag,'TempDens','Prof','PchipPotTemp',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);
%{
        if numel(data.intStart) > 1
            load(['./Data/dataMask','Relative',num2str(min(data.intStart)),'_',num2str(max(data.intStart)),dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),'_w',num2str(windowSize),'.mat']);    
        else
%}
            load(['./Data/dataMask',verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);    
%        end
        intEnd = data.intEnd;
        clear data;
    end
    maskJohn = ncread('./RG_climatology/RG_ArgoClim_Temperature_2016.nc','BATHYMETRY_MASK',[1 1 25],[Inf Inf 1]);
    maskJohn(maskJohn == 0) = 1;
    maskJohn = [NaN*ones(360,25) maskJohn NaN*ones(360,25)];
    mask = maskJohn .* dataMask;

    cp0 = gsw_cp0; % as in McDougall 2003
    rho0 = 1030;


    % Conversion: Deg / m
    [latGrid,longGrid] = meshgrid(linspace(-89.5,89.5,180),linspace(20.5,379.5,360));
    latDistGrid = 1 ./ distance(latGrid - 0.5, longGrid, latGrid + 0.5, longGrid,...
                            referenceEllipsoid('GRS80', 'm'));
    longDistGrid = 1 ./ distance(latGrid, longGrid - 0.5, latGrid, longGrid + 0.5,...
                            referenceEllipsoid('GRS80', 'm'));

    %% Movie goes here
    % Save directory
    if is2step
        destFolder = [typeTag,fluxType,responseTag,verticalSelection,dataYear,adjustNumTag,absoluteTag,windowTypeTag,'_Eq',num2str(eqBorder),...
                 '/','Pre_','Anomaly_',kernelType,'_',windowSizeTag]
        srcFolder = ['anomaly_',typeTag,fluxType,responseTag,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeFullTag,'_Eq',num2str(eqBorder),'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag];
    else
        profileTag = []; %''_Profile_YF3''
        destFolder = [typeTag,responseTag,verticalSelection,dataYear,adjustNumTag,absoluteTag,windowTypeTag, meanTag,profileTag,...
                 '/','Pre_','Anomaly_',kernelType,'_',windowSizeTag,'_month',num2str(month)]
        srcFolder = ['anomaly_',typeTag,responseTag,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeFullTag,'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag];
%        srcFolder = ['anomaly_',typeTag,responseTag,'_w',num2str(windowSize),'_Eq',num2str(eqBorder),'_',kernelType,'_',verticalSelection];
    end
    mkdir(['./Figures/',destFolder])

    nNaN = zeros(size(latGrid));
    meanPredGrid = zeros(size(latGrid));
    for iYear = startYear:endYear
        for iMonth = 1:12

            if isDeriv
                % load(['./Results/',srcFolder,'/anomaly_',responseTag,'_w',num2str(windowSize),'_',kernelType,targetVar,'Deriv_',verticalSelection,...
                % '18/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,targetVar,'Deriv_',num2str(iMonth,'%02d'),'_',num2str(iYear),'.mat']);
                load(['./Results/',srcFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,targetVar,'Deriv_',num2str(iMonth,'%02d'),'_',num2str(iYear),'.mat']);
               % switch targetVar
               %     case 'lat'
               %         predGrid = predGrid .* latDistGrid;
               %     case 'lon'
               %         predGrid = predGrid .* longDistGrid;
               % end
            else
                if is2step
                    % Equitorial Mask
                    eqmask = ~(latGrid < eqBorder & latGrid > - eqBorder) .* 1;
                    eqmask(eqmask == 0) = NaN;
                    mask = mask .* eqmask;
                    
                    load(['./Results/',srcFolder,'/anomaly',typeTag,fluxType,responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'_',...
                            num2str(iMonth,'%02d'),'_',num2str(iYear),'.mat']);
%{
                    switch fluxType
                        case 'heat'
                            predGrid = gsw_cp0 .* mask .* predGrid.* 1e-15;
                        case 'mass'
                            predGrid = mask .* predGrid;
                        case 'vol'
                            predGrid = mask .* predGrid.* 1e-6;
                    end
%}
                else
                    load(['./Results/',srcFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'_',...
                            num2str(iMonth,'%02d'),'_',num2str(iYear),'.mat']);    
                end
            end

            % Suppress NaN not in datamask
            isNaNpred = isnan(predGrid) & ~isnan(mask);
            predGrid(isNaNpred) = 0;
            
            nNaN = nNaN + isNaNpred;
            meanPredGrid = meanPredGrid + predGrid;
        end
    end
  
    meanPredGrid = meanPredGrid ./ ((endYear - startYear +1).*12 - nNaN);

    if isDeriv
        save(['./Results/',srcFolder,'/MeanAnomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,targetVar,'Deriv','.mat'], 'meanPredGrid');
    else
        if is2step
            save(['./Results/',srcFolder,'/MeanAnomaly',typeTag,fluxType,responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'.mat'], 'meanPredGrid');
        else
            save(['./Results/',srcFolder,'/MeanAnomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'.mat'], 'meanPredGrid');    
        end
    end

    %% PLOTTING
    if isPlot
        if isAbsolute && ~isIntFlux
            presLabel = num2str(intStart);
        else
            presLabel = [num2str(intStart),'-',num2str(intEnd)];
        end

        if isAbsolute
            presLabel = ['Absolute ', presLabel];
        end

        % Check NaNs
        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');
        
        pos = get(gca,'position');
        set(gca,'position',pos + [-0.13 -0.025 0.2 0])

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 22.5 15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))

        pos = get(gcf,'position');
        set(gcf,'position',pos + [0 0 0 -4])

        clf;

        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');
                 
        pos = get(gca,'position');
        set(gca,'position',pos + [-0.13 -0.025 0.2 0])
        
        surfm(latGrid,longGrid, nNaN./ ((endYear - startYear +1).*12));
        colorbar()
        load coast;
        plotm(lat,long,'k');
        
        title('NaN Ratio');    
        print('-depsc2',['./Figures/',destFolder,'/',targetLabel,'MeanMonthlyAnomalyNaNRatio','.eps']);

        % Plot meanPredGrid
        figure;
        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 22.5 15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))

        pos = get(gcf,'position');
        set(gcf,'position',pos + [0 0 0 -4])

        clf;

        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');
                 
        pos = get(gca,'position');
        set(gca,'position',pos + [-0.13 -0.025 0.2 0])
        
        ylimit = quantile(mask(:).*meanPredGrid(:), [0.005 0.995]);
        switch responseTag
            case 'Flux'
                surfm(latGrid,longGrid,mask.*meanPredGrid); 
                title(['Monthly Averaged ',targetLabel,' ', fluxType,' Transport anomaly, ',presLabel,'dBar, ']);    
                caxis(ylimit);
            case 'Dens'
                surfm(latGrid,longGrid,mask.*meanPredGrid); 
                
                if isDeriv
                    switch verticalSelection
                        case 'MidMeso'
                            title(['Monthly Averaged ',targetLabel, ' derivative Density anomaly, 600-1000 m']);
                            switch targetVar
                                case 'lat'
                                    caxis([-15, 15]);
                                case 'lon'
                                    caxis([-20, 20]);
                            end                        
                        case 'FullMeso'
                            title(['Monthly Averaged ',targetLabel, ' derivative Density anomaly, 200-1000 m']);
                            switch targetVar
                                case 'lat'
                                    caxis([-45, 45]);
                                case 'lon'
                                    caxis([-45, 45]);
                            end 
                        otherwise % 'Relative'
                            title(['Monthly Averaged ',targetLabel, ' derivative Dynamic Height anomaly, ',presLabel,'dBar']);
                            caxis(ylimit);
                    end
                else
                    switch verticalSelection
                        case 'MidMeso'
                            title(['Monthly Averaged ','Density anomaly, 600-1000 m']);
                            caxis([-30, 30]);
                        case 'FullMeso'
                            title(['Monthly Averaged ','Density anomaly, 200-1000 m']);
                            caxis([-60, 60]);
                        otherwise % 'Relative'
                            title(['Monthly Averaged ','Dynamic Height anomaly, ',presLabel,'dBar']);
                            caxis(ylimit);
                    end
                end
    %            caxis([quantile(predGrid(:), 0.01), quantile(predGrid(:), 0.99)])
            case 'Temp'
                surfm(latGrid,longGrid,cp0.*rho0.*mask.*meanPredGrid);
                title(['Ocean heat content anomaly, 0-2000 m, ']);
                %caxis([-400,400]);
                caxis([-1.5e9,1.5e9]);
        end
           
        load coast;
        plotm(lat,long,'k');
        
        cLims = caxis;  
        cb = colorbar;
        cb.FontSize = 14;
        if is2step
            switch fluxType
                case 'heat'
                    cb.Label.String = 'Heat Transport [PW]';
                case 'mass'
                    cb.Label.String = 'Mass Transport [kg / s]';
                case 'vol'
                    cb.Label.String = 'Volume Transport [Sv]';
            end
        else
            if isDeriv
    %                cb.Label.String = 'Density/deg [kg / m^3 deg]';        
                cb.Label.String = 'm/s^2';        
            else
    %                cb.Label.String = 'Density, kg/m^3';        
                cb.Label.String = 'm^2/s^2';        
            end
        end
        colormap(darkb2r(cLims(1),cLims(2)));

        drawnow;


        if isDeriv
            print('-depsc2',['./Figures/',destFolder,'/',targetLabel,'MeanMonthlyAnomaly','.eps']);
        else
            print('-depsc2',['./Figures/',destFolder,'/','MeanMonthlyAnomaly','.eps']);
        end
    end
end
