function [dateRange, plotMatrix] = plotAnomaliesHov_core(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, lonDeg)

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
          isSpherical = true;
          windowTypeTag = 'spherical'
          windowSizeMargined = windowSizeKrig * 2;
        otherwise
          isSpherical = false;
          windowTypeTag = []
          windowSizeMargined = windowSizeKrig;
    end

    if isAdjusted
      adjustTag = 'Adjusted';
      adjustNumTag = ['Adjusted', num2str(nAdjust)];
      windowSizeTag = windowSizeFullTag;
    else
      adjustTag = [];
      adjustNumTag = [];
      isManualAdjust = true;
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
    nYear = endYear - startYear + 1;

    switch targetVar
        case 'lon'
            targetLabel = 'Meridional';
        case 'lat'
            targetLabel = 'Zonal';
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
            data = load(['./Data/','int','TempDens','Prof','PchipPotTemp',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'_w',num2str(windowSize),'.mat']);
            intEnd = data.intEnd;
            clear data;
        end
        load(['./Data/dataMask',typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowType,'_w',windowSizeTag,'.mat']);
%        load(['./Data/dataMask','target','Temp',verticalSelection,dataYear,'_',num2str(minNumberOfObs),'.mat']);
    else
        data = load(['./Data/',typeTag,'TempDens','Prof','PchipPotTemp',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'_w',num2str(windowSize),'.mat']);
        load(['./Data/dataMask',verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),'_w',num2str(windowSize),'.mat']);    
        intEnd = data.intEnd;
        clear data;
    end
    maskJohn = ncread('./RG_climatology/RG_ArgoClim_Temperature_2016.nc','BATHYMETRY_MASK',[1 1 25],[Inf Inf 1]);
    maskJohn(maskJohn == 0) = 1;
    maskJohn = [NaN*ones(360,25) maskJohn NaN*ones(360,25)];
    mask = maskJohn .* dataMask;

    cp0 = gsw_cp0; % as in McDougall 2003
    rho0 = 1030;

    %% Movie goes here
    % Save directory
    if is2step
        srcFolder = ['anomaly_',typeTag,fluxType,responseTag,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeFullTag,'_Eq',num2str(eqBorder),'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag]
        if isAdjusted || isManualAdjust
          meanPred = load(['./Results/',srcFolder,'/MeanAnomaly',typeTag,fluxType,responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'.mat']);
        end
    else
        profileTag = []; %''_Profile_YF3''
        srcFolder = ['anomaly_',typeTag,responseTag,adjustTag,absoluteTag,'_w',num2str(windowSize),'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d')];
%        srcFolder = ['anomaly_',typeTag,responseTag,'_w',num2str(windowSize),'_Eq',num2str(eqBorder),'_',kernelType,'_',verticalSelection];
    end
    
    % Area computation
    [latGrid,longGrid] = meshgrid(linspace(-89.5,89.5,180),linspace(20.5,379.5,360));
    latDistGrid = distance(latGrid - 0.5, longGrid, latGrid + 0.5, longGrid,...
                            referenceEllipsoid('WGS84', 'm'));
    longDistGrid = distance(latGrid, longGrid - 0.5, latGrid, longGrid + 0.5,...
                            referenceEllipsoid('WGS84', 'm'));
    areaGrid = areaquad(latGrid - 0.5, longGrid -0.5, latGrid + 0.5, longGrid+0.5,...
                            referenceEllipsoid('WGS84', 'm'));
    areaGrid = areaGrid ./ latDistGrid; % Per Latitude!
    
    % Convert total area
    meanPred.meanPredGrid = meanPred.meanPredGrid .* areaGrid;

    predGridCube = zeros([size(latGrid), 12*nYear]);
    dateRange = zeros(12*nYear, 1);

    cnt = 1;
    for iYear = startYear:endYear
        for iMonth = 1:12

            if isDeriv
                S = load(['./Results/',srcFolder,'/anomaly_',responseTag,'_w',num2str(windowSize),'_',kernelType,targetVar,'Deriv_',verticalSelection,...
                '18/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,targetVar,'Deriv_',num2str(iMonth,'%02d'),'_',num2str(iYear),'.mat']);
                S.predGrid = S.predGrid .* areaGrid;
                predGrid = S.predGrid;

               % switch targetVar
               %     case 'lat'
               %         predGrid = predGrid .* latDistGrid;
               %     case 'lon'
               %         predGrid = predGrid .* longDistGrid;
               % end
            else
                if is2step
                    S = load(['./Results/',srcFolder,'/anomaly',typeTag,fluxType,responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'_',...
                            num2str(iMonth,'%02d'),'_',num2str(iYear),'.mat']);
                    S.predGrid = S.predGrid .* areaGrid;
                    
                    % Compute OHT
                    if isAdjusted || isManualAdjust
                        S.predGrid = S.predGrid - meanPred.meanPredGrid;
                    end
                    predGrid = gsw_cp0 .* S.predGrid;%.* 1e-15;
                    
                    % Equitorial Mask
                    eqmask = ~(latGrid < eqBorder & latGrid > - eqBorder) .* 1;
                    eqmask(eqmask == 0) = NaN;
                    mask = mask .* eqmask;
                else
                    S = load(['./Results/',srcFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'_',...
                            num2str(iMonth,'%02d'),'_',num2str(iYear),'.mat']);    
                    predGrid = S.predGrid;
                    if isAdjusted
                        S.predGrid = S.predGrid - meanPred.meanPredGrid;
                    end
                end
            end

            predGridCube(:,:,cnt) = mask.*predGrid;
            dateRange(cnt) = S.midJulDay;
            cnt = cnt + 1;
        end
    end
%    dateRangeTime = datetime(dateRange,'ConvertFrom','datenum');    
    if(numel(lonDeg) > 1)
        
        lonIdx = find(longGrid(:,1) > lonDeg(1) & longGrid(:,1) < lonDeg(2));
    
        % Check mask range
        chkMask = longDistGrid(lonIdx, :) .* mask(lonIdx,:);
        chkMask(isnan(chkMask)) = false;
        
        [~, fc] = find(chkMask, 1, 'first');
        [~, lc] = find(chkMask, 1, 'last');
        
        fc = max(fc, 75);
        lc = min(lc, 180-75+1);
        noMaskIdx = fc:lc;
        predGridCube(isnan(predGridCube)) = 0;
    
        plotMatrix = squeeze(predGridCube(lonIdx,noMaskIdx,:));

        if strcmp(targetVar, 'lat')
            plotMatrix = squeeze(sum(plotMatrix, 1)) ./ sum(chkMask(:, noMaskIdx), 1)'; % Averaged
        else
            plotMatrix = squeeze(sum(plotMatrix, 1));        
        end
            
    else
        lonIdx = find(longGrid(:,1) == lonDeg);
        
        % Check mask range
        chkMask = mask(lonIdx,:);
        chkMask(isnan(chkMask)) = false;
        noMaskIdx = find(chkMask, 1, 'first'):find(chkMask, 1, 'last');
        
        plotMatrix = squeeze(predGridCube(lonIdx,noMaskIdx,:));    
    end
        
end