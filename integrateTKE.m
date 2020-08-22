%% time series of the Total Kinetic Energy
addpath(genpath('../gsw_matlab'));

Params_LatFlux_Step1
iterEM = 3
isAbsolute = true
eqBorder = 3

intStartList = [10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800]

%function integrateTKE(intStartList, kernelType, month, typeTag, responseTag, dataYear, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute)
 
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
          windowSizeMargined = windowSizeCov * 2;
        otherwise
          windowTypeTag = []
          windowSizeMargined = windowSizeCov;
    end

    if isempty(isAdjusted)
      isAdjusted = false;
    end
    if isempty(isAbsolute)
      isAbsolute = false;
    end

    if isAdjusted
      adjustTag = 'Adjusted';
      adjustNumTag = ['Adjusted', num2str(nAdjust)];
      if nAdjust == 1
          adjustPrevNumTag = [];
      else
          adjustPrevNumTag = ['Adjusted', num2str(nAdjust-1)];
      end
      windowSizeTag = windowSizeFullTag;
    else
      adjustTag = [];
      adjustNumTag = [];
      adjustPrevNumTag = [];
    end

    if isAbsolute
      absoluteTag = 'Absolute';
    else
      absoluteTag = [];
    end   

    if isempty(iterEM) || (iterEM == 0)
        EMTag = []
        prevEMTag = []
    elseif iterEM == 1
        EMTag = ['EM',num2str(iterEM)]
        prevEMTag = []
    else
        EMTag = ['EM',num2str(iterEM)]
        prevEMTag = ['EM',num2str(iterEM-1)]
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

    nTargetPres = numel(intStartList);
    intStart = min(intStartList);
    intEnd = max(intStartList);
    verticalSelection = strcat('Relative', num2str(intStart));
    intTag = [num2str(intStart),'_', num2str(intEnd)];
    gridz = linspace(intStart, intEnd, 5000); % Fine grid, takes about 13 mins to run

    yearRange = strsplit(dataYear, '_');
    startYear = str2num(yearRange{2});
    endYear = str2num(yearRange{3});
    nYear = endYear - startYear + 1;

    % Load Data mask
    if is2step
        if ~isIntFlux
            data = load(['./Data/','int','TempDens','Prof','PchipPotTemp',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'_w',num2str(windowSizeMean),'.mat']);
            intEnd = data.intEnd;
            clear data;
        end
        load(['./Data/dataMask',typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),'_w',num2str(windowSizeMean),'.mat']);
    else
        % Most conservative mask
        verticalSelection = strcat('Relative', num2str(intEnd));
        load(['./Data/dataMask',verticalSelection,dataYear,adjustTag,[],'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);    
    end
    maskJohn = ncread('./RG_climatology/RG_ArgoClim_Temperature_2016.nc','BATHYMETRY_MASK',[1 1 25],[Inf Inf 1]);
    maskJohn(maskJohn == 0) = 1;
    maskJohn = [NaN*ones(360,25) maskJohn NaN*ones(360,25)];

    % Equitorial Mask
    [latGrid,longGrid] = meshgrid(linspace(-89.5,89.5,180),linspace(20.5,379.5,360));
    nGrid = numel(latGrid);
    eqmask = ~(latGrid < eqBorder & latGrid > - eqBorder) .* 1;
    eqmask(eqmask == 0) = NaN;

    % Determine Region from the DataMask
    mask = maskJohn .* dataMask;
    mask = mask .* eqmask;

    tag = 'PchipPotTemp';
    YFtag = []
    meanTag = 'NoTrend';
    if is2step
        destFolder = [typeTag,fluxType,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_Eq',num2str(eqBorder),...
                 '/','Pre_','Anomaly_',kernelType,'_',num2str(windowSize)];
        srcFolder = ['anomaly_',typeTag,fluxType,responseTag,adjustTag,absoluteTag,'_w',num2str(windowSize),'_Eq',num2str(eqBorder),'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d')]
%        srcFolder = ['anomaly_',typeTag,fluxType,responseTag,adjustTag,absoluteTag,'_w',num2str(windowSize),'_Eq',num2str(eqBorder),'_',kernelType,'_',verticalSelection];

        if isAdjusted
          meanPred = load(['./Results/',srcFolder,'/MeanAnomaly',typeTag,fluxType,responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'.mat']);
        end

        % MeanField
        load(['./Results/meanField',typeTag, fluxType,responseTag,meanTag,tag,verticalSelection,dataYear,'AdjustedAbsolute','_w',num2str(windowSize),'_',num2str(minNumberOfObs),'_Eq',num2str(eqBorder),'.mat'])
    else
        % Load Mean Field
        nHar = 6;
        iLat = 1 + (2*nHar) + 1;
        iLong = iLat + 1;
        refPres = 900

        if isAbsolute
            refPres = 900
            load(['./Misc/AGVA/','meanZonVelRef',num2str(refPres),'.mat'], 'meanRefZonVel');
            load(['./Misc/AGVA/','meanMerVelRef',num2str(refPres),'.mat'], 'meanRefMerVel');
        else
            meanRefMerVel = zeros(size(latGrid));
            meanRefZonVel = zeros(size(latGrid));
        end

        destFolder = [typeTag,responseTag,intTag,dataYear,adjustTag,absoluteTag, meanTag,...
                 '/','Pre_','Anomaly_',kernelType,'_',windowSizeFullTag,'_month',num2str(month)]
%        srcFolder = ['anomaly_',typeTag,responseTag,'_w',num2str(windowSize),'_Eq',num2str(eqBorder),'_',kernelType,'_',verticalSelection];
    end
    mkdir(['./Figures'])
    mkdir(['./Figures/',typeTag,responseTag,intTag,dataYear,adjustTag,absoluteTag, meanTag])
    mkdir(['./Figures/', destFolder])


    distLatGrid = 1 ./ distance(latGrid - 0.5, longGrid, latGrid + 0.5, longGrid,...
            referenceEllipsoid('GRS80', 'm'));
    distLonGrid = 1 ./ distance(latGrid, longGrid - 0.5, latGrid, longGrid + 0.5,...
            referenceEllipsoid('GRS80', 'm'));

    % Integration scale factor 1/g
    scaleGridInt = zeros([nTargetPres, size(latGrid)]);
    for presIdx = 1:nTargetPres
        scaleGridInt(presIdx, :, :) = 1 ./ gsw_grav(latGrid, intStartList(presIdx));%targetFlux.profScaleInt;
    end

    predGridCube = zeros([size(latGrid), 12*nYear]);
    dateRange = zeros(12*nYear, 1);
    cnt = 1;
    for iYear = startYear:endYear
        for iMonth = 1:12
            disp([num2str(iYear),' / ', num2str(iMonth)])

            % Initialize
            TkeAbsGrid = zeros([nTargetPres, size(latGrid)]);
            for presIdx = 1:nTargetPres
                verticalSelection = strcat('Relative', num2str(intStartList(presIdx)));
                % Load Mean Estimate
                load(['./Results/meanField',responseTag,meanTag,tag,verticalSelection,dataYear,windowTypeTag,'_w',windowSizeTag,'_',num2str(minNumberOfObs),YFtag,EMTag,'.mat']);

                % Load Mean Krig from Prev Adjustment
                if isAdjusted
                    srcFolder = ['anomaly_',typeTag,responseTag,adjustPrevNumTag,windowTypeTag,'_w',windowSizeFullTag,'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag];
                    % Load meanPredGrid
                    meanLatpred = load(['./Results/',srcFolder,'/MeanAnomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'lat','Deriv','.mat']);
                    meanLonpred = load(['./Results/',srcFolder,'/MeanAnomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'lon','Deriv','.mat']);
                else
                    meanLatpred = struct('meanPredGrid', zeros(size(latGrid)));
                    meanLonpred = struct('meanPredGrid', zeros(size(latGrid)));
                end

                % Load Kriging
                srcFolder = ['anomaly_',typeTag,responseTag,adjustNumTag,[],windowTypeTag,'_w',windowSizeFullTag,'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag];
                SLat = load(['./Results/',srcFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'lat','Deriv_',num2str(iMonth,'%02d'),'_',num2str(iYear),'.mat']);
                SLon = load(['./Results/',srcFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'lon','Deriv_',num2str(iMonth,'%02d'),'_',num2str(iYear),'.mat']);

                %% CAUTION: Typically Kriging is unstable than the Mean, either because parameter is unstable or no profile. For now, we replace missing to zero.
                SLat.predGrid(isnan(SLat.predGrid)) = 0;
                SLon.predGrid(isnan(SLat.predGrid)) = 0;
                
                % Compute TKE
                zonRelVel = - (betaGrid(:,:,iLat) + meanLatpred.meanPredGrid + SLat.predGrid) .* distLatGrid ./ gsw_f(latGrid);
                merRelVel = (betaGrid(:,:,iLong) + meanLonpred.meanPredGrid + SLon.predGrid) .* distLonGrid ./ gsw_f(latGrid);
                TkeAbsGrid(presIdx, :, :) = (meanRefZonVel + zonRelVel).^ 2 + (meanRefMerVel + merRelVel) .^ 2;
            end

            % Compute Integration
            TkeAbsIntGrid = zeros(size(latGrid));
            for iGrid = 1:nGrid
                [iGridSub1, iGridSub2] = ind2sub(size(latGrid), iGrid);
                if isnan(mask(iGrid)) || all(isnan(TkeAbsGrid(:, iGridSub1, iGridSub2)))
                    TkeAbsIntGrid(iGrid) = NaN;
                    continue;
                end

                pchip_int_tke = interp1(intStartList, TkeAbsGrid(:, iGridSub1, iGridSub2) .* scaleGridInt(:, iGridSub1, iGridSub2), gridz, 'pchip');
                TkeAbsIntGrid(iGrid)  = trapz(gridz, pchip_int_tke);
            end

            predGridCube(:,:,cnt) = TkeAbsIntGrid;
            dateRange(cnt) = SLat.midJulDay;
            cnt = cnt + 1;
        end
    end
    
    %% Plot
    figure;

    dateRangeTime = datetime(dateRange,'ConvertFrom','datenum');
    dateMonth = datevec(dateRangeTime);
    dateMonth = dateMonth(:,2);
    
    predGridCube(isnan(predGridCube)) = 0;
    predGridAggr = squeeze(sum(predGridCube, [1,2]));
    plot(dateRangeTime(1:(end-1)), predGridAggr(1:(end-1)), 'LineWidth',2)
    
    ax1 = gca;
    ax1.XTick = [dateRangeTime(dateMonth ==1); dateRangeTime(end)+20];
    datetick('x','yyyy','keepticks');
    ax1.XGrid = 'on';
 
    ax1.XAxis.MinorTick = 'on';
    ax1.XAxis.MinorTickValues = dateRangeTime;
    
    title('Absolute TKE');
    print('-depsc2',['./Figures/',destFolder,'/','TKE',EMOutTag,'.eps']);

%end