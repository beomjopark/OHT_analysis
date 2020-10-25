function localMLESpaceTimeSeason(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isProgressbar, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth)
%% Fit Local GP with specified kernelType using MLE
    if nargin < 11 || isempty(isStandardize)
        isStandardize = false;
    end
    if nargin < 13 || isempty(eqBorder)
        isEqBorder = false;
    else
        isEqBorder = true;
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

    if isFullMonth && ~isempty(EMTag)
        EMTag = ['Full', EMTag]
    end


    yearRange = strsplit(dataYear, '_');
    startYear = str2double(yearRange{2});
    endYear = str2double(yearRange{3});
    nYear = endYear - startYear + 1;

    [latGrid,longGrid] = meshgrid(linspace(-89.5,89.5,180),linspace(20.5,379.5,360));
    nGrid = numel(latGrid);

    thetasOpt = zeros(size(latGrid));
    thetaLatOpt = zeros(size(latGrid));
    thetaLongOpt = zeros(size(latGrid));
    thetatOpt = zeros(size(latGrid));
    sigmaOpt = zeros(size(latGrid));
    nll = zeros(size(latGrid));
    nResGrid = zeros(size(latGrid));
    exitFlags = zeros(size(latGrid));
    isFminError = zeros(size(latGrid));

    % Load Data mask

    if is2step
        if isnumeric(verticalSelection) % intlat/intlon
            targetPres = verticalSelection;
            presString = [num2str(min(targetPres)),'_',num2str(max(targetPres))];
            load(['./Data/dataMask',typeTag,responseTag,presString,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);
        else
            % For each target Pressure
            load(['./Data/dataMask',typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);
        end
    else
%{
        if numel(intStart) > 1
            load(['./Data/dataMask','Relative',num2str(min(intStart)),'_',num2str(max(intStart)),dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);    
        else
%}
        load(['./Data/dataMask',verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);    
%        end
%        clear dat;
    end
    maskJohn = ncread('./RG_climatology/RG_ArgoClim_Temperature_2016.nc','BATHYMETRY_MASK',[1 1 25],[Inf Inf 1]);
    maskJohn(maskJohn == 0) = 1;
    maskJohn = [NaN*ones(360,25) maskJohn NaN*ones(360,25)];
    mask = maskJohn .* dataMask;

    if strcmp(windowType, 'spherical')
        % Determine reference distance
        refDist = distance(0, 180, 0+windowSizeCov, 180,...
                                referenceEllipsoid('GRS80', 'm'))
    end
    
    % Discard previous itearIdx, if it exists
    %{
    fileName = ['iterIdxLocalMLESpaceTime',kernelType,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(windowSize),'.txt'];
    fileID = fopen(fileName,'w');
    fclose(fileID);
    %}
    if isProgressbar
        parfor_progress(nGrid);
    end
    

    tic;
    parfor iGrid = 1:nGrid
    %for iGrid = 1:nGrid
        fprintf('%d\n', iGrid);
%{
        fileID = fopen(fileName,'a');
        fprintf(fileID,'%d \n',iGrid);
        fclose(fileID);
%}

        predLat = latGrid(iGrid);
        predLong = longGrid(iGrid);

        latMin = predLat - windowSizeMargined;
        latMax = predLat + windowSizeMargined;
        longMin = predLong - windowSizeMargined;
        longMax = predLong + windowSizeMargined;

        if isEqBorder
            if abs(predLat) <= eqBorder
                fprintf('iGrid: %d is in Equator Border\n', iGrid);
                thetasOpt(iGrid) = NaN;
                thetaLatOpt(iGrid) = NaN;
                thetaLongOpt(iGrid) = NaN;
                thetatOpt(iGrid) = NaN;
                sigmaOpt(iGrid) = NaN;
                nll(iGrid) = NaN;
            
                if isProgressbar
                    parfor_progress;
                end
                continue;
            end
            if predLat > 0
                latMin = max([latMin, eqBorder]);
            else
                latMax = min([latMax, - eqBorder]);
            end
        end

        % Do not compute mean if outside land/datamask
    %    [longIdx,latIdx] = ind2sub(sizeGrid, iGrid);
        if isnan(mask(iGrid)) 
            thetasOpt(iGrid) = NaN;
            thetaLatOpt(iGrid) = NaN;
            thetaLongOpt(iGrid) = NaN;
            thetatOpt(iGrid) = NaN;
            sigmaOpt(iGrid) = NaN;
            nll(iGrid) = NaN;
            
            if isProgressbar
                parfor_progress;
            end
            continue;             
        end
        
        profLatAggr = cell(1,nYear);
        profLongAggr = cell(1,nYear);
        profJulDayAggr = cell(1,nYear);
        responseResAggr = cell(1,nYear);        
        for iYear = startYear:endYear
            S = load(['./Data/Extended/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(month,'%02d'),'_',num2str(iYear),'_extended.mat']);
            
            profLat3Months = S.profLatAggr3Months;
            profLong3Months = S.profLongAggr3Months;
            profJulDay3Months = S.profJulDayAggr3Months;
            switch typeTag
                case 'int'
                    responseRes3Months = S.intDensRes3Months;
                case {'intlat', 'intlon'}
                    responseRes3Months = S.intFluxRes3Months;
                otherwise % lat lon
                    responseRes3Months = S.FluxRes3Months;
            end

            idx = find(profLat3Months > latMin & profLat3Months < latMax & profLong3Months > longMin & profLong3Months < longMax);
            switch windowType
                case 'spherical'
                  is_in_circle = distance(predLat, predLong, profLat3Months(idx), profLong3Months(idx),...
                                    referenceEllipsoid('GRS80', 'm')) < refDist;
                  idx = idx(is_in_circle);
                case 'box_var'
                  idx = find(profLatAggrSel > latMin & profLatAggrSel < latMax & profLongAggrSel > longMin & profLongAggrSel < longMax);
            end

            profLatAggr{iYear-startYear+1} = profLat3Months(idx)';
            profLongAggr{iYear-startYear+1} = profLong3Months(idx)';
            profJulDayAggr{iYear-startYear+1} = profJulDay3Months(idx)';
            responseResAggr{iYear-startYear+1} = responseRes3Months(idx)';
        end
        
        nResGrid(iGrid) = sum(cellfun(@length, responseResAggr));
        
        if nResGrid(iGrid) == 0 % No observations in the window
            fprintf('iGrid: %d has no Grid\n', iGrid);
            thetasOpt(iGrid) = NaN;
            thetaLatOpt(iGrid) = NaN;
            thetaLongOpt(iGrid) = NaN;
            thetatOpt(iGrid) = NaN;
            sigmaOpt(iGrid) = NaN;
            nll(iGrid) = NaN;
            
            if isProgressbar
                parfor_progress;
            end
            continue;
        end

        fun = @(params) negLogLikSpaceTime_chol(params,...
                    profLatAggr,profLongAggr,profJulDayAggr,responseResAggr,...
                    kernelType);
        
%        switch kernelType
%           case 'Matern'
        if isStandardize
            logThetasInit = log(5^2);
            logSigmaInit =  log(10);
        else
            if is2step
                logThetasInit = log(10^2);
                logSigmaInit =  log(10);
%{
                switch fluxType
                case 'mass'
                case 'vol'
                    logThetasInit = log(10^2);
                    logSigmaInit =  log(10);
                end
%}
            else
                switch responseTag
                    case 'Dens'
                        logThetasInit = log(400^2);
                        logSigmaInit =  log(100);
                    otherwise
                        logThetasInit = log(400^2);
                        logSigmaInit =  log(100);
                end
            end
        end
        logThetaLatInit = log(5);
        logThetaLongInit = log(5);
        logThetatInit = log(5);
%            case 'ExpGeom'
%                logThetasInit = log(400^2);
%                logThetaLatInit = log(5);
%                logThetaLongInit = log(5);
%                logThetatInit = log(5);
%                logSigmaInit = log(100);
%        end
        
        opts = optimoptions(@fminunc,'Display','notify','Algorithm','quasi-newton','MaxFunctionEvaluations',1500);
        try
            [paramOpt,nll(iGrid),exitFlags(iGrid),~] = fminunc(fun,[logThetasInit, logThetaLatInit, logThetaLongInit, logThetatInit, logSigmaInit],opts);
            thetasOpt(iGrid) = exp(paramOpt(1));
            thetaLatOpt(iGrid) = exp(paramOpt(2));
            thetaLongOpt(iGrid) = exp(paramOpt(3));
            thetatOpt(iGrid) = exp(paramOpt(4));
            sigmaOpt(iGrid) = exp(paramOpt(5));
        catch
            fprintf('Grid Point%d is culprit', iGrid);
            thetasOpt(iGrid) = NaN;
            thetaLatOpt(iGrid) = NaN;
            thetaLongOpt(iGrid) = NaN;
            thetatOpt(iGrid) = NaN;
            sigmaOpt(iGrid) = NaN;
            nll(iGrid) = NaN;
            isFminError(iGrid) = true;
        end

        if isProgressbar
            parfor_progress;
        end
    end
    toc;


    if isProgressbar
        parfor_progress(0);
    end

    if is2step
        saveName = ['./Results/localMLESpaceTime',kernelType,typeTag, fluxType,responseTag,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),EMTag,'.mat']
    else
        saveName = ['./Results/localMLESpaceTime',kernelType,responseTag,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,EMTag,'.mat']
    end
    save(saveName,...
        'latGrid','longGrid','thetasOpt','thetaLatOpt','thetaLongOpt','thetatOpt','sigmaOpt','nll','nResGrid','exitFlags','isFminError','nYear');

end
