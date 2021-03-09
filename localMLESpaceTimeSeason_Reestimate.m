function localMLESpaceTimeSeason_Reestimate(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isProgressbar, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth)
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

    if is2step
        saveName = ['./Results/localMLESpaceTime',kernelType,typeTag, fluxType,responseTag,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),EMTag,'.mat']
    else
        saveName = ['./Results/localMLESpaceTime',kernelType,responseTag,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,EMTag,'.mat']
    end
    S = load(saveName, ...
            'latGrid','longGrid','thetasOpt','thetaLatOpt','thetaLongOpt','thetatOpt','sigmaOpt','nll','nResGrid','exitFlags','isFminError','nYear')
    
    latGrid = S.latGrid;
    longGrid = S.longGrid;
    nGrid = numel(latGrid);

    thetasOpt = S.thetasOpt;
    thetaLatOpt = S.thetaLatOpt;
    thetaLongOpt = S.thetaLongOpt;
    thetatOpt = S.thetatOpt;
    sigmaOpt = S.sigmaOpt;
    nll = S.nll;
    nResGrid = S.nResGrid;
    exitFlags = S.exitFlags;
    isFminError = S.isFminError;

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
        dat = load(['./Data/',typeTag,'TempDens','Prof','PchipPotTemp',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);
        intStart = dat.intStart;

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
        
    
    fprintf('Reestimate %d observations from parfor error\n', sum(isFminError(:)))

    % Spurious check
    reGrid = find(mask.*thetatOpt > 365)';
    fprintf('Reestimate %d observations of spurious error\n', numel(reGrid));
    reGrid = union(find(isFminError)', reGrid);

    % Discard previous itearIdx, if it exists
    % fileName = ['iterIdxLocalMLESpaceTime',kernelType,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(windowSize),'.txt'];
    % fileID = fopen(fileName,'w');
    % fclose(fileID);

    if isProgressbar
        parfor_progress(nGrid);
    end

    tic;
    for iGrid = reGrid%find(isFminError) %26605:26610%1:nGrid
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

        profLatAggr = cell(1,nYear);
        profLongAggr = cell(1,nYear);
        profJulDayAggr = cell(1,nYear);
        responseResAggr = cell(1,nYear);

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
        
        for iYear = startYear:endYear
            S = load(['./Data/Extended/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(month,'%02d'),'_',num2str(iYear),'_extended.mat']);
            
            profLat3Months = S.profLatAggr3Months;
            profLong3Months = S.profLongAggr3Months;
            profJulDay3Months = S.profJulDayAggr3Months;
            switch typeTag
                case 'int'
                    if strcmp(responseTag, 'Temp')
                        responseRes3Months = S.targetTempRes3Months;
                    else
                        responseRes3Months = S.intDensRes3Months;
                    end
                case {'intlat', 'intlon'}
                    responseRes3Months = S.intFluxRes3Months;
                otherwise
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
        
        nResGrid(iGrid) = sum(cellfun(@length,responseResAggr));
        
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

        fun = @(params) negLogLikSpaceTime(params,...
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
                    logThetasInit = log(10^2);
                    logSigmaInit =  log(10);
                case 'vol'
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
        if isFminError(iGrid)
            try
                [paramOpt,nll(iGrid),exitFlags(iGrid),~] = fminunc(fun,[logThetasInit, logThetaLatInit, logThetaLongInit, logThetatInit, logSigmaInit],opts);
                fprintf('Grid Point %d reoptimized\n', iGrid);

                thetasOpt(iGrid) = exp(paramOpt(1));
                thetaLatOpt(iGrid) = exp(paramOpt(2));
                thetaLongOpt(iGrid) = exp(paramOpt(3));
                thetatOpt(iGrid) = exp(paramOpt(4));
                sigmaOpt(iGrid) = exp(paramOpt(5));
                isFminError(iGrid) = false;
            catch
                fprintf('Grid Point%d is still problematic\n', iGrid);

                thetasOpt(iGrid) = NaN;
                thetaLatOpt(iGrid) = NaN;
                thetaLongOpt(iGrid) = NaN;
                thetatOpt(iGrid) = NaN;
                sigmaOpt(iGrid) = NaN;
                nll(iGrid) = NaN;
                isFminError(iGrid) = true;
            end
        else
            % Spurious estimation cases
            try
                [paramOptCon,nllCon,exitFlagsCon,~] = fmincon(fun,[logThetasInit, logThetaLatInit, logThetaLongInit, logThetatInit, logSigmaInit],...
                    [],[],[],[],[],[Inf;Inf;Inf;log(365*nYear);Inf],[],...
                    optimoptions(@fmincon,'Display','notify'));
                if nllCon < nll(iGrid)
                    fprintf('Grid Point %d reoptimized with constraint\n', iGrid);
                    paramOpt = paramOptCon;
                    nll(iGrid) = nllCon;
                    exitFlags(iGrid) = exitFlagsCon;
                end
                thetasOpt(iGrid) = exp(paramOpt(1));
                thetaLatOpt(iGrid) = exp(paramOpt(2));
                thetaLongOpt(iGrid) = exp(paramOpt(3));
                thetatOpt(iGrid) = exp(paramOpt(4));
                sigmaOpt(iGrid) = exp(paramOpt(5));
            catch
                fprintf('Grid Point %d is not reoptimized\n', iGrid);                
            end
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
