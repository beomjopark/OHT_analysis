function computeAnomaliesSeasonSpaceTime(kernelType, month, typeTag, responseTag, verticalSelection, targetPres, dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM)

    if nargin < 12 || isempty(isStandardize)
        isStandardize = false;
    end
    if nargin < 14 || isempty(eqBorder)
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
            isSpherical = true;
            windowTypeTag = 'spherical'
            windowSizeMargined = windowSizeKrig * 2;
        otherwise
            isSpherical = false;
            windowTypeTag = []
            windowSizeMargined = windowSizeKrig;
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

    if is2step
        fitMLE = load(['./Results/localMLESpaceTime',kernelType,typeTag, fluxType,responseTag,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',num2str(windowSize),'_Eq',num2str(eqBorder),EMTag,'.mat']);
    else
        switch numel(month)
          case 12  % Full month provided
              thetasOpt = cell(numel(month), 1);
              thetaLatOpt = cell(numel(month), 1);
              thetaLongOpt = cell(numel(month), 1);
              thetatOpt = cell(numel(month), 1);
              sigmaOpt = cell(numel(month), 1);
              nResGrid = cell(numel(month), 1);

              fitMLE = load(['./Results/localMLESpaceTime',kernelType,responseTag,verticalSelection,'Season_',num2str(1,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,EMTag,'.mat']);
              isFminErrorTemp = zeros([numel(month), size(fitMLE.latGrid)]);

              for iMonth = month
                  fitMLE = load(['./Results/localMLESpaceTime',kernelType,responseTag,verticalSelection,'Season_',num2str(iMonth,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,EMTag,'.mat']);
                  latGrid = fitMLE.latGrid;
                  longGrid = fitMLE.longGrid;
                  thetasOpt{iMonth} = fitMLE.thetasOpt;
                  thetaLatOpt{iMonth} = fitMLE.thetaLatOpt;
                  thetaLongOpt{iMonth} = fitMLE.thetaLongOpt;
                  thetatOpt{iMonth} = fitMLE.thetatOpt;
                  sigmaOpt{iMonth} = fitMLE.sigmaOpt;
                  nResGrid{iMonth} = fitMLE.nResGrid;
                  isFminErrorTemp(iMonth,:,:) = (isnan(thetasOpt{iMonth}) | isnan(thetaLatOpt{iMonth}) | isnan(thetaLongOpt{iMonth}) | isnan(thetatOpt{iMonth}) | isnan(sigmaOpt{iMonth}));
              end
              isFminError = squeeze(any(isFminErrorTemp, 1));
              clear isFminErrorTemp
              clear fitMLE
          case 1
              fitMLE = load(['./Results/localMLESpaceTime',kernelType,responseTag,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,EMTag,'.mat']);
              latGrid = fitMLE.latGrid;
              longGrid = fitMLE.longGrid;
              thetasOpt = fitMLE.thetasOpt;
              thetaLatOpt = fitMLE.thetaLatOpt;
              thetaLongOpt = fitMLE.thetaLongOpt;
              thetatOpt = fitMLE.thetatOpt;
              sigmaOpt = fitMLE.sigmaOpt;
              nResGrid = fitMLE.nResGrid;
%              nll = fitMLE.nll;
              isFminError = (isnan(thetasOpt) | isnan(thetaLatOpt) | isnan(thetaLongOpt) | isnan(thetatOpt) | isnan(sigmaOpt));
              clear fitMLE
          otherwise
            error('Not implemented');
        end
    end


    % If rescale needed
    if isStandardize
        if is2step
            resDat = load(['./Data/',typeTag, fluxType,responseTag,'Res',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),EMTag,'.mat']);
        else
            if isAdjusted
                filterTag = ['Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag];%[];
            else
                filterTag = ['Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag];
            end
            resDat = load(['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,filterTag,EMTag,'.mat']);
        end
        stdRes = resDat.stdRes;
        clear resDat
    end

    % Compute the anomalies for prediction points in target
    nTargetPres = numel(targetPres);
    if nTargetPres > 1
        target = load(['./Data/targetTempProf','PchipPotTemp',num2str(min(targetPres)),'_',num2str(max(targetPres)),dataYear,'.mat'], 'profLongAggrSel', 'profLatAggrSel','profJulDayAggrSel');
    else
        target = load(['./Data/targetTempProf','PchipPotTemp',num2str(targetPres),dataYear,'.mat'], 'profLongAggrSel', 'profLatAggrSel','profJulDayAggrSel');    
    end

    %% Compute Anomalies at the rounded grid points
    nProf = length(target.profLongAggrSel);
    predRes = zeros(1, nProf);
    predVarianceRes = zeros(1, nProf);

    profLatAggrSelRounded = roundHalf(target.profLatAggrSel);
    profLongAggrSelRounded = roundHalf(target.profLongAggrSel);
    profJulDayAggrSel = target.profJulDayAggrSel;
    
    temp = datevec(profJulDayAggrSel);
    profYearAggrSel = temp(:,1)';
    profMonthAggrSel = temp(:,2)';

    clear temp;
    clear target;

    % Defined for index Search
    latLinspace = latGrid(1,:);
    longLinspace = longGrid(:,1);
    sizeGrid = size(latGrid);

    % Create Destination Folder
    if is2step
        destFolder = ['anomaly_',typeTag,fluxType,responseTag,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeFullTag,'_Eq',num2str(eqBorder),'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag];
    else
        destFolder = ['anomaly_',typeTag,responseTag,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeFullTag,'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag];
    end
    mkdir(['./Results/',destFolder])


    if strcmp(windowType, 'spherical')
        % Determine reference distance
        refDist = distance(0, 180, 0+windowSizeKrig, 180,...
                                referenceEllipsoid('GRS80', 'm'))
    end

    % responseName
    switch typeTag
        case 'int'
          if strcmp(responseTag, 'Temp')
              responseName = 'targetTempRes3Months';
%              responseRes3Months = S.targetTempRes3Months;
          else
              responseName = 'intDensRes3Months';            
%              responseRes3Months = S.intDensRes3Months;
          end
        case {'intlat', 'intlon'}
            responseName = 'intFluxRes3Months';
%            responseRes3Months = S.intFluxRes3Months;
        otherwise %latlon
            responseName = 'FluxRes3Months';
%            responseRes3Months = S.FluxRes3Months;
    end

    symPDopts.POSDEF = true;
    symPDopts.SYM = true;

    for iYear = startYear:endYear
        for iMonth = 1:12
            fprintf("%d / %d\n", iYear, iMonth);
            if isFullMonth
                nResGrid_mon = nResGrid{iMonth};
                thetasOpt_mon = thetasOpt{iMonth};
                thetaLatOpt_mon = thetaLatOpt{iMonth};
                thetaLongOpt_mon = thetaLongOpt{iMonth};
                thetatOpt_mon = thetatOpt{iMonth};
                sigmaOpt_mon = sigmaOpt{iMonth};
            else
                nResGrid_mon = nResGrid;
                thetasOpt_mon = thetasOpt;
                thetaLatOpt_mon = thetaLatOpt;
                thetaLongOpt_mon = thetaLongOpt;
                thetatOpt_mon = thetatOpt;
                sigmaOpt_mon = sigmaOpt;
            end

            % Filter Relevant Year Month Profile
            isTime = ((profYearAggrSel == iYear) & (profMonthAggrSel == iMonth));
            nGrid = sum(isTime);
            if ~nGrid
              continue;
            end

            profLatAggrSelCur = profLatAggrSelRounded(isTime);
            profLongAggrSelCur = profLongAggrSelRounded(isTime);
            profJulDayAggrSelCur = profJulDayAggrSel(isTime);

            tic;
            predResTemp = zeros(1, nGrid);
            predVarResTemp = zeros(1, nGrid);
            parfor iProf = 1:nGrid

                predLat = profLatAggrSelCur(iProf);  %latGrid(iProf);
                predLong = profLongAggrSelCur(iProf);  %longGrid(iProf);
                predJulDay = profJulDayAggrSelCur(iProf);

                % Find matching Grid Index
                iLat = find(latLinspace == predLat);
                iLong = find(longLinspace == predLong);
                iGrid = sub2ind(sizeGrid, iLong, iLat);
                
                nResGrid_cur = nResGrid_mon(iGrid);
                thetasOpt_cur = thetasOpt_mon(iGrid);
                thetaLatOpt_cur = thetaLatOpt_mon(iGrid);
                thetaLongOpt_cur = thetaLongOpt_mon(iGrid);
                thetatOpt_cur = thetatOpt_mon(iGrid);
                sigmaOpt_cur = sigmaOpt_mon(iGrid);

                latMin = predLat - windowSizeMargined;
                latMax = predLat + windowSizeMargined;
                longMin = predLong - windowSizeMargined;
                longMax = predLong + windowSizeMargined;

                if isEqBorder
                    if abs(predLat) <= eqBorder
                        predResTemp(iProf) = NaN;
                        continue;
                    end
                    if predLat > 0
                        latMin = max([latMin, eqBorder]);
                    else
                        latMax = min([latMax, - eqBorder]);
                    end
                end

                S = load(['./Data/Extended/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(iMonth,'%02d'),'_',num2str(iYear),'_extended.mat']);

                profLatAggr3Months = S.profLatAggr3Months;
                profLongAggr3Months = S.profLongAggr3Months;
                profJulDayAggr3Months = S.profJulDayAggr3Months;
                responseRes3Months = S.(responseName);


                idx = find(profLatAggr3Months > latMin & profLatAggr3Months < latMax & profLongAggr3Months > longMin & profLongAggr3Months < longMax);
                if isSpherical
                      is_in_circle = distance(predLat, predLong, profLatAggr3Months(idx), profLongAggr3Months(idx),...
                                        referenceEllipsoid('GRS80', 'm')) < refDist;
                      idx = idx(is_in_circle);
                else
                    idx = find(profLatAggr3Months > latMin & profLatAggr3Months < latMax & profLongAggr3Months > longMin & profLongAggr3Months < longMax);
                end
                
                if isempty(idx) || (nResGrid_cur == 0)
                    predResTemp(iProf) = NaN;
                    continue;
                end

                profLatAggrSel = profLatAggr3Months(idx);
                profLongAggrSel = profLongAggr3Months(idx);
                profJulDayAggrSel = profJulDayAggr3Months(idx);
                responseResSel = responseRes3Months(idx)';

                nRes = length(responseResSel);
                
                if isDeriv
                    covObs = feval(['spaceTimeCovariance',kernelType,'_vec'],...
                        profLatAggrSel, profLongAggrSel, profJulDayAggrSel,...
                        thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur);
                    covGridObs = feval(['spaceTimeCovariance',kernelType,'Deriv'],...
                        predLat, predLong, predJulDay,...
                        profLatAggrSel,profLongAggrSel,profJulDayAggrSel,...
                        thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur,...
                        targetVar);
                else
                    covObs = feval(['spaceTimeCovariance',kernelType,'_vec'],...
                        profLatAggrSel, profLongAggrSel, profJulDayAggrSel,...
                        thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur);
                    covGridObs = feval(['spaceTimeCovariance',kernelType],...
                        predLat, predLong, predJulDay,...
                        profLatAggrSel,profLongAggrSel,profJulDayAggrSel,...
                        thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur);            
                end
                
                covObs = covObs + sigmaOpt_cur.^2 * eye(nRes);
                try
                    predResTemp(iProf) = covGridObs * linsolve(covObs, responseResSel, symPDopts);
                    if isDeriv
                        switch targetVar
                            case 'lat'
                                thetaScale2 = thetaLatOpt_cur.^2;
                            otherwise
                                thetaScale2 = thetaLongOpt_cur.^2;
                        end
                        predVarResTemp(iProf) = 3 * thetasOpt_cur ./ thetaScale2  - covGridObs*linsolve(covObs, covGridObs', symPDopts);
                    else
                        predVarResTemp(iProf) = thetasOpt_cur + sigmaOpt_cur^2 - covGridObs*linsolve(covObs, covGridObs', symPDopts);
                    end
                catch
                    predResTemp(iProf) = NaN;
                    predVarResTemp(iProf) = NaN;
                end
            end
            toc;

            predRes(isTime) = predResTemp;
            predVarianceRes(isTime) = predVarResTemp;
        end 
    end

    if isStandardize
        predRes = predRes .* stdRes;
        predVarianceRes = predVarianceRes .* (stdRes^2);
    end

    if isDeriv
        save(['./Results/',destFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,targetVar,'Deriv_Profile','.mat'],...
          'predRes','predVarianceRes','profLatAggrSelRounded','profLongAggrSelRounded','profJulDayAggrSel');
    else
        if is2step
            save(['./Results/',destFolder,'/anomaly',typeTag,fluxType,responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'_Profile','.mat'],...
              'predRes','predVarianceRes','profLatAggrSelRounded','profLongAggrSelRounded','profJulDayAggrSel');                    
        else
            save(['./Results/',destFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'_Profile','.mat'],...
              'predRes','predVarianceRes','profLatAggrSelRounded','profLongAggrSelRounded','profJulDayAggrSel');
        end
    end

end