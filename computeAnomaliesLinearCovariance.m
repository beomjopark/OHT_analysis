function computeAnomaliesLinearCovariance(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM)

    %% Designed for Single Layer DUACS * Velocity situation
    %% Purpose here is to compute k K^-1 [Var Data] K^-1 k
    %% Difference to computeAnomaliesSeasonSpaceTime is we create covariance matrix from the 

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

    switch numel(month)
      case 12  % Full month provided
          thetasOpt = cell(numel(month), 1);
          thetaLatOpt = cell(numel(month), 1);
          thetaLongOpt = cell(numel(month), 1);
          thetatOpt = cell(numel(month), 1);
          sigmaOpt = cell(numel(month), 1);
          nResGrid = cell(numel(month), 1);
          if is2step
            srcName = ['./Results/localMLESpaceTime',kernelType,typeTag, fluxType,responseTag,verticalSelection,'Season_',num2str(1,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),EMTag,'.mat']
          else
            srcName = ['./Results/localMLESpaceTime',kernelType,responseTag,verticalSelection,'Season_',num2str(1,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,EMTag,'.mat']
          end
          fitMLE = load(srcName);       
          isFminErrorTemp = zeros([numel(month), size(fitMLE.latGrid)]);

          for iMonth = month
              if is2step
                srcName = ['./Results/localMLESpaceTime',kernelType,typeTag, fluxType,responseTag,verticalSelection,'Season_',num2str(iMonth,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),EMTag,'.mat']
              else
                srcName = ['./Results/localMLESpaceTime',kernelType,responseTag,verticalSelection,'Season_',num2str(iMonth,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,EMTag,'.mat']
              end
              fitMLE = load(srcName);             

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
          if is2step
            srcName = ['./Results/localMLESpaceTime',kernelType,typeTag, fluxType,responseTag,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),EMTag,'.mat']
            velResponseTag = 'DUACS';
            srcVelName =  ['./Results/localMLESpaceTime',kernelType,velResponseTag,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w', windowSizeTag,EMTag,'.mat']
          else
            srcName = ['./Results/localMLESpaceTime',kernelType,responseTag,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,EMTag,'.mat']
          end
          fitMLE = load(srcName);
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

          fitMLE1 =  load(srcVelName);
          isFminError = isFminError | (isnan(fitMLE1.thetasOpt) | isnan(fitMLE1.thetaLatOpt) | isnan(fitMLE1.thetaLongOpt) | isnan(fitMLE1.thetatOpt) | isnan(fitMLE1.sigmaOpt));
      otherwise
        error('Not implemented');
    end


    % If rescale needed
    if isStandardize
        if is2step
            res1Dat = load(['./Data/','int',velResponseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,EMTag,'.mat'])
            stdRes1 = res1Dat.stdRes;
            clear res1Dat

            resDat = load(['./Data/',typeTag, fluxType,responseTag,'Res',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),EMTag,'.mat'])
        else
            resDat = load(['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,EMTag,'.mat']);
        end
        stdRes = resDat.stdRes;
        clear resDat
    end

    % Create Destination Folder
    if is2step
        destFolder = ['anomaly_',typeTag,fluxType,responseTag,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeFullTag,'_Eq',num2str(eqBorder),'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag]
    else
        destFolder = ['anomaly_',typeTag,responseTag,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeFullTag,'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag]
    end
    mkdir(['./Results/',destFolder])


    if strcmp(windowType, 'spherical')
        % Determine reference distance
        refDist = distance(0, 180, 0+windowSizeKrig, 180,...
                                referenceEllipsoid('WGS84', 'm'))
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
    Lopts.LT = true;
    Uopts.UT = true;

    for iYear = startYear:endYear
        for iMonth = 1:12

            fprintf("%d / %d\n", iYear, iMonth);

            predYear = iYear;
            predMonth = iMonth;

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

                thetasOpt_f1mon = fitMLE1.thetasOpt;
                thetaLatOpt_f1mon = fitMLE1.thetaLatOpt;
                thetaLongOpt_f1mon = fitMLE1.thetaLongOpt;
                thetatOpt_f1mon = fitMLE1.thetatOpt;
                sigmaOpt_f1mon = fitMLE1.sigmaOpt;                
            end


            nGrid = numel(latGrid);

            midJulDay = datenum(predYear,predMonth,15,12,0,0); % NB: It might be more correct to use 16 insted of 15 here

%            predGrid = zeros(size(latGrid));
            predVarianceGrid = zeros(size(latGrid));
            tic;
            parfor iGrid = 1:nGrid
                
                nResGrid_cur = nResGrid_mon(iGrid);
                thetasOpt_cur = thetasOpt_mon(iGrid);
                thetaLatOpt_cur = thetaLatOpt_mon(iGrid);
                thetaLongOpt_cur = thetaLongOpt_mon(iGrid);
                thetatOpt_cur = thetatOpt_mon(iGrid);
                sigmaOpt_cur = sigmaOpt_mon(iGrid);

                % f1MLE
                thetasOpt_f1cur = thetasOpt_f1mon(iGrid);
                thetaLatOpt_f1cur = thetaLatOpt_f1mon(iGrid);
                thetaLongOpt_f1cur = thetaLongOpt_f1mon(iGrid);
                thetatOpt_f1cur = thetatOpt_f1mon(iGrid);
                sigmaOpt_f1cur = sigmaOpt_f1mon(iGrid);

                if (isnan(thetasOpt_cur) || isnan(thetaLatOpt_cur) || isnan(thetaLongOpt_cur) || isnan(thetatOpt_cur) || isnan(sigmaOpt_cur))
%                  predGrid(iGrid) = NaN;
                  continue;
                end

                predLat = latGrid(iGrid);
                predLong = longGrid(iGrid);

                latMin = predLat - windowSizeMargined;
                latMax = predLat + windowSizeMargined;
                longMin = predLong - windowSizeMargined;
                longMax = predLong + windowSizeMargined;

                if isEqBorder
                    if abs(predLat) <= eqBorder
%                        predGrid(iGrid) = NaN;
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
                targetTempProf3Months = S.targetTempProf3Months;


                idx = find(profLatAggr3Months > latMin & profLatAggr3Months < latMax & profLongAggr3Months > longMin & profLongAggr3Months < longMax);
                if isSpherical
                  is_in_circle = distance(predLat, predLong, profLatAggr3Months(idx), profLongAggr3Months(idx),...
                                    referenceEllipsoid('WGS84', 'm')) < refDist;
                  idx = idx(is_in_circle);
                else
                  idx = find(profLatAggr3Months > latMin & profLatAggr3Months < latMax & profLongAggr3Months > longMin & profLongAggr3Months < longMax);
                end
                
                if isempty(idx) || (nResGrid_cur == 0)
%                    predGrid(iGrid) = NaN;
                    continue;
                end

                profLatAggrSel = profLatAggr3Months(idx);
                profLongAggrSel = profLongAggr3Months(idx);
                profJulDayAggrSel = profJulDayAggr3Months(idx);
                responseResSel = responseRes3Months(idx)';
                targetTempProfSel = targetTempProf3Months(idx);

                nRes = length(responseResSel);

                %% FirstStep  Profile profileVariance: Think of isDeriv = TRUE
                covObs = feval(['spaceTimeCovariance',kernelType,'_vec'],...
                    profLatAggrSel, profLongAggrSel, profJulDayAggrSel,...
                    thetasOpt_f1cur,thetaLatOpt_f1cur,thetaLongOpt_f1cur,thetatOpt_f1cur);
                covGridObs = feval(['spaceTimeCovariance',kernelType,'Deriv'],...
                    profLatAggrSel,profLongAggrSel,profJulDayAggrSel,...
                    profLatAggrSel,profLongAggrSel,profJulDayAggrSel,...
                    thetasOpt_f1cur,thetaLatOpt_f1cur,thetaLongOpt_f1cur,thetatOpt_f1cur,...
                    targetVar);
                covObs(1:nRes+1:end) = covObs(1:nRes+1:end) + sigmaOpt_f1cur.^2;

                % SPD check
                [L, spdFlag] = chol(covObs, 'lower');
                if spdFlag
                    covObs = nearestSPD(covObs);
                    [L, spdFlag] = chol(covObs, 'lower');
                end

                v = linsolve(L, covGridObs', Lopts);
                cov1step = feval(['spaceTimeCovariance',kernelType,'Hess'],...
                    profLatAggrSel,profLongAggrSel,profJulDayAggrSel,...
                    thetasOpt_f1cur,thetaLatOpt_f1cur,thetaLongOpt_f1cur,thetatOpt_f1cur,...
                    targetVar) - (v' * v);

                if isStandardize
                    cov1step = cov1step .* (stdRes1 .^ 2);
                end

                switch typeTag
                    case {'lat', 'ESAlat'}
                        distProf = 1 ./ distance(profLatAggrSel - 0.5, profLongAggrSel, profLatAggrSel + 0.5, profLongAggrSel,...
                                        referenceEllipsoid('WGS84', 'm'));
                    case {'lon', 'ESAlon'}
                        distProf = 1 ./ distance(profLatAggrSel, profLongAggrSel - 0.5, profLatAggrSel, profLongAggrSel + 0.5,...
                                        referenceEllipsoid('WGS84', 'm'));
                end                    

                scaleFactor = targetTempProfSel .* distProf .* gsw_grav(profLatAggrSel) ./ gsw_f(profLatAggrSel) ;
                cov1Obs = (scaleFactor' * scaleFactor) .* cov1step;
                % SPD check
                [U1, spdFlag] = chol(cov1Obs, 'upper');
                if spdFlag
                    cov1Obs = nearestSPD(cov1Obs);
                    [U1, spdFlag] = chol(cov1Obs, 'upper');
                end



                % Propagation step
%{
                if isDeriv
                    covGridObs = feval(['spaceTimeCovariance',kernelType,'Deriv'],...
                        predLat, predLong, midJulDay,...
                        profLatAggrSel,profLongAggrSel,profJulDayAggrSel,...
                        thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur,...
                        targetVar);
                else
%}
                covObs = feval(['spaceTimeCovariance',kernelType,'_vec'],...
                    profLatAggrSel, profLongAggrSel, profJulDayAggrSel,...
                    thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur);
                covGridObs = feval(['spaceTimeCovariance',kernelType],...
                    predLat, predLong, midJulDay,...
                    profLatAggrSel,profLongAggrSel,profJulDayAggrSel,...
                    thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur);
%                end
                covObs(1:nRes+1:end) = covObs(1:nRes+1:end) + sigmaOpt_cur.^2;

                % SPD check
                [~, spdFlag] = chol(covObs, 'lower');
                if spdFlag
                    covObs = nearestSPD(covObs);
                end

%{
                if isDeriv
                    if strcmp(targetVar, 'lat')
                      thetaScale2 = thetaLatOpt_cur.^2;
                    else % lon
                      thetaScale2 = thetaLongOpt_cur.^2;
                    end
%                    predVarianceGrid(iGrid) = 3 * thetasOpt_cur ./ thetaScale2  - covGridObs*linsolve(covObs, covGridObs', symPDopts);
                    v = linsolve(L, covGridObs', Lopts);
                    predVarianceGrid(iGrid) = 3 * thetasOpt_cur ./ thetaScale2  - (v' * v);
                else
%}
                v = U1 * linsolve(covObs, covGridObs', symPDopts);
                predVarianceGrid(iGrid) = v' * v;
  %              end
            end
            toc;

            if isStandardize
                predVarianceGrid = predVarianceGrid .* (stdRes^2);
            end
            if isDeriv
%{
                saveName = ['./Results/',destFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,targetVar,'Deriv_',...
                    num2str(predMonth,'%02d'),'_',num2str(predYear),'.mat'];
%}
            else
                if is2step
                    saveName = ['./Results/',destFolder,'/anomaly',typeTag,fluxType,responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'LinCov_',...
                        num2str(predMonth,'%02d'),'_',num2str(predYear),'.mat'];
                else
                    saveName = ['./Results/',destFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'LinCov_',...
                            num2str(predMonth,'%02d'),'_',num2str(predYear),'.mat'];
                end
            end
            save(saveName,'predVarianceGrid','latGrid','longGrid','midJulDay');

        end 
    end

end