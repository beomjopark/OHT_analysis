%% Finding the convolutional kernel
function computeAnomaliesSeasonSpaceTime_Kernel(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM)

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
          windowTypeTag = 'spherical'
          windowSizeMargined = windowSizeKrig * 2;
        otherwise
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

    % Create Destination Folder
    if is2step
        destFolder = ['anomaly_',typeTag,fluxType,responseTag,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeFullTag,'_Eq',num2str(eqBorder),'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag];
    else
        destFolder = ['anomaly_',typeTag,responseTag,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeFullTag,'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag,'Kernel'];
    end
    mkdir(['./Results/',destFolder])


    if strcmp(windowType, 'spherical')
        % Determine reference distance
        refDist = distance(0, 180, 0+windowSizeKrig, 180,...
                                referenceEllipsoid('GRS80', 'm'))
    end


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
            end


            nGrid = numel(latGrid);

            midJulDay = datenum(predYear,predMonth,15,12,0,0); % NB: It might be more correct to use 16 insted of 15 here

            % Ficticious Window
            nInjectDay = 10;
            latGridFlat = latGrid(:);
            longGridFlat = longGrid(:);

%            predVarianceGrid = zeros(size(latGrid));
            tic;
%            parfor iGrid = 1:nGrid
            for iGrid = int64(linspace(1, nGrid, 100))%1:nGrid
                predGrid = sparse(size(latGrid,1), size(latGrid,2));

                nResGrid_cur = nResGrid_mon(iGrid);
                thetasOpt_cur = thetasOpt_mon(iGrid);
                thetaLatOpt_cur = thetaLatOpt_mon(iGrid);
                thetaLongOpt_cur = thetaLongOpt_mon(iGrid);
                thetatOpt_cur = thetatOpt_mon(iGrid);
                sigmaOpt_cur = sigmaOpt_mon(iGrid);

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
                switch typeTag
                    case 'int'
                        responseRes3Months = S.intDensRes3Months;
                    case 'lat'
                        responseRes3Months = S.latFluxRes3Months;
                    case 'lon'
                        responseRes3Months = S.lonFluxRes3Months;
                    otherwise %intlatlon
                        responseRes3Months = S.intFluxRes3Months;
                end

                idx = find(profLatAggr3Months > latMin & profLatAggr3Months < latMax & profLongAggr3Months > longMin & profLongAggr3Months < longMax);
                switch windowType
                    case 'spherical'
                      is_in_circle = distance(predLat, predLong, profLatAggr3Months(idx), profLongAggr3Months(idx),...
                                        referenceEllipsoid('GRS80', 'm')) < refDist;
                      idx = idx(is_in_circle);
                    case 'box_var'
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

                %% Add ficticious centerpoint
                predStartYear = predYear;
                predEndYear = predYear;
                predStartMonth = predMonth - 1;
                predEndMonth = predMonth + 1;
            
                if predYear == startYear && predMonth == 1
                    predStartMonth = predMonth;
                elseif predYear == endYear && predMonth == 12
                    predEndMonth = predMonth;
                elseif predMonth == 1
                    predStartYear = predYear - 1;
                    predStartMonth = 12;
                elseif predMonth == 12
                    predEndYear = predYear + 1;
                    predEndMonth = 1;
                else
                    % Do nothing
                end
                
                injectDay = linspace(datenum(predStartYear,predStartMonth,1,0,0,0),...
                    datenum(predEndYear,predEndMonth,eomday(predEndYear,predEndMonth),23,59,0),...
                    nInjectDay);
                profLatAggrSel = [profLatAggrSel predLat.*ones(1, nInjectDay)];
                profLongAggrSel = [profLongAggrSel predLong.*ones(1, nInjectDay)];
                profJulDayAggrSel = [profJulDayAggrSel injectDay];
                responseResSel = [zeros(size(responseResSel)); ones(nInjectDay, 1)];

                nRes = length(responseResSel);


                %% Indexing for grid
                gridIdx = find(latGridFlat > latMin & latGridFlat < latMax & longGridFlat > longMin & longGridFlat < longMax);
                switch windowType
                    case 'spherical'
                      is_in_circle = distance(predLat, predLong, latGridFlat(gridIdx), longGridFlat(gridIdx),...
                                        referenceEllipsoid('GRS80', 'm')) < refDist;
                      gridIdx = gridIdx(is_in_circle);
                    case 'box_var'
                        gridIdx = find(latGridFlat > latMin & latGridFlat < latMax & longGridFlat > longMin & longGridFlat < longMax);
                end
                
                if isDeriv
                    covObs = feval(['spaceTimeCovariance',kernelType,'_vec'],...
                        profLatAggrSel, profLongAggrSel, profJulDayAggrSel,...
                        thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur);
                    covGridObs = feval(['spaceTimeCovariance',kernelType,'Deriv'],...
                        latGridFlat(gridIdx), longGridFlat(gridIdx), midJulDay,...
                        profLatAggrSel,profLongAggrSel,profJulDayAggrSel,...
                        thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur,...
                        targetVar);
                else
                    covObs = feval(['spaceTimeCovariance',kernelType,'_vec'],...
                        profLatAggrSel, profLongAggrSel, profJulDayAggrSel,...
                        thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur);
                    covGridObs = feval(['spaceTimeCovariance',kernelType],...
                        latGridFlat(gridIdx), longGridFlat(gridIdx), midJulDay,...
                        profLatAggrSel,profLongAggrSel,profJulDayAggrSel,...
                        thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur);
                end
                
                predGrid(gridIdx) = covGridObs * ((covObs + sigmaOpt_cur.^2 * eye(nRes)) \ responseResSel);

                if isStandardize
                    predGrid = predGrid .* stdRes;
                    %predVarianceGrid = predVarianceGrid .* (stdRes^2);
                end
                
                if isDeriv
                    parsave(['./Results/',destFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,targetVar,'Deriv_',...
                    num2str(predMonth,'%02d'),'_',num2str(predYear),'_',num2str(iGrid),'.mat'], iGrid, predGrid);
                else
                    if is2step
                        parsave(['./Results/',destFolder,'/anomaly',typeTag,fluxType,responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'_',...
                            num2str(predMonth,'%02d'),'_',num2str(predYear),'_',num2str(iGrid),'.mat'], iGrid, predGrid);                    
                    else
                        parsave(['./Results/',destFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'_',...
                                num2str(predMonth,'%02d'),'_',num2str(predYear),'_',num2str(iGrid),'.mat'], iGrid, predGrid);
                    end
                end

                
                
                % Outlier clearance
%                if abs(predGrid(iGrid)) > 10
%                    [dup_row, dup_col] = find(triu(squareform(pdist(covObs)) < 1e-05 - eye(nRes)));
%                    idxSet = [];
%                    excludeSet = [];
%                    
%                    idxSet = dup_row(1);
%                    excludeSet = dup_col(1);
%                    for iRow = 1:nRes
%                        if any()
%                        idxSet = [idxSet, dup_row(iRow)]
%                    end
%                    unique(dup_row)
%            
%                    covObsRed = covObs([1:4, 6:7, 10:end], [1:4, 6:7, 10:end]);
%                    nResRed = size(covObsRed, 1);
%                    predGrid(iGrid) = covGridObs([1:4, 6:7, 10:end]) * ((covObsRed + sigmaOpt(iGrid)^2*eye(nResRed)) \ responseResSel([1:4, 6:7, 10:end]));
%                end
%{                
                if isDeriv
                    switch targetVar
                        case 'lat'
                            thetaScale2 = thetaLatOpt_cur.^2;
                        otherwise
                            thetaScale2 = thetaLongOpt_cur.^2;
                    end
                    predVarianceGrid(iGrid) = 3 * thetasOpt_cur ./ thetaScale2  - covGridObs*((covObs + sigmaOpt_cur^2*eye(nRes))\(covGridObs'));
                else
                    predVarianceGrid(iGrid) = thetasOpt_cur + sigmaOpt_cur^2 - covGridObs*((covObs + sigmaOpt_cur^2*eye(nRes))\(covGridObs'));
                end
%}
            end
            toc;
            
            %{
            if isStandardize
                predGrid = predGrid .* stdRes;
                predVarianceGrid = predVarianceGrid .* (stdRes^2);
            end
            if isDeriv
                save(['./Results/',destFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,targetVar,'Deriv_',...
                    num2str(predMonth,'%02d'),'_',num2str(predYear),'.mat'],'predGrid','predVarianceGrid','latGrid','longGrid','midJulDay');
            else
                if is2step
                    save(['./Results/',destFolder,'/anomaly',typeTag,fluxType,responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'_',...
                        num2str(predMonth,'%02d'),'_',num2str(predYear),'.mat'],'predGrid','predVarianceGrid','latGrid','longGrid','midJulDay');                    
                else
                    save(['./Results/',destFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'_',...
                            num2str(predMonth,'%02d'),'_',num2str(predYear),'.mat'],'predGrid','predVarianceGrid','latGrid','longGrid','midJulDay');
                end
            end
            %}
        end 
    end

end