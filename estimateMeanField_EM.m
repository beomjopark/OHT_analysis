function estimateMeanField_EM(kernelType, month, meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs,...
                           is2step, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM)
%% Estimate MeanField with local regression.
%% Choice of function set by meanTag
  if isempty(eqBorder)
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
      windowSizeMargined = windowSizeMean * 2;
    otherwise
      windowTypeTag = []
      windowSizeMargined = windowSizeMean;
  end


  % Additional Parameter
  tag = 'PchipPotTemp';
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
    prevEMTag = []
  elseif iterEM == 1
    EMTag = ['EM',num2str(iterEM)]
    prevEMTag = []
  else
    EMTag = ['EM',num2str(iterEM)]
    prevEMTag = ['EM',num2str(iterEM-1)]
  end


  isFullMonth = (numel(month) == 12)
  
  if isFullMonth && ~isempty(EMTag)
    EMTag = ['Full', EMTag]
  end
  if isFullMonth && ~isempty(prevEMTag)
    prevEMTag = ['Full', prevEMTag]
  end

  if strcmp(dataYear, '_2007_2018')
      midJulDay = datenum(2013,01,01,12,0,0); %datestr(median([min(profJulDayAggrSel) max(profJulDayAggrSel)]))
  else
      midJulDay = datenum(2012,01,01,12,0,0); %734869.5;
  end

  % Load Profile and dataMask
  if is2step
    if isnumeric(verticalSelection) % intlat/intlon
        targetPres = verticalSelection;
        presString = [num2str(min(targetPres)),'_',num2str(max(targetPres))];
        srcName = ['./Data/',typeTag,fluxType,responseTag,'Prof',tag,presString,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'.mat']
        load(srcName);

        maskName = ['./Data/dataMask',typeTag,responseTag,presString,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'.mat']
        load(maskName);
    else
        % For each target Pressure
        srcName = ['./Data/',typeTag,fluxType,responseTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'.mat']
        load(srcName);
        
        maskName = ['./Data/dataMask',typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'.mat']
        load(maskName);
    end
  else
      if strcmp(responseTag, 'Sal')
          load(['./Data/',typeTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'.mat'])
      elseif strcmp(responseTag, 'DUACS') || strcmp(responseTag, 'ESA') || strcmp(responseTag, 'SOSITemp') % ESA TEMP
          load(['./Data/',typeTag,responseTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'.mat']);
      else % 'Temp', 'Dens'
          load(['./Data/',typeTag,'TempDens','Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'.mat']);
      end
      load(['./Data/dataMask',verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'.mat']);    
      %{
      if numel(intStart) > 1
        load(['./Data/dataMask','Relative',num2str(min(intStart)),'_',num2str(max(intStart)),dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);    
      else
      %}
%      end
  end

  [latGrid,longGrid] = meshgrid(linspace(-89.5,89.5,180), linspace(20.5,379.5,360));
  nGrid = numel(latGrid);

  % Load MLE
  yearRange = strsplit(dataYear, '_');
  startYear = str2num(yearRange{2});
  endYear = str2num(yearRange{3});
  if isempty(iterEM) || (iterEM == 0)
    thetasOpt = NaN(size(latGrid));
    thetaLatOpt = NaN(size(latGrid));
    thetaLongOpt = NaN(size(latGrid));
    thetatOpt = NaN(size(latGrid));
    sigmaOpt = NaN(size(latGrid));
    isFminError = zeros(size(latGrid), "logical");
  else
    switch numel(month)
      case 12  % Full month provided
          thetasOpt = cell(numel(month), 1);
          thetaLatOpt = cell(numel(month), 1);
          thetaLongOpt = cell(numel(month), 1);
          thetatOpt = cell(numel(month), 1);
          sigmaOpt = cell(numel(month), 1);
          isFminErrorTemp = zeros([numel(month), size(latGrid)]);
          for iMonth = month
              if is2step
                srcName = ['./Results/localMLESpaceTime',kernelType,typeTag, fluxType,responseTag,verticalSelection,'Season_',num2str(iMonth,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),prevEMTag,'.mat']
              else
                srcName = ['./Results/localMLESpaceTime',kernelType,responseTag,verticalSelection,'Season_',num2str(iMonth,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,prevEMTag,'.mat']
              end
              fitMLE = load(srcName);
              latGrid = fitMLE.latGrid;
              longGrid = fitMLE.longGrid;
              thetasOpt{iMonth} = fitMLE.thetasOpt;
              thetaLatOpt{iMonth} = fitMLE.thetaLatOpt;
              thetaLongOpt{iMonth} = fitMLE.thetaLongOpt;
              thetatOpt{iMonth} = fitMLE.thetatOpt;
              sigmaOpt{iMonth} = fitMLE.sigmaOpt;
              isFminErrorTemp(iMonth,:,:) = (isnan(thetasOpt{iMonth}) | isnan(thetaLatOpt{iMonth}) | isnan(thetaLongOpt{iMonth}) | isnan(thetatOpt{iMonth}) | isnan(sigmaOpt{iMonth}));
          end
          isFminError = squeeze(any(isFminErrorTemp, 1));
          clear isFminErrorTemp
      case 1
          if is2step
            if isnumeric(verticalSelection)
              srcName = ['./Results/localMLESpaceTime',kernelType,typeTag, fluxType,responseTag,presString,'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),prevEMTag,'.mat']
            else
              srcName = ['./Results/localMLESpaceTime',kernelType,typeTag, fluxType,responseTag,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),prevEMTag,'.mat']
            end
          else
            srcName = ['./Results/localMLESpaceTime',kernelType,responseTag,verticalSelection,'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,prevEMTag,'.mat']
          end
          fitMLE = load(srcName);
          latGrid = fitMLE.latGrid;
          longGrid = fitMLE.longGrid;
          thetasOpt = fitMLE.thetasOpt;
          thetaLatOpt = fitMLE.thetaLatOpt;
          thetaLongOpt = fitMLE.thetaLongOpt;
          thetatOpt = fitMLE.thetatOpt;
          sigmaOpt = fitMLE.sigmaOpt;
%{
          nll = fitMLE.nll;
          nResGrid = fitMLE.nResGrid;
%}
          isFminError = (isnan(thetasOpt) | isnan(thetaLatOpt) | isnan(thetaLongOpt) | isnan(thetatOpt) | isnan(sigmaOpt));
      otherwise
        error('Not implemented');
    end
  end

  % Setup dataMask
  maskJohn = ncread('./RG_climatology/RG_ArgoClim_Temperature_2016.nc','BATHYMETRY_MASK',[1 1 25],[Inf Inf 1]);
  maskJohn(maskJohn == 0) = 1;
  maskJohn = [NaN*ones(360,25) maskJohn NaN*ones(360,25)];
  mask = maskJohn .* dataMask;

  %% Enable wrap around by duplicating boundary data
  leftBoundaryIdx = find(profLongAggrSel <= 20 + windowSizeMargined);
  rightBoundaryIdx = find(profLongAggrSel >= 380 - windowSizeMargined);
  profLongAggrSel = [profLongAggrSel profLongAggrSel(leftBoundaryIdx) + 360 profLongAggrSel(rightBoundaryIdx) - 360];
  profLatAggrSel = [profLatAggrSel profLatAggrSel(leftBoundaryIdx) profLatAggrSel(rightBoundaryIdx)];
  profJulDayAggrSel = [profJulDayAggrSel profJulDayAggrSel(leftBoundaryIdx) profJulDayAggrSel(rightBoundaryIdx)];

  parsedJulDay = datevec(profJulDayAggrSel);
  profYearAggrSel = parsedJulDay(:,1);
  profMonthAggrSel = parsedJulDay(:,2);

  if is2step
      switch typeTag
          case {'intlat', 'intlon'}
              responseProf = [intFluxProf intFluxProf(leftBoundaryIdx) intFluxProf(rightBoundaryIdx)];
          otherwise %{'lat', 'lon'}
              responseProf = [profFluxAggrSel profFluxAggrSel(leftBoundaryIdx) profFluxAggrSel(rightBoundaryIdx)];
      end
  else
      switch responseTag
          case 'Temp'
              responseProf = [intTempProf intTempProf(leftBoundaryIdx) intTempProf(rightBoundaryIdx)];
          case 'Dens'
              responseProf = [targetDynhProf targetDynhProf(leftBoundaryIdx) targetDynhProf(rightBoundaryIdx)];
%              responseProf = [intDensProf intDensProf(leftBoundaryIdx) intDensProf(rightBoundaryIdx)];
%              responseProf = [intDensRes intDensRes(leftBoundaryIdx) intDensRes(rightBoundaryIdx)];
          case 'DUACS'
              responseProf = [targetADTProf targetADTProf(leftBoundaryIdx) targetADTProf(rightBoundaryIdx)];
          case {'ESA', 'SOSITemp'}
              responseProf = [targetSSTProf targetSSTProf(leftBoundaryIdx) targetSSTProf(rightBoundaryIdx)];
          case 'Sal'
              responseProf = [targetSalProfPchip targetSalProfPchip(leftBoundaryIdx) targetSalProfPchip(rightBoundaryIdx)];
      end
  end

  %% Calculate mean field using a moving window

  betaGrid = zeros([nGrid, 18]);
  switch meanTag 
      case 'NoTrend'
          betaLength = 18; %[betaWindow; 0; 0]; % Put zero coefficients in place of the time trend terms
          betaGrid = zeros([nGrid, betaLength]);
      case 'NoTrendVelSeas3'
          betaLength = 24;
          betaGrid = zeros([nGrid, betaLength]);
      case 'TrendDens'
          betaLength = 21;
          betaGrid = zeros([nGrid, betaLength]);
      case 'TrendDens2'
          betaLength = 22;
          betaGrid = zeros([nGrid, betaLength]);
  end

  if strcmp(windowType, 'spherical')
    % Determine reference distance
    refDist = distance(0, 180, 0+windowSizeMean, 180,...
                            referenceEllipsoid('WGS84', 'm'))
  end
  
  tic;
  parfor iGrid = 1:nGrid
  %randIdx = randi(nGrid);
  %randIdx = 9090;
      if ~mod(iGrid, floor(nGrid/20))
          disp([int2str(iGrid), '/', int2str(nGrid)]);
      end
      
      latSel = latGrid(iGrid);
      longSel = longGrid(iGrid);

      latMin = latSel - windowSizeMargined;
      latMax = latSel + windowSizeMargined;
      longMin = longSel - windowSizeMargined;
      longMax = longSel + windowSizeMargined;

      if isEqBorder
          if abs(latSel) <= eqBorder             
%              [iGridSub1, iGridSub2] = ind2sub(size(latGrid), iGrid);
              betaGrid(iGrid, :) = NaN;
              continue;
          end
          if latSel > 0
              latMin = max([latMin, eqBorder]);
          else
              latMax = min([latMax, - eqBorder]);
          end
      end

      idx = find(profLatAggrSel > latMin & profLatAggrSel < latMax & profLongAggrSel > longMin & profLongAggrSel < longMax);
      switch windowType
        case 'spherical'
          is_in_circle = distance(latSel, longSel, profLatAggrSel(idx), profLongAggrSel(idx),...
                            referenceEllipsoid('WGS84', 'm')) < refDist;
          idx = idx(is_in_circle);
        case 'box_var'
          idx = find(profLatAggrSel > latMin & profLatAggrSel < latMax & profLongAggrSel > longMin & profLongAggrSel < longMax);          
      end

      % Need at least 20 data points to estimate the regression coefficients, also do not compute mean if outside land/datamask
      if (length(idx) < 20) || isnan(mask(iGrid))
          disp(['lat: ',num2str(latSel),', lon: ',num2str(longSel),' unreliable grid point'])
%          [iGridSub1, iGridSub2] = ind2sub(size(latGrid), iGrid);
          betaGrid(iGrid, :) = NaN;
          continue;
      end
      
      profJulDayAggrWindow = profJulDayAggrSel(idx)';
      profYearAggrWindow = profYearAggrSel(idx)';
      profMonthAggrWindow = profMonthAggrSel(idx)';
      profLatAggrWindow = profLatAggrSel(idx)';
      profLongAggrWindow = profLongAggrSel(idx)';
      responseProfWindow = responseProf(idx)';
        
%{
      figure;
      histogram(responseProfWindow)
      print('-depsc2',['./Figures/',destFolder,'/','resid_',num2str(latSel),'_',num2str(longSel),meanTag,tag,verticalSelection,'_',num2str(windowSize),'.eps']);
      close all
      continue;
%}
      
      profYearDayAggrWindow = fromJulDayToYearDay(profJulDayAggrWindow);
      profYearLengthAggrWindow = yearLength(profJulDayAggrWindow);
      profYearDayRatioWindow = profYearDayAggrWindow ./ profYearLengthAggrWindow;
      
      % Setup Design Matrix
      switch meanTag
          case 'NoTrend'
              XWindow = [ones(length(profJulDayAggrWindow), 1) ...
                 sin(2*pi*1*profYearDayRatioWindow) cos(2*pi*1*profYearDayRatioWindow) ...
                 sin(2*pi*2*profYearDayRatioWindow) cos(2*pi*2*profYearDayRatioWindow) ...
                 sin(2*pi*3*profYearDayRatioWindow) cos(2*pi*3*profYearDayRatioWindow) ...
                 sin(2*pi*4*profYearDayRatioWindow) cos(2*pi*4*profYearDayRatioWindow) ...
                 sin(2*pi*5*profYearDayRatioWindow) cos(2*pi*5*profYearDayRatioWindow) ...
                 sin(2*pi*6*profYearDayRatioWindow) cos(2*pi*6*profYearDayRatioWindow) ...
                 (profLatAggrWindow-latSel) (profLongAggrWindow-longSel) (profLatAggrWindow-latSel).*(profLongAggrWindow-longSel) ...
                 (profLatAggrWindow-latSel).^2 (profLongAggrWindow-longSel).^2];
          case 'NoTrendVelSeas3'
              XWindow = [ones(length(profJulDayAggrWindow), 1) ...
                 sin(2*pi*1*profYearDayRatioWindow) cos(2*pi*1*profYearDayRatioWindow) ...
                 sin(2*pi*2*profYearDayRatioWindow) cos(2*pi*2*profYearDayRatioWindow) ...
                 sin(2*pi*3*profYearDayRatioWindow) cos(2*pi*3*profYearDayRatioWindow) ...
                 (profLatAggrWindow-latSel) (profLongAggrWindow-longSel) (profLatAggrWindow-latSel).*(profLongAggrWindow-longSel) ...
                 (profLatAggrWindow-latSel).^2 (profLongAggrWindow-longSel).^2 ...
                 (profLatAggrWindow-latSel).*sin(2*pi*1*profYearDayRatioWindow) ...
                 (profLatAggrWindow-latSel).*cos(2*pi*1*profYearDayRatioWindow) ...
                 (profLatAggrWindow-latSel).*sin(2*pi*2*profYearDayRatioWindow) ...
                 (profLatAggrWindow-latSel).*cos(2*pi*2*profYearDayRatioWindow) ...                 
                 (profLatAggrWindow-latSel).*sin(2*pi*3*profYearDayRatioWindow) ...
                 (profLatAggrWindow-latSel).*cos(2*pi*3*profYearDayRatioWindow) ...                 
                 (profLongAggrWindow-longSel).*sin(2*pi*1*profYearDayRatioWindow) ...
                 (profLongAggrWindow-longSel).*cos(2*pi*1*profYearDayRatioWindow) ...
                 (profLongAggrWindow-longSel).*sin(2*pi*2*profYearDayRatioWindow) ...
                 (profLongAggrWindow-longSel).*cos(2*pi*2*profYearDayRatioWindow) ...                 
                 (profLongAggrWindow-longSel).*sin(2*pi*3*profYearDayRatioWindow) ...
                 (profLongAggrWindow-longSel).*cos(2*pi*3*profYearDayRatioWindow) ...                 
                 ];
          case 'TrendDens'
              XWindow = [ones(length(profJulDayAggrWindow),1) ...
                 sin(2*pi*1*profYearDayRatioWindow) cos(2*pi*1*profYearDayRatioWindow) ...
                 sin(2*pi*2*profYearDayRatioWindow) cos(2*pi*2*profYearDayRatioWindow) ...
                 sin(2*pi*3*profYearDayRatioWindow) cos(2*pi*3*profYearDayRatioWindow) ...
                 sin(2*pi*4*profYearDayRatioWindow) cos(2*pi*4*profYearDayRatioWindow) ...
                 sin(2*pi*5*profYearDayRatioWindow) cos(2*pi*5*profYearDayRatioWindow) ...
                 sin(2*pi*6*profYearDayRatioWindow) cos(2*pi*6*profYearDayRatioWindow) ...
                 (profLatAggrWindow-latSel) (profLongAggrWindow-longSel) (profLatAggrWindow-latSel).*(profLongAggrWindow-longSel) ...
                 (profLatAggrWindow-latSel).^2 (profLongAggrWindow-longSel).^2 ...
                 (profJulDayAggrWindow-midJulDay)...
                 (profLatAggrWindow-latSel).*(profJulDayAggrWindow-midJulDay) (profLongAggrWindow-longSel).*(profJulDayAggrWindow-midJulDay)];
          case 'TrendDens2'
              XWindow = [ones(length(profJulDayAggrWindow),1) ...
                 sin(2*pi*1*profYearDayRatioWindow) cos(2*pi*1*profYearDayRatioWindow) ...
                 sin(2*pi*2*profYearDayRatioWindow) cos(2*pi*2*profYearDayRatioWindow) ...
                 sin(2*pi*3*profYearDayRatioWindow) cos(2*pi*3*profYearDayRatioWindow) ...
                 sin(2*pi*4*profYearDayRatioWindow) cos(2*pi*4*profYearDayRatioWindow) ...
                 sin(2*pi*5*profYearDayRatioWindow) cos(2*pi*5*profYearDayRatioWindow) ...
                 sin(2*pi*6*profYearDayRatioWindow) cos(2*pi*6*profYearDayRatioWindow) ...
                 (profLatAggrWindow-latSel) (profLongAggrWindow-longSel) (profLatAggrWindow-latSel).*(profLongAggrWindow-longSel) ...
                 (profLatAggrWindow-latSel).^2 (profLongAggrWindow-longSel).^2 ...
                 (profJulDayAggrWindow-midJulDay)...
                 (profLatAggrWindow-latSel).*(profJulDayAggrWindow-midJulDay) (profLongAggrWindow-longSel).*(profJulDayAggrWindow-midJulDay)...
                 (profJulDayAggrWindow-midJulDay).^2];
      end

      if (iterEM == 0) || (isFminError(iGrid))
        % Fit linear regression
        betaWindow = XWindow \ responseProfWindow;
      else
        % Construct covariance
        if isFullMonth
          thetasOpt_cur = cellfun(@(x) x(iGrid), thetasOpt);
          thetaLatOpt_cur = cellfun(@(x) x(iGrid), thetaLatOpt);
          thetaLongOpt_cur = cellfun(@(x) x(iGrid), thetaLongOpt);
          thetatOpt_cur = cellfun(@(x) x(iGrid), thetatOpt);
          sigmaOpt_cur = cellfun(@(x) x(iGrid), sigmaOpt);
        else
          thetasOpt_cur = thetasOpt(iGrid);
          thetaLatOpt_cur = thetaLatOpt(iGrid);
          thetaLongOpt_cur = thetaLongOpt(iGrid);
          thetatOpt_cur = thetatOpt(iGrid);
          sigmaOpt_cur = sigmaOpt(iGrid);
        end

        % accumulator
        quadMat = zeros(betaLength, betaLength);
        crossMat = zeros(betaLength, 1);
        for iYear = startYear:endYear
          isIyear = (profYearAggrWindow == iYear);

          if iYear == startYear
              isPrevYear = false(size(profYearAggrWindow));
          else
              isPrevYear = (profYearAggrWindow == (iYear-1));
          end

%{
          if iYear == endYear
            isFutureYear = false(size(profYearAggrWindow));
          else
            isFutureYear = (profYearAggrWindow == (iYear+1));
          end
%}

          % Window for 3 month within year
          % Firstyear, iMonth = 1:3
          if iYear == startYear
            monthWindow = [1, 2, 3];
            iMonthCenter = 2;
            isImonthWindow = isIyear & ismember(profMonthAggrWindow, monthWindow);
            isWindow = isImonthWindow;

            if isFullMonth
              covObs = feval(['spaceTimeCovariance',kernelType,'_vec'],...
                    profLatAggrWindow(isWindow), profLongAggrWindow(isWindow), profJulDayAggrWindow(isWindow),...
                    thetasOpt_cur(iMonthCenter),thetaLatOpt_cur(iMonthCenter),thetaLongOpt_cur(iMonthCenter),thetatOpt_cur(iMonthCenter));
              covObs(1:(sum(isWindow)+1):end) = thetasOpt_cur(iMonthCenter) + sigmaOpt_cur(iMonthCenter) .^ 2;
            else
              covObs = feval(['spaceTimeCovariance',kernelType,'_vec'],...
                  profLatAggrWindow(isWindow), profLongAggrWindow(isWindow), profJulDayAggrWindow(isWindow),...
                  thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur);
              covObs(1:(sum(isWindow)+1):end) = thetasOpt_cur + sigmaOpt_cur .^ 2;
            end

            KinvEta = covObs \ XWindow(isWindow,:);  % isWindow * betaLength
            quadMat = quadMat + (XWindow(isWindow,:)' * KinvEta);
            crossMat = crossMat + (KinvEta' * responseProfWindow(isWindow));
          end

          for iMonth = 1:12
            if (iYear == startYear) && (iMonth < 4)
              continue;
            end

            if iMonth == 1
              monthWindow = 1;
              extendMonthWindow = [11, 12];
              iMonthCenter = 12;
            elseif iMonth == 2
              monthWindow = 1:2;
              extendMonthWindow = 12;
              iMonthCenter = iMonth - 1;
            else
              monthWindow = iMonth + (-2:0);%mod(iMonth + (-2:0), 12);
              extendMonthWindow = [];
              iMonthCenter = iMonth - 1;
            end
            
            isImonthWindow = isIyear & ismember(profMonthAggrWindow, monthWindow);
            isWindow = isImonthWindow | (isPrevYear & ismember(profMonthAggrWindow, extendMonthWindow));
            isWindowImonth = isWindow & (profMonthAggrWindow == iMonth);
            
            isImonth = (profMonthAggrWindow(isWindow) == iMonth); % iMonth within Window

            
            if isFullMonth
              covObs = feval(['spaceTimeCovariance',kernelType,'_vec'],...
                    profLatAggrWindow(isWindow), profLongAggrWindow(isWindow), profJulDayAggrWindow(isWindow),...
                    thetasOpt_cur(iMonthCenter),thetaLatOpt_cur(iMonthCenter),thetaLongOpt_cur(iMonthCenter),thetatOpt_cur(iMonthCenter));
              covObs(1:(sum(isWindow)+1):end) = thetasOpt_cur(iMonthCenter) + sigmaOpt_cur(iMonthCenter) .^ 2;
            else
              covObs = feval(['spaceTimeCovariance',kernelType,'_vec'],...
                  profLatAggrWindow(isWindow), profLongAggrWindow(isWindow), profJulDayAggrWindow(isWindow),...
                  thetasOpt_cur,thetaLatOpt_cur,thetaLongOpt_cur,thetatOpt_cur);
              covObs(1:(sum(isWindow)+1):end) = thetasOpt_cur + sigmaOpt_cur .^ 2;
            end

            A = covObs(isImonth, ~isImonth) / covObs(~isImonth, ~isImonth); % isImonth * ~isImonth
            etaDelta = XWindow(isWindowImonth,:) - A * XWindow(~isWindowImonth & isWindow, :);
            KinvEta = (covObs(isImonth, isImonth) - A * covObs(~isImonth, isImonth)) \ etaDelta;
            quadMat = quadMat + (etaDelta' * KinvEta);
            crossMat = crossMat + (KinvEta' * (responseProfWindow(isWindowImonth) - A * responseProfWindow(~isWindowImonth & isWindow)));
      %      disp(min(eig(covObs(isWindow, isWindow))))
          end
        end
        
        betaWindow = quadMat \ crossMat;
      end

      betaGrid(iGrid, :) = betaWindow;
  end
  toc;

  betaGrid = reshape(betaGrid, [size(latGrid), betaLength]);

  if is2step
      if isnumeric(verticalSelection)
        saveName = ['./Results/meanField',typeTag,fluxType,responseTag,meanTag,tag,presString,dataYear,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_',num2str(minNumberOfObs),'_Eq',num2str(eqBorder),EMTag,'.mat']
      else
        saveName = ['./Results/meanField',typeTag,fluxType,responseTag,meanTag,tag,verticalSelection,dataYear,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_',num2str(minNumberOfObs),'_Eq',num2str(eqBorder),EMTag,'.mat']
      end
      save(saveName,...
            'betaGrid','latGrid','longGrid','midJulDay', 'fluxType');
  else
      saveName = ['./Results/meanField',responseTag,meanTag,tag,verticalSelection,dataYear,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_',num2str(minNumberOfObs),EMTag,'.mat']
      save(saveName,...
          'betaGrid','latGrid','longGrid','midJulDay');
  end

end