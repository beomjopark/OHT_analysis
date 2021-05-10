function adjustResidual(month, kernelType, meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isPlot, isStandardize, fluxType, eqBorder, nAdjust, iterEM, isFullMonth)
  %% Subtract MeanField from the observation
  if nargin < 10
    isStandardize = false;
  end
  if ~isStandardize
      stdRes = 1;
  else
      fprintf("Standardized!")
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

  adjustNumTag = ['Adjusted', num2str(nAdjust+1)]; % Increase adjustment counter
  if nAdjust == 0
    adjustPrevNumTag = [];
  else
    adjustPrevNumTag = ['Adjusted', num2str(nAdjust)];
  end

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

  EMOutTag = EMTag;
  
  tag = 'PchipPotTemp';
  YFtag = [];%'_YF3NN'
  if is2step
    % load(['./Data/',typeTag,fluxType,responseTag,'Prof',tag,verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'_w',num2str(windowSize),'.mat']);
    load(['./Results/meanField',typeTag,fluxType,responseTag,meanTag,tag,verticalSelection,dataYear,adjustTag,[],windowTypeTag,'_w',windowSizeTag,'_',num2str(minNumberOfObs),'_Eq',num2str(eqBorder),EMTag,'.mat'], 'latGrid', 'longGrid');
  else
    load(['./Results/meanField',responseTag,meanTag,tag,verticalSelection,dataYear,adjustPrevNumTag,[],windowTypeTag,'_w',windowSizeTag,'_',num2str(minNumberOfObs),EMTag,'.mat'], 'latGrid', 'longGrid');
  end
  
  if is2step
    load(['./Data/',typeTag, fluxType,responseTag,'Res',verticalSelection,dataYear,adjustPrevNumTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'_Eq',num2str(eqBorder),'.mat']);

    % Recover back into original scale
    switch typeTag
      case 'lat'
          Res = latFluxRes .* stdRes;
      case 'lon'
          Res = lonFluxRes .* stdRes;
      otherwise
          Res = intFluxRes .* stdRes;
    end

    % Load meanPredGrid
    srcFolder = ['anomaly_',typeTag,fluxType,responseTag,adjustPrevNumTag,windowTypeTag,'_w',num2str(windowSize),'_Eq',num2str(eqBorder),'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d')];
    load(['./Results/',srcFolder,'/MeanAnomaly',typeTag,fluxType,responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'.mat']);
  else
    load(['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustPrevNumTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,EMTag,'.mat']);
    switch responseTag
        case 'Temp'
            Res = targetTempRes .* stdRes;
        case 'ESA'
            Res = targetSSTRes .* stdRes;
        case 'Dens'
            Res = intDensRes .* stdRes;
        case 'DUACS'
            Res = targetADTRes .* stdRes;
        case 'Sal'
            Res = targetSalRes .* stdRes;
    end

    % Load meanPredGrid
    srcFolder = ['anomaly_',typeTag,responseTag,adjustPrevNumTag,[],windowTypeTag,'_w',windowSizeFullTag,'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag]
    load(['./Results/',srcFolder,'/MeanAnomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'.mat']);
  end
 
  %% Subtract mean from the closest grid
  profLatAggrSelRounded = roundHalf(profLatAggrSel);
  profLongAggrSelRounded = roundHalf(profLongAggrSel);
  nProf = length(profLatAggrSelRounded);

%  meanProf = zeros(size(intDensRes));
  for iProf = 1:nProf
      if ~mod(iProf, floor(nProf/20))
          disp([int2str(iProf), '/', int2str(nProf)]);
      end

      profLat = profLatAggrSelRounded(iProf);
      profLong = profLongAggrSelRounded(iProf);

      iLat = find(abs(latGrid(1,:) - profLat) < 1e-03);
      iLong = find(abs(longGrid(:,1) - profLong) < 1e-03);
      Res(iProf) = Res(iProf) - meanPredGrid(iLong, iLat);
  end

  % Save Data
  if isStandardize
    stdRes = std(Res, 'omitnan');
  end

  if is2step
    saveName = ['./Data/',typeTag, fluxType,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeFullTag,EMTag,'_Eq',num2str(eqBorder),'.mat']
    switch typeTag
      case 'lat'
          if isStandardize
            latFluxRes = Res ./ stdRes;
          end
          save(saveName,...
      'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','latFluxRes', 'stdRes');
      case 'lon'
          if isStandardize
            lonFluxRes = Res ./ stdRes;
          end
          save(saveName,...
      'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','lonFluxRes', 'stdRes');
      otherwise
          if isStandardize
            intFluxRes = Res ./ stdRes;
          end
          save(saveName,...
      'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','intFluxRes', 'stdRes');
    end
  else
    saveName = ['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeFullTag,EMTag,'.mat']
    switch responseTag
        case 'Temp'
            if isStandardize
              targetTempRes = Res ./ stdRes;
            end
            save(saveName,...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','targetTempRes', 'stdRes');
        case 'Dens'
            if isStandardize
              intDensRes = Res ./ stdRes;
            end
            save(saveName,...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','intDensRes', 'stdRes');
        case 'ESA'
            if isStandardize
              targetSSTRes = Res ./ stdRes;
            end
            save(saveName,...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','targetSSTRes', 'stdRes');            
        case 'DUACS'
            if isStandardize
              targetADTRes = Res ./ stdRes;
            end
            save(saveName,...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','targetADTRes', 'stdRes');
        case 'Sal'
            if isStandardize
              targetSalRes = Res ./ stdRes;
            end
            save(saveName,...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','targetSalRes', 'stdRes');
    end
  end
  if isPlot
    % Create Folder
    if is2step
        destFolder = [typeTag,fluxType,responseTag,verticalSelection,dataYear,adjustNumTag,[],windowTypeTag,'_Eq',num2str(eqBorder),...
                 '/','Post_','Anomaly_',kernelType,'_',windowSizeTag,EMTag]
    else
        profileTag = []; %''_Profile_YF3''
        destFolder = [typeTag,responseTag,verticalSelection,dataYear,adjustNumTag,[],windowTypeTag, meanTag,profileTag,...
                 '/','Post_','Anomaly_',kernelType,'_',windowSizeTag,EMTag,'_month',num2str(month)]
    end
    mkdir(['./Figures/',destFolder])


    %% Plot integrated conservative temperature residuals
    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');

    load coast;
    plotm(lat,long,'k');

    if is2step
      % Border Mask
      mask = ~(profLatAggrSel < eqBorder & profLatAggrSel > - eqBorder) .* 1;
      mask(mask == 0) = NaN;
      switch typeTag
        case 'lat'
            scatterm(profLatAggrSel,profLongAggrSel,10,mask.*latFluxRes,'.');
            temp = mask.*latFluxRes;
        case 'lon'
            scatterm(profLatAggrSel,profLongAggrSel,10,mask.*lonFluxRes,'.');
            temp = mask.*lonFluxRes;
        otherwise
            scatterm(profLatAggrSel,profLongAggrSel,10,mask.*intFluxRes,'.');
            temp = mask.*intFluxRes;
      end
    else
      switch responseTag
          case 'Temp'
              scatterm(profLatAggrSel,profLongAggrSel,10,targetTempRes,'.');
              temp = targetTempRes;
          case 'Dens'
              scatterm(profLatAggrSel,profLongAggrSel,10,intDensRes,'.');
              temp = intDensRes;
      end
    end

    cutoff = max(abs([quantile(temp(:),0.01), quantile(temp(:),0.99)]));
    caxis([-cutoff, cutoff]);
    cLims = caxis;
    colormap(darkb2r(cLims(1), cLims(2)));

    colorbar;
    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    if is2step
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,fluxType,responseTag,verticalSelection,dataYear,'Res_',num2str(windowSize),'.eps']);
    else
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,responseTag,verticalSelection,dataYear,'Res_',num2str(windowSize),'.eps']);
    end

  end

end
