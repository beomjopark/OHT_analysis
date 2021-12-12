function subtractMeanSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isPlot, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth)
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

  if isFullMonth && ~isempty(EMTag)
    EMTag = ['Full', EMTag]
  end
  if isFullMonth && ~isempty(prevEMTag)
    prevEMTag = ['Full', prevEMTag]
  end

  YFtag = []%'_YF3NN'

  if is2step
    dataName = ['./Data/',typeTag,fluxType,responseTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'.mat']
    srcName = ['./Results/meanField',typeTag,fluxType,responseTag,meanTag,tag,verticalSelection,dataYear,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_',num2str(minNumberOfObs),'_Eq',num2str(eqBorder),EMTag,'.mat']
  else
    switch responseTag
        case {'Temp', 'Dens'}
          dataName = ['./Data/',typeTag,'TempDens','Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']
          srcName = ['./Results/','meanField',responseTag,meanTag,tag,verticalSelection,dataYear,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_',num2str(minNumberOfObs),YFtag,EMTag,'.mat']
%            load(['./Data/','pre_thresh/',typeTag,'TempDens','Prof',tag,verticalSelection,dataYear,adjustNumTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat']);
%            load(['./Results/','prethresh/','meanField',responseTag,meanTag,tag,verticalSelection,dataYear,adjustNumTag,absoluteTag,windowTypeTag,'_w',num2str(windowSize),'_',num2str(minNumberOfObs),YFtag,'.mat']);

            if exist('is_thresh_filter', 'var')
                disp('The Profiles are thresholded.')
            end
        case {'DUACS', 'ESA', 'SOSITemp'}
            dataName = ['./Data/',typeTag,responseTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']
            srcName = ['./Results/','meanField',responseTag,meanTag,tag,verticalSelection,dataYear,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_',num2str(minNumberOfObs),EMTag,'.mat']
        case 'Sal'
            dataName = ['./Data/',typeTag,'Prof',tag,verticalSelection,'Filtered_',num2str(minNumberOfObs),'.mat']
            srcName = ['./Results/meanField',responseTag,meanTag,tag,verticalSelection,'_',windowSizeTag,'_',num2str(minNumberOfObs),'.mat']
    end
  end
  load(dataName);
  load(srcName);
  % Forgot to change the sign
%{
  if strcmp(fluxType, 'vol') && strcmp(typeTag, 'lat')
      profLatFluxAggrSel = - profLatFluxAggrSel;
      betaGrid = - betaGrid;
  end
%}

  %% Subtract mean
  profLatAggrSelRounded = roundHalf(profLatAggrSel);
  profLongAggrSelRounded = roundHalf(profLongAggrSel);
  nProf = length(profLatAggrSelRounded);

%  profRes = zeros(1, nProf);
  if is2step
    switch typeTag
      case {'intlat', 'intlon'}
        profResponse = intFluxProf;
      otherwise
        profResponse = profFluxAggrSel;
    end
  else
    switch responseTag
        case 'Temp'
          profResponse = intTempProf;
        case 'Dens'
          profResponse = targetDynhProf;
        case 'DUACS'
          profResponse = targetADTProf;
        case {'ESA', 'SOSITemp'}
          profResponse = targetSSTProf;          
        case 'Sal'
          profResponse = targetSalProfPchip;
    end
  end


  profYearDayRatio = fromJulDayToYearDay(profJulDayAggrSel) ./ yearLength(profJulDayAggrSel);
  ResponseProfHat = zeros(1, nProf);
  for iProf = 1:nProf
      if ~mod(iProf, floor(nProf/20))
          disp([int2str(iProf), '/', int2str(nProf)]);
      end

      profLat = profLatAggrSelRounded(iProf);
      profLong = profLongAggrSelRounded(iProf);

      iLat = find(abs(latGrid(1,:) - profLat) < 1e-03);
      iLong = find(abs(longGrid(:,1) - profLong) < 1e-03);

      betaTemp = squeeze(betaGrid(iLong,iLat,:));
      switch meanTag
          case 'NoTrend'
%{
              fun = @(x) betaTemp(1) ...
                 + betaTemp(2) * sin(2*pi*1*profYearDayRatio(iProf))' + betaTemp(3) * cos(2*pi*1*profYearDayRatio(iProf))' ...
                 + betaTemp(4) * sin(2*pi*2*profYearDayRatio(iProf))' + betaTemp(5) * cos(2*pi*2*profYearDayRatio(iProf))' ...
                 + betaTemp(6) * sin(2*pi*3*profYearDayRatio(iProf))' + betaTemp(7) * cos(2*pi*3*profYearDayRatio(iProf))' ...
                 + betaTemp(8) * sin(2*pi*4*profYearDayRatio(iProf))' + betaTemp(9) * cos(2*pi*4*profYearDayRatio(iProf))' ...
                 + betaTemp(10) * sin(2*pi*5*profYearDayRatio(iProf))' + betaTemp(11) * cos(2*pi*5*profYearDayRatio(iProf))' ...
                 + betaTemp(12) * sin(2*pi*6*profYearDayRatio(iProf))' + betaTemp(13) * cos(2*pi*6*profYearDayRatio(iProf))';% ...                 %+ betaTemp(19) * (x - midJulDay) + betaTemp(20) * (x - midJulDay).^2;
%}

            Xdesign = [1 sin(2*pi*1*profYearDayRatio(iProf)) cos(2*pi*1*profYearDayRatio(iProf)) ...
             sin(2*pi*2*profYearDayRatio(iProf)) cos(2*pi*2*profYearDayRatio(iProf)) ...
             sin(2*pi*3*profYearDayRatio(iProf)) cos(2*pi*3*profYearDayRatio(iProf)) ...
             sin(2*pi*4*profYearDayRatio(iProf)) cos(2*pi*4*profYearDayRatio(iProf)) ...
             sin(2*pi*5*profYearDayRatio(iProf)) cos(2*pi*5*profYearDayRatio(iProf)) ...
             sin(2*pi*6*profYearDayRatio(iProf)) cos(2*pi*6*profYearDayRatio(iProf)) ...
             (profLatAggrSel(iProf) - profLat) (profLongAggrSel(iProf) - profLong) ...
             (profLatAggrSel(iProf) - profLat).*(profLongAggrSel(iProf) - profLong)...
             (profLatAggrSel(iProf) - profLat).^2 (profLongAggrSel(iProf) - profLong).^2];
             fun = @(x) Xdesign * betaTemp;
          case 'NoTrendVelSeas3'
%{
              fun = @(x) betaTemp(1) ...
                 + betaTemp(2) * sin(2*pi*1*profYearDayRatio(iProf))' + betaTemp(3) * cos(2*pi*1*profYearDayRatio(iProf))' ...
                 + betaTemp(4) * sin(2*pi*2*profYearDayRatio(iProf))' + betaTemp(5) * cos(2*pi*2*profYearDayRatio(iProf))' ...
                 + betaTemp(6) * sin(2*pi*3*profYearDayRatio(iProf))' + betaTemp(7) * cos(2*pi*3*profYearDayRatio(iProf))';% ...
                 %+ betaTemp(19) * (x - midJulDay) + betaTemp(20) * (x - midJulDay).^2;
%}

            Xdesign = [1 sin(2*pi*1*profYearDayRatio(iProf)) cos(2*pi*1*profYearDayRatio(iProf)) ...
             sin(2*pi*2*profYearDayRatio(iProf)) cos(2*pi*2*profYearDayRatio(iProf)) ...
             sin(2*pi*3*profYearDayRatio(iProf)) cos(2*pi*3*profYearDayRatio(iProf)) ...
             (profLatAggrSel(iProf) - profLat) (profLongAggrSel(iProf) - profLong) ...
             (profLatAggrSel(iProf) - profLat).*(profLongAggrSel(iProf) - profLong)...
             (profLatAggrSel(iProf) - profLat).^2 (profLongAggrSel(iProf) - profLong).^2 ...
             (profLatAggrSel(iProf) - profLat).*sin(2*pi*1*profYearDayRatio(iProf))...
             (profLatAggrSel(iProf) - profLat).*cos(2*pi*1*profYearDayRatio(iProf)) ...
             (profLatAggrSel(iProf) - profLat).*sin(2*pi*2*profYearDayRatio(iProf))...
             (profLatAggrSel(iProf) - profLat).*cos(2*pi*2*profYearDayRatio(iProf)) ...
             (profLatAggrSel(iProf) - profLat).*sin(2*pi*3*profYearDayRatio(iProf))...
             (profLatAggrSel(iProf) - profLat).*cos(2*pi*3*profYearDayRatio(iProf)) ...
             (profLongAggrSel(iProf) - profLong).*sin(2*pi*1*profYearDayRatio(iProf))...
             (profLongAggrSel(iProf) - profLong).*cos(2*pi*1*profYearDayRatio(iProf)) ...
             (profLongAggrSel(iProf) - profLong).*sin(2*pi*2*profYearDayRatio(iProf))...
             (profLongAggrSel(iProf) - profLong).*cos(2*pi*2*profYearDayRatio(iProf)) ...
             (profLongAggrSel(iProf) - profLong).*sin(2*pi*3*profYearDayRatio(iProf))...
             (profLongAggrSel(iProf) - profLong).*cos(2*pi*3*profYearDayRatio(iProf))...
             ];
             fun = @(x) Xdesign * betaTemp;

          case 'Trend'
              fun = @(x) betaTemp(1) ...
                 + betaTemp(2) * sin(2*pi*1*profYearDayRatio(iProf))' + betaTemp(3) * cos(2*pi*1*profYearDayRatio(iProf))' ...
                 + betaTemp(4) * sin(2*pi*2*profYearDayRatio(iProf))' + betaTemp(5) * cos(2*pi*2*profYearDayRatio(iProf))' ...
                 + betaTemp(6) * sin(2*pi*3*profYearDayRatio(iProf))' + betaTemp(7) * cos(2*pi*3*profYearDayRatio(iProf))' ...
                 + betaTemp(8) * sin(2*pi*4*profYearDayRatio(iProf))' + betaTemp(9) * cos(2*pi*4*profYearDayRatio(iProf))' ...
                 + betaTemp(10) * sin(2*pi*5*profYearDayRatio(iProf))' + betaTemp(11) * cos(2*pi*5*profYearDayRatio(iProf))' ...
                 + betaTemp(12) * sin(2*pi*6*profYearDayRatio(iProf))' + betaTemp(13) * cos(2*pi*6*profYearDayRatio(iProf))' ...
                 + betaTemp(19) * (x - midJulDay);% + betaTemp(20) * (x - midJulDay).^2;
          case 'Trend2'
              fun = @(x) betaTemp(1) ...
                 + betaTemp(2) * sin(2*pi*1*profYearDayRatio(iProf))' + betaTemp(3) * cos(2*pi*1*profYearDayRatio(iProf))' ...
                 + betaTemp(4) * sin(2*pi*2*profYearDayRatio(iProf))' + betaTemp(5) * cos(2*pi*2*profYearDayRatio(iProf))' ...
                 + betaTemp(6) * sin(2*pi*3*profYearDayRatio(iProf))' + betaTemp(7) * cos(2*pi*3*profYearDayRatio(iProf))' ...
                 + betaTemp(8) * sin(2*pi*4*profYearDayRatio(iProf))' + betaTemp(9) * cos(2*pi*4*profYearDayRatio(iProf))' ...
                 + betaTemp(10) * sin(2*pi*5*profYearDayRatio(iProf))' + betaTemp(11) * cos(2*pi*5*profYearDayRatio(iProf))' ...
                 + betaTemp(12) * sin(2*pi*6*profYearDayRatio(iProf))' + betaTemp(13) * cos(2*pi*6*profYearDayRatio(iProf))' ...
                 + betaTemp(19) * (x - midJulDay) + betaTemp(20) * (x - midJulDay).^2;
      end

      ResponseProfHat(iProf) = fun(profJulDayAggrSel(iProf));
  end
  profRes = profResponse - ResponseProfHat;

  % filter NaN from datamask
  isNanRes = isnan(ResponseProfHat);
  disp(['Data Mask Filtered ', num2str(sum(~isNanRes)), ' out of ', num2str(numel(ResponseProfHat))]);

  profLatAggrSelNan = profLatAggrSel(isNanRes);
  profLongAggrSelNan = profLongAggrSel(isNanRes);

  profLatAggrSel = profLatAggrSel(~isNanRes);
  profLongAggrSel = profLongAggrSel(~isNanRes);
  profJulDayAggrSel = profJulDayAggrSel(~isNanRes);

  if is2step
    switch typeTag
      case {'intlat', 'intlon'}
          intFluxRes = profRes(~isNanRes);
          if isStandardize
            stdRes = std(intFluxRes, 'omitnan');
            intFluxRes = intFluxRes ./ stdRes;
          end
          save(['./Data/',typeTag, fluxType,responseTag,'Res',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),EMTag,'.mat'],...
      'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','intFluxRes', 'stdRes', 'isNanRes');
      otherwise
          FluxRes = profRes(~isNanRes);
          if isStandardize
            stdRes = std(FluxRes, 'omitnan');
            FluxRes = FluxRes ./ stdRes;
          end
          if isempty(targetTempProf)
            save(['./Data/',typeTag, fluxType,responseTag,'Res',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),EMTag,'.mat'],...
              'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','FluxRes', 'stdRes', 'isNanRes');
          else
            targetTempProf = targetTempProf(~isNanRes);
            save(['./Data/',typeTag, fluxType,responseTag,'Res',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),EMTag,'.mat'],...
              'profLatAggrSel','profLongAggrSel','profJulDayAggrSel', 'targetTempProf', 'FluxRes', 'stdRes', 'isNanRes');
          end
    end
  else
    switch responseTag
        case 'Temp'
            targetTempRes = profRes(~isNanRes);          
            if isStandardize
              stdRes = std(targetTempRes, 'omitnan');
              targetTempRes = targetTempRes ./ stdRes;
            end
            save(['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,EMTag,'.mat'],...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','targetTempRes', 'stdRes', 'isNanRes');
        case 'Dens'
            intDensRes = profRes(~isNanRes);          
            if isStandardize
              stdRes = std(intDensRes, 'omitnan');
              intDensRes = intDensRes ./ stdRes;
            end
            save(['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,EMTag,'.mat'],...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','intDensRes', 'stdRes', 'isNanRes');
        case 'DUACS'
            targetADTRes = profRes(~isNanRes);          
            if isStandardize
              stdRes = std(targetADTRes, 'omitnan');
              targetADTRes = targetADTRes ./ stdRes;
            end
            save(['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,EMTag,'.mat'],...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','targetADTRes', 'stdRes', 'isNanRes');
        case {'ESA', 'SOSITemp'}
            targetSSTRes = profRes(~isNanRes);          
            if isStandardize
              stdRes = std(targetSSTRes, 'omitnan');
              targetSSTRes = targetSSTRes ./ stdRes;
            end
            save(['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,EMTag,'.mat'],...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','targetSSTRes', 'stdRes', 'isNanRes');
        case 'Sal'
            targetSalRes = profRes(~isNanRes);          
            if isStandardize
              stdRes = std(targetSalRes, 'omitnan');
              targetSalRes = targetSalRes ./ stdRes;
            end
            save(['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),EMTag,'.mat'],...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','targetSalRes', 'stdRes', 'isNanRes');          
    end
  end

  if isPlot
    % Create Folder
    if is2step
      destFolder = [typeTag,fluxType,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_Eq',num2str(eqBorder)];
    else
      destFolder = [typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,meanTag,YFtag,'/mean',meanTag,windowTypeTag,windowSizeTag];
    end
    mkdir(['./Figures/',destFolder])

    %% Check filtered parts
    %rng(12345);
    plotIdx = 1:nProf; %randperm(nProf);

    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');

    load coast;
    plotm(lat,long,'k');

    scatterm(profLatAggrSelNan,profLongAggrSelNan,10,'r.');

    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))    
    if is2step
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,fluxType,responseTag,verticalSelection,dataYear,windowTypeTag,'NaNProf.eps']);
    else
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,responseTag,verticalSelection,dataYear,windowTypeTag,'NaNProf.eps']);
    end    
    
    %% Plot integrated conservative temperatures
    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');

    load coast;
    plotm(lat,long,'k');
    
    profResponseClean = profResponse(~isNanRes);
    if is2step
        % Border Mask
        mask = ~(profLatAggrSel < eqBorder & profLatAggrSel > - eqBorder) .* 1;
        mask(mask == 0) = NaN;
        
        scatterm(profLatAggrSel,profLongAggrSel,10,mask.*profResponseClean,'.');
        caxis([quantile(profResponseClean(:), 0.01), quantile(profResponseClean(:), 0.99)]);
        cLims = caxis;
        colormap(darkb2r(cLims(1),cLims(2)));
    else
        scatterm(profLatAggrSel,profLongAggrSel,10,profResponseClean,'.');
        %scatterm(profLatAggrSel,profLongAggrSel,10,targetTempProfPchip,'.');    
        caxis([quantile(profResponseClean(:), 0.01), quantile(profResponseClean(:), 0.99)]);
        colormap(parula(100));
    end

    colorbar;

    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    if is2step
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,fluxType,responseTag,verticalSelection,dataYear,windowTypeTag,'Prof.eps']);
    else
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,responseTag,verticalSelection,dataYear,windowTypeTag,'Prof.eps']);
    end
      


    %% Plot fitted integrated conservative temperatures
    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');

    load coast;
    plotm(lat,long,'k');
    profResClean = ResponseProfHat(~isNanRes);

    if is2step
        % Border Mask
        mask = ~(profLatAggrSel < eqBorder & profLatAggrSel > - eqBorder) .* 1;
        mask(mask == 0) = NaN;
        
        scatterm(profLatAggrSel,profLongAggrSel,10,mask.*profResClean,'.');
        caxis([quantile(profResClean(:), 0.01), quantile(profResClean(:), 0.99)]);
        cLims = caxis;
        colormap(darkb2r(cLims(1),cLims(2)));
    else
        scatterm(profLatAggrSel,profLongAggrSel,10,profResClean,'.');
        %scatterm(profLatAggrSel,profLongAggrSel,10,targetTempProfPchip,'.');    
        caxis([quantile(profResClean(:), 0.01), quantile(profResClean(:), 0.99)]);
        colormap(parula(100));
    end
    %colormap(darkb2r(cLims(1),cLims(2)));
    %caxis([1.844e06,1.848e06]);
    colorbar;

    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    if is2step
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,fluxType,responseTag,verticalSelection,dataYear,'ProfHat_',windowTypeTag,num2str(windowSize),'.eps']);
    else
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,responseTag,verticalSelection,dataYear,'ProfHat_',windowTypeTag,num2str(windowSize),'.eps']);
    end


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

%    colormap(parula(100));
    colorbar;

    cutoff = max(abs([quantile(temp(:),0.01), quantile(temp(:),0.99)]));
    caxis([-cutoff, cutoff]);
    cLims = caxis;
    colormap(darkb2r(cLims(1), cLims(2)));

    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    if is2step
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,fluxType,responseTag,verticalSelection,dataYear,'Res_',windowTypeTag,num2str(windowSize),'.eps']);
    else
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,responseTag,verticalSelection,dataYear,'Res_',windowTypeTag,num2str(windowSize),'.eps']);
    end

    %% Plot Transport : Deprecated. See plotMeanField
    %{
    if is2step
      figure;
      handle = worldmap('World');
      setm(handle, 'Origin', [0 200 0]);
      tightmap;
      mlabel('off');
      plabel('off');

      load coast;
      plotm(lat,long,'k');

      isIntFlux = (numel(typeTag) > 3);
      if isIntFlux %intlat / int lon
        typeTag = typeTag(4:end);
        % Actually integrated need gravity pre multiplied...
        intMid = mean([intStart, intEnd]);
        gravityProf =  gsw_grav(profLatAggrSel, intMid);
      else
        gravityProf = gsw_grav(profLatAggrSel, intStart);
      end

      switch fluxType
        case 'heat'
          profHtGrid = gsw_cp0 .* ResponseProfHat ./ gravityProf;
          profHtGrid = mask.*profHtGrid.* 1e-15;
          surfm(latGrid,longGrid,profHtGrid);
          caxis([quantile(profHtGrid(:), 0.003), quantile(profHtGrid(:), 0.997)]);

          cb = colorbar;
          cb.Label.String = 'Heat Transport [PW]';  
        case 'mass'
          profMtGrid = ResponseProfHat ./ gravityProf;
          profMtGrid = mask.*profMtGrid;
          surfm(latGrid,longGrid,profMtGrid);
          caxis([quantile(profMtGrid(:), 0.003), quantile(profMtGrid(:), 0.997)]);

          cb = colorbar;
          cb.Label.String = 'Mass Transport [kg / s]';  
        case 'vol'
          profVtGrid = ResponseProfHat ./ gravityProf;
          profVtGrid = mask.*profVtGrid.* 1e-6;
          surfm(latGrid,longGrid,profVtGrid);
          caxis([quantile(profVtGrid(:), 0.003), quantile(profVtGrid(:), 0.997)]);

          cb = colorbar;
          cb.Label.String = 'Volume Transport [Sv]';          
      end
      cb.FontSize = 14;

      cLims = caxis;
      colormap(darkb2r(cLims(1),cLims(2)));

      set(gcf,'units','centimeters')
      set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
      set(gcf,'paperunits',get(gcf,'units')) 
      set(gcf,'paperpos',get(gcf,'pos'))
      switch typeTag
          case 'lat'
              print('-depsc2',['./Figures/',destFolder,'/','ZhtProf',verticalSelection,dataYear,'Hat.eps']); %lat
          case 'lon'
              print('-depsc2',['./Figures/',destFolder,'/','MhtProf',verticalSelection,dataYear,'Hat.eps']); %lon
      end
    end
    %}
  end

end
