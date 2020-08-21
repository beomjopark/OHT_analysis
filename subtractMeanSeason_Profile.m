function subtractMeanSeason_Profile(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowSize, minNumberOfObs, is2step, isPlot, isStandardize, fluxType, eqBorder, h_YFParam)
  %% Subtract MeanField from the observation
  if nargin < 10
    isStandardize = false;
  end
  if ~isStandardize
      stdRes = 1;
  else
      fprintf("Standardized!")
  end

  if ~isempty(h_YFParam)
    profileTag = ['_Profile_YF',num2str(h_YFParam)];
  else
    profileTag = '_Profile';
  end

  tag = 'PchipPotTemp';

  if is2step
    load(['./Data/',typeTag,fluxType,responseTag,'Prof',tag,verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'_w',num2str(windowSize),'.mat']);
    load(['./Results/meanField',typeTag,fluxType,responseTag,meanTag,tag,verticalSelection,dataYear,'_w',num2str(windowSize),'_',num2str(minNumberOfObs),'_Eq',num2str(eqBorder),profileTag,'.mat']);
  else
    switch responseTag
        case 'Temp'
            switch verticalSelection
                case 'MidMeso'
                    load(['./Data/',typeTag,responseTag,'Prof600','Filtered_',num2str(minNumberOfObs),'.mat']);
                    load(['./Results/meanField',responseTag,meanTag,'target600','_',num2str(windowSize),'_',num2str(minNumberOfObs),profileTag,'.mat']);
                case 'FullMeso'
                    load(['./Data/',typeTag,responseTag,'Prof200','Filtered_',num2str(minNumberOfObs),'.mat']);
                    load(['./Results/meanField',responseTag,meanTag,'target200','_',num2str(windowSize),'_',num2str(minNumberOfObs),profileTag,'.mat']);
            end
        case 'Dens'
            load(['./Data/',typeTag,'TempDens','Prof',tag,verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'_w',num2str(windowSize),'.mat']);
            load(['./Results/meanField',responseTag,meanTag,tag,verticalSelection,dataYear,'_w',num2str(windowSize),'_',num2str(minNumberOfObs),profileTag,'.mat']);
        case 'Sal'
            load(['./Data/',typeTag,'Prof',tag,verticalSelection,'Filtered_',num2str(minNumberOfObs),'.mat']);
            load(['./Results/meanField',responseTag,meanTag,tag,verticalSelection,'_',num2str(windowSize),'_',num2str(minNumberOfObs),profileTag,'.mat']);
    end
  end
  
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

  if is2step
    switch typeTag
      case 'lat'
        latFluxRes = zeros(1,nProf);
      case 'lon'
        lonFluxRes = zeros(1,nProf); 
      otherwise
        intFluxRes = zeros(1,nProf);
    end
  else
    switch responseTag
        case 'Temp'
            nProf = length(targetTempProfPchip);
            targetTempRes = zeros(1,nProf);
        case 'Dens'
            nProf = length(targetDynhProf);
            intDensRes = zeros(1,nProf);
        case 'Sal'
            nProf = length(targetSalProfPchip);
            targetSalRes = zeros(1,nProf);         
    end
  end

  profYearDayRatio = fromJulDayToYearDay(profJulDayAggrSel) ./ yearLength(profJulDayAggrSel);

  for iProf = 1:nProf
      if ~mod(iProf, floor(nProf/20))
          disp([int2str(iProf), '/', int2str(nProf)]);
      end

      %{
      targetLat = -40;
      targetLon = 40;
      if ~(profLatAggrSel(iProf) < (targetLat+windowSize) && profLatAggrSel(iProf) > (targetLat-windowSize) ...
        && profLongAggrSel(iProf) < (targetLon+windowSize) && profLongAggrSel(iProf) > (targetLon-windowSize) )
          betaTemp = NaN(1, size(betaGrid, 2));
          continue;
      else
      %}
      betaTemp = betaGrid(iProf,:);
%      end
%{

      profLat = profLatAggrSelRounded(iProf);
      profLong = profLongAggrSelRounded(iProf);

      iLat = find(latGrid(1,:) == profLat);
      iLong = find(longGrid(:,1) == profLong);
%}
      
      switch meanTag
          case 'NoTrend'
              fun = @(x) betaTemp(1) ...
                 + betaTemp(2) .* sin(2*pi*1*profYearDayRatio(iProf))' + betaTemp(3) .* cos(2*pi*1*profYearDayRatio(iProf))' ...
                 + betaTemp(4) .* sin(2*pi*2*profYearDayRatio(iProf))' + betaTemp(5) .* cos(2*pi*2*profYearDayRatio(iProf))' ...
                 + betaTemp(6) .* sin(2*pi*3*profYearDayRatio(iProf))' + betaTemp(7) .* cos(2*pi*3*profYearDayRatio(iProf))' ...
                 + betaTemp(8) .* sin(2*pi*4*profYearDayRatio(iProf))' + betaTemp(9) .* cos(2*pi*4*profYearDayRatio(iProf))' ...
                 + betaTemp(10) .* sin(2*pi*5*profYearDayRatio(iProf))' + betaTemp(11) .* cos(2*pi*5*profYearDayRatio(iProf))' ...
                 + betaTemp(12) .* sin(2*pi*6*profYearDayRatio(iProf))' + betaTemp(13) .* cos(2*pi*6*profYearDayRatio(iProf))';% ...
                 %+ betaTemp(19) * (x - midJulDay) + betaTemp(20) * (x - midJulDay).^2;

%{
              fun = @(x,y) betaTemp(1) ...
                 + betaTemp(2) * sin(2*pi*1*y)' + betaTemp(3) * cos(2*pi*1*y)' ...
                 + betaTemp(4) * sin(2*pi*2*y)' + betaTemp(5) * cos(2*pi*2*y)' ...
                 + betaTemp(6) * sin(2*pi*3*y)' + betaTemp(7) * cos(2*pi*3*y)' ...
                 + betaTemp(8) * sin(2*pi*4*y)' + betaTemp(9) * cos(2*pi*4*y)' ...
                 + betaTemp(10) * sin(2*pi*5*y)' + betaTemp(11) * cos(2*pi*5*y)' ...
                 + betaTemp(12) * sin(2*pi*6*y)' + betaTemp(13) * cos(2*pi*6*y)';% ...
                 %+ betaTemp(19) * (x - midJulDay) + betaTemp(20) * (x - midJulDay).^2;
%}
          case 'NoTrendVelSeas3'
              fun = @(x) betaTemp(1) ...
                 + betaTemp(2) * sin(2*pi*1*profYearDayRatio(iProf))' + betaTemp(3) * cos(2*pi*1*profYearDayRatio(iProf))' ...
                 + betaTemp(4) * sin(2*pi*2*profYearDayRatio(iProf))' + betaTemp(5) * cos(2*pi*2*profYearDayRatio(iProf))' ...
                 + betaTemp(6) * sin(2*pi*3*profYearDayRatio(iProf))' + betaTemp(7) * cos(2*pi*3*profYearDayRatio(iProf))';% ...
                 %+ betaTemp(19) * (x - midJulDay) + betaTemp(20) * (x - midJulDay).^2;
          case 'Trend'
              fun = @(x) betaTemp(1) ...
                 + betaTemp(2) * sin(2*pi*1*fromJulDayToYearDay(x)./yearLength(x))' + betaTemp(3) * cos(2*pi*1*fromJulDayToYearDay(x)./yearLength(x))' ...
                 + betaTemp(4) * sin(2*pi*2*fromJulDayToYearDay(x)./yearLength(x))' + betaTemp(5) * cos(2*pi*2*fromJulDayToYearDay(x)./yearLength(x))' ...
                 + betaTemp(6) * sin(2*pi*3*fromJulDayToYearDay(x)./yearLength(x))' + betaTemp(7) * cos(2*pi*3*fromJulDayToYearDay(x)./yearLength(x))' ...
                 + betaTemp(8) * sin(2*pi*4*fromJulDayToYearDay(x)./yearLength(x))' + betaTemp(9) * cos(2*pi*4*fromJulDayToYearDay(x)./yearLength(x))' ...
                 + betaTemp(10) * sin(2*pi*5*fromJulDayToYearDay(x)./yearLength(x))' + betaTemp(11) * cos(2*pi*5*fromJulDayToYearDay(x)./yearLength(x))' ...
                 + betaTemp(12) * sin(2*pi*6*fromJulDayToYearDay(x)./yearLength(x))' + betaTemp(13) * cos(2*pi*6*fromJulDayToYearDay(x)./yearLength(x))' ...
                 + betaTemp(19) * (x - midJulDay);% + betaTemp(20) * (x - midJulDay).^2;
          case 'Trend2'
              fun = @(x) betaTemp(1) ...
                 + betaTemp(2) * sin(2*pi*1*fromJulDayToYearDay(x)./yearLength(x))' + betaTemp(3) * cos(2*pi*1*fromJulDayToYearDay(x)./yearLength(x))' ...
                 + betaTemp(4) * sin(2*pi*2*fromJulDayToYearDay(x)./yearLength(x))' + betaTemp(5) * cos(2*pi*2*fromJulDayToYearDay(x)./yearLength(x))' ...
                 + betaTemp(6) * sin(2*pi*3*fromJulDayToYearDay(x)./yearLength(x))' + betaTemp(7) * cos(2*pi*3*fromJulDayToYearDay(x)./yearLength(x))' ...
                 + betaTemp(8) * sin(2*pi*4*fromJulDayToYearDay(x)./yearLength(x))' + betaTemp(9) * cos(2*pi*4*fromJulDayToYearDay(x)./yearLength(x))' ...
                 + betaTemp(10) * sin(2*pi*5*fromJulDayToYearDay(x)./yearLength(x))' + betaTemp(11) * cos(2*pi*5*fromJulDayToYearDay(x)./yearLength(x))' ...
                 + betaTemp(12) * sin(2*pi*6*fromJulDayToYearDay(x)./yearLength(x))' + betaTemp(13) * cos(2*pi*6*fromJulDayToYearDay(x)./yearLength(x))' ...
                 + betaTemp(19) * (x - midJulDay) + betaTemp(20) * (x - midJulDay).^2;
      end

      ResponseProfHat = fun(profJulDayAggrSel(iProf));%, profYearDayRatio(iProf));

      if is2step
        switch typeTag
          case 'lat'
            latFluxRes(iProf) = profLatFluxAggrSel(iProf) - ResponseProfHat;
          case 'lon'
            lonFluxRes(iProf) = profLonFluxAggrSel(iProf) - ResponseProfHat;
          otherwise
            intFluxRes(iProf) = intFluxProf(iProf) - ResponseProfHat;
        end
      else
        switch responseTag
            case 'Temp'
              targetTempRes(iProf) = targetTempProfPchip(iProf) - ResponseProfHat;
            case 'Dens'
              intDensRes(iProf) = targetDynhProf(iProf) - ResponseProfHat;
            case 'Sal'
              targetSalRes(iProf) = targetSalProfPchip(iProf) - ResponseProfHat;
        end
      end
  end

  if is2step
    switch typeTag
      case 'lat'
          if isStandardize
            stdRes = std(latFluxRes, 'omitnan');
            latFluxRes = latFluxRes ./ stdRes;
          end
          save(['./Data/',typeTag, fluxType,responseTag,'Res',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'_w',num2str(windowSize),'_Eq',num2str(eqBorder),'.mat'],...
      'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','latFluxRes', 'stdRes');
      case 'lon'
          if isStandardize
            stdRes = std(lonFluxRes, 'omitnan');
            lonFluxRes = lonFluxRes ./ stdRes;
          end
          save(['./Data/',typeTag, fluxType,responseTag,'Res',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'_w',num2str(windowSize),'_Eq',num2str(eqBorder),'.mat'],...
      'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','lonFluxRes', 'stdRes');
      otherwise
          if isStandardize
            stdRes = std(intFluxRes, 'omitnan');
            intFluxRes = intFluxRes ./ stdRes;
          end
          save(['./Data/',typeTag, fluxType,responseTag,'Res',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'_w',num2str(windowSize),'_Eq',num2str(eqBorder),'.mat'],...
      'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','intFluxRes', 'stdRes');
    end
  else
    switch responseTag
        case 'Temp'
            if isStandardize
              stdRes = std(targetTempRes, 'omitnan');
              targetTempRes = targetTempRes ./ stdRes;
            end
%            save(['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'.mat'],...
%        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','targetTempRes', 'stdRes');
        case 'Dens'
            if isStandardize
              stdRes = std(intDensRes, 'omitnan');
              intDensRes = intDensRes ./ stdRes;
            end
            save(['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'_w',num2str(windowSize),'.mat'],...
        'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel','intDensRes', 'stdRes');
        case 'Sal'
            if isStandardize
              stdRes = std(targetSalRes, 'omitnan');
              targetSalRes = targetSalRes ./ stdRes;
            end
%            save(['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'.mat'],...
%        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel','targetSalRes', 'stdRes');          
    end
  end

  if isPlot
    % Create Folder
    if is2step
      destFolder = [typeTag,fluxType,responseTag,verticalSelection,dataYear,'_Eq',num2str(eqBorder),profileTag];
    else
      destFolder = [typeTag,responseTag,verticalSelection,dataYear,meanTag,profileTag];
    end
    mkdir(['./Figures/',destFolder])

    %% Plot integrated conservative temperatures
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

    if is2step
        % Border Mask
        mask = ~(profLatAggrSel < eqBorder & profLatAggrSel > - eqBorder) .* 1;
        mask(mask == 0) = NaN;
        
%        targetLat = -40;
%        targetLon = 40;
      
%        boxMask = ~(profLatAggrSel < (targetLat+windowSize) && profLatAggrSel > (targetLat-windowSize) ...
%        && profLongAggrSel < (targetLon+windowSize) && profLongAggrSel > (targetLon-windowSize) );
%        boxMask(boxMask == 0) = NaN;
        
%        mask = mask.*boxMask;

        switch typeTag
            case 'lat'
                scatterm(profLatAggrSel,profLongAggrSel,10,mask.*profLatFluxAggrSel,'.');
                caxis([quantile(profLatFluxAggrSel(:), 0.01), quantile(profLatFluxAggrSel(:), 0.99)]);
            case 'lon'
                scatterm(profLatAggrSel,profLongAggrSel,10,mask.*profLonFluxAggrSel,'.');
                caxis([quantile(profLonFluxAggrSel(:), 0.01), quantile(profLonFluxAggrSel(:), 0.99)]);
            otherwise
                scatterm(profLatAggrSel,profLongAggrSel,10,mask.*intFluxProf,'.');
                caxis([quantile(intFluxProf(:), 0.01), quantile(intFluxProf(:), 0.99)]);
        end
        cLims = caxis;
        colormap(darkb2r(cLims(1),cLims(2)));
    else
        scatterm(profLatAggrSel,profLongAggrSel,10,targetDynhProf,'.');
        %scatterm(profLatAggrSel,profLongAggrSel,10,targetTempProfPchip,'.');    
        caxis([quantile(targetDynhProf(:), 0.01), quantile(targetDynhProf(:), 0.99)]);
        colormap(parula(100));
    end

    colorbar;

    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    if is2step
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,fluxType,responseTag,verticalSelection,dataYear,'Prof.eps']);
    else
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,responseTag,verticalSelection,dataYear,'Prof.eps']);
    end
      


    %% Plot fitted integrated conservative temperatures
    if is2step
      switch typeTag
        case 'lat'
          ResponseProfHat = profLatFluxAggrSel - stdRes.*latFluxRes;
        case 'lon'
          ResponseProfHat = profLonFluxAggrSel - stdRes.*lonFluxRes;
        otherwise
          ResponseProfHat = intFluxProf - stdRes.*intFluxRes;
      end
      ResponseProfHat = mask .* ResponseProfHat;
    else
      switch responseTag
          case 'Temp'
            ResponseProfHat = targetTempProfPchip - stdRes.*targetTempRes;
          case 'Dens'
            ResponseProfHat = targetDynhProf - stdRes.*intDensRes;
      end
    end

    % Random permutation of the plotting order
    %rng(12345);
    figure;
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');

    load coast;
    plotm(lat,long,'k');
    scatterm(profLatAggrSel,profLongAggrSel,10,ResponseProfHat,'.');    
    caxis([quantile(ResponseProfHat(:), 0.01), quantile(ResponseProfHat(:), 0.99)]);

    if is2step
      cLims = caxis;
      colormap(darkb2r(cLims(1),cLims(2)));
    else
        colormap(parula(100));
    end
    colorbar;

    %colormap(darkb2r(cLims(1),cLims(2)));
    %caxis([1.844e06,1.848e06]);

    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    if is2step
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,fluxType,responseTag,verticalSelection,dataYear,'ProfHat_',num2str(windowSize),'.eps']);
    else
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,responseTag,verticalSelection,dataYear,'ProfHat_',num2str(windowSize),'.eps']);
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

    colormap(parula(100));
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
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,fluxType,responseTag,verticalSelection,dataYear,'Res_',num2str(windowSize),'.eps']);
    else
      print('-depsc2',['./Figures/',destFolder,'/',typeTag,responseTag,verticalSelection,dataYear,'Res_',num2str(windowSize),'.eps']);
    end


    if exist('keepRatio', 'var')
        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

            % Border Mask         
    %        targetLat = -40;
    %        targetLon = 40;
          
    %        boxMask = ~(profLatAggrSel < (targetLat+windowSize) && profLatAggrSel > (targetLat-windowSize) ...
    %        && profLongAggrSel < (targetLon+windowSize) && profLongAggrSel > (targetLon-windowSize) );
    %        boxMask(boxMask == 0) = NaN;
            
    %        mask = mask.*boxMask;

            scatterm(profLatAggrSel,profLongAggrSel,10,keepRatio,'.');
            caxis([quantile(keepRatio(:), 0.01), quantile(keepRatio(:), 0.99)]);
            colormap(parula(100));

        colorbar;

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        print('-depsc2',['./Figures/',destFolder,'/',typeTag,responseTag,verticalSelection,dataYear,'YFRatio_',num2str(windowSize),'.eps']);
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
