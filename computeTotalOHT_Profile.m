function computeTotalOHT_Profile(responseTag, verticalSelection, targetPres, intStart, dataYear, windowType, windowSize, minNumberOfObs, meanTag, kernelType, month, targetVar, isPlot, isAdjusted, isAbsolute, nAdjust, iterEM)
%% Compute the in-situ OHT used for is2step

    refPres = 900

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
      adjustPrevNumTag = [];
%      windowSizeTag = windowSizeFullTag;
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

    switch meanTag
        case 'NoTrend'
            nHar = 6;
            iLatCoef = 1 + (2*nHar) + 1;
            iLongCoef = iLatCoef + 1;
            iLatLongCoef = iLongCoef + 1;
            iLat2Coef = iLatLongCoef + 1;
            iLong2Coef = iLat2Coef + 1;
        case 'NoTrendVelSeas3'
            nHar = 3;
            iLatCoef = 1 + (2*nHar) + 1;
            iLongCoef = iLatCoef + 1;
            iLatCoefHar = iLatCoef + 4 + (1:(2*nHar));
            iLongCoefHar = iLatCoefHar + 2*nHar;
    end

    yearRange = strsplit(dataYear, '_');
    startYear = str2double(yearRange{2});
    endYear = str2double(yearRange{3});

    %% Load Target Temp Data
    nTargetPres = numel(targetPres);
    if nTargetPres > 1
        target = load(['./Data/targetTempProf',tag,num2str(min(targetPres)),'_',num2str(max(targetPres)),dataYear,'.mat']);
        % Additional load for density Data
        targetDens = load(['./Data/targetDensProf',tag,num2str(min(targetPres)),'_',num2str(max(targetPres)),dataYear,'.mat']);
        targetDensProf = targetDens.targetDensProf;
        clear targetDens;
    else
        target = load(['./Data/targetTempProf',tag,num2str(targetPres),dataYear,'.mat']);    
    end
    profLatAggrSel = target.profLatAggrSel;
    profLongAggrSel = target.profLongAggrSel;
    profJulDayAggrSel = target.profJulDayAggrSel;
    targetTempProf = target.targetTempProf;

    clear target;


    %% Load Fitted Mean and Anomaly
    % MeanField itself is not adjusted: will be added later
    load(['./Results/','meanField',responseTag,meanTag,tag,verticalSelection,dataYear,[],[],windowTypeTag,'_w',windowSizeTag,'_',num2str(minNumberOfObs),EMTag,'.mat']);
    % Adjusted pred Anomaly
    destFolder = ['anomaly_','int',responseTag,adjustNumTag,[],windowTypeTag,'_w',windowSizeFullTag,'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag];

    load(['./Results/',destFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,targetVar,'Deriv_Profile','.mat'], 'predRes');

    %% CHECK
    predRes(isnan(predRes)) = 0;

    if isAdjusted
        % Fetch pre-adjusted mean Anomaly: meanPredGrid
        srcFolder = ['anomaly_','int',responseTag,adjustPrevNumTag,[],windowTypeTag,'_w',windowSizeFullTag,'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag];
        load(['./Results/',srcFolder,'/MeanAnomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,targetVar,'Deriv','.mat'], 'meanPredGrid');
    else
        meanPredGrid = zeros(size(latGrid));
    end

    % ZONAL CHANGE
    if strcmp(targetVar, 'lat')
        predRes = - predRes;
        meanPredGrid = - meanPredGrid;
    end

    if isAbsolute
        load(['../Misc/AGVA/','meanZonVelRef',num2str(refPres),'.mat'], 'meanRefZonVel');
        load(['../Misc/AGVA/','meanMerVelRef',num2str(refPres),'.mat'], 'meanRefMerVel');
    end

    % Conversion: Deg / m
    switch targetVar
        case 'lat'
            distGrid = 1 ./ distance(latGrid - 0.5, longGrid, latGrid + 0.5, longGrid,...
                            referenceEllipsoid('WGS84', 'm'));
        case 'lon'
            distGrid = 1 ./ distance(latGrid, longGrid - 0.5, latGrid, longGrid + 0.5,...
                            referenceEllipsoid('WGS84', 'm'));
    end

    nProf = length(profLongAggrSel);
    profLatAggrSelRounded = roundHalf(profLatAggrSel);
    profLongAggrSelRounded = roundHalf(profLongAggrSel);
    profDerivSel = zeros(size(predRes));

    if isAbsolute
        profMeanRefVel = zeros(size(predRes));
    end

    latLinspace = latGrid(1,:);
    longLinspace = longGrid(:,1);
    sizeGrid = size(latGrid);

    for iProf = 1:nProf
        if ~mod(iProf, floor(nProf/20))
            disp([int2str(iProf), '/', int2str(nProf)]);
        end

        % Choose the nearest gridpoint for meanField prediction
        predLat = profLatAggrSelRounded(iProf);
        predLong = profLongAggrSelRounded(iProf);

        iLat = find(abs(latLinspace - predLat) < 1e-03);
        iLong = find(abs(longLinspace - predLong) < 1e-03);
    %    iGrid = sub2ind(sizeGrid, find(longLinspace == predLong), find(latLinspace == predLat));
        
        % Predict meanField from the nearest grid
        switch targetVar
            case 'lat'
                betaDeriv = betaGrid(iLong,iLat,iLatCoef) +...
                            (profLongAggrSel(iProf) - predLong) .* betaGrid(iLong,iLat,iLatLongCoef) + ...
                            2 .* (profLatAggrSel(iProf) - predLat) .* betaGrid(iLong,iLat,iLat2Coef);
                betaDeriv = - betaDeriv; % Zonal sign change
            case 'lon'
                betaDeriv = betaGrid(iLong,iLat,iLongCoef) +...
                            (profLatAggrSel(iProf) - predLat) .* betaGrid(iLong,iLat,iLatLongCoef) + ...
                            2 .* (profLongAggrSel(iProf) - predLong) .* betaGrid(iLong,iLat,iLong2Coef);
        end
        if isAdjusted
            betaDeriv = betaDeriv + meanPredGrid(iLong, iLat);
        end
        
        profDerivSel(iProf) = (betaDeriv + predRes(iProf)) .* distGrid(iLong, iLat);

        % Precompute
        if isAbsolute
            switch targetVar
                case 'lat'
                    profMeanRefVel(iProf) = meanRefZonVel(iLong, iLat);
                case 'lon'
                    profMeanRefVel(iProf) = meanRefMerVel(iLong, iLat);
            end
        end
    end

    % filter out NaN whose observation in target but also in int DataMask
%    drop = isnan(profDerivSel);
    if isAbsolute
%        fprintf("Drop %d Additional from Abs Vel\n", sum(~drop & isnan(profMeanRefVel)));
%        drop = drop | isnan(profMeanRefVel);
        drop = isnan(profMeanRefVel);
        profMeanRefVel = profMeanRefVel(~drop);
    end
    fprintf("Drop %d NaN Observations\n", sum(drop))
    profLatAggrSel = profLatAggrSel(~drop);
    profLongAggrSel = profLongAggrSel(~drop);
    profJulDayAggrSel = profJulDayAggrSel(~drop);
    targetTempProf = targetTempProf(targetPres==intStart, ~drop);
    targetDensProf = targetDensProf(targetPres==intStart, ~drop);
    profDerivSel = profDerivSel(~drop);
    %profDerivSel(profLatAggrSel < (-58.5 + 5) & profLatAggrSel > (-58.5 - 5) & profLongAggrSel < (329.5 + 5) & profLongAggrSel > (329.5-5))

    %% Calculate flux
    %Convert Kelvin
    targetTempProf = targetTempProf + 273.16;
    profVelAggr = profDerivSel ./ gsw_f(profLatAggrSel);
    if isAbsolute
        profVelAggr = profVelAggr + profMeanRefVel;
    end

    % Inflation factor when integrating wrt pressure
    profScaleInt = 1 ./ targetDensProf ./ gsw_grav(profLatAggrSel, intStart);

    profHeatFluxAggr = targetTempProf .* targetDensProf .* profVelAggr;
    profMassFluxAggr = targetDensProf .* profVelAggr;

    % For integrating in pressure, do not need to multiply insitu density
    profVolFluxAggrInt = profVelAggr ./ targetDensProf; % volume
    profHeatFluxAggrInt = targetTempProf .* profVelAggr; % Heat
%            profMassFluxAggrInt = profVelAggr; % Mass flux will be integrated from DerivSel


    saveName = ['./Data/',targetVar,'FluxProf',verticalSelection,dataYear,adjustTag,absoluteTag,'.mat'];
    if isAbsolute
        save(saveName,...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
        'profHeatFluxAggr', 'profMassFluxAggr', 'profVelAggr',...
        'profHeatFluxAggrInt', 'profDerivSel', 'profVolFluxAggrInt', 'profScaleInt',...
        'intStart', 'isAdjusted', 'isAbsolute', 'profMeanRefVel');
    else
        save(saveName,...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
        'profHeatFluxAggr', 'profMassFluxAggr', 'profVelAggr',...
        'profHeatFluxAggrInt', 'profDerivSel', 'profVolFluxAggrInt','profScaleInt',...
        'intStart', 'isAdjusted', 'isAbsolute');
    end


    if isPlot
        % Save directory
        destFolder = [targetVar,'Flux',verticalSelection,dataYear,adjustTag,absoluteTag];
        mkdir(['./Figures/',destFolder])

   
        % Equitorial Mask
        mask = ~(profLatAggrSel < 1 & profLatAggrSel > - 1) .* 1;
        mask(mask == 0) = NaN;
        
        %% Plot Target Temp and Density
        % Temperature
        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');
       
        scatterm(profLatAggrSel,profLongAggrSel,10, targetTempProf,'.');
        caxis(quantile(targetTempProf, [0.01, 0.99]));

        cb = colorbar;
        cb.Label.String = 'Temperature [K]';
        cb.FontSize = 14;

        colormap(parula);

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))        
        print('-depsc2',['./Figures/',destFolder,'/targetTempProf',num2str(intStart),dataYear,'.eps']); %lon

        % Density
        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');
       
        scatterm(profLatAggrSel,profLongAggrSel,10, targetDensProf,'.');
        caxis(quantile(targetDensProf, [0.01, 0.99]));

        cb = colorbar;
        cb.Label.String = 'Density [kg/m]';
        cb.FontSize = 14;

        colormap(parula);

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))        
        print('-depsc2',['./Figures/',destFolder,'/targetDensProf',num2str(intStart),dataYear,'.eps']); %lon
        
        %% Plot Velocity and Flux
        idx = 1:(numel(profLatAggrSel)); %randperm(nProf);

        % Plot Velocity
        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

        profVelAggrSel = mask.*profVelAggr(idx);
        scatterm(profLatAggrSel(idx),profLongAggrSel(idx),10, profVelAggrSel(idx),'.');
        caxis([quantile(profVelAggrSel, 0.01), quantile(profVelAggrSel, 0.99)]);

        cb = colorbar;
        cb.Label.String = 'Velocity [m/s]';
        cb.FontSize = 14;

        cLims = caxis;
%        colormap(darkb2r(cLims(1),cLims(2)));

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        switch targetVar
            case 'lat'
                print('-depsc2',['./Figures/',destFolder,'/ZvelProf',num2str(intStart),dataYear,'.eps']); %lat
            case 'lon'
                print('-depsc2',['./Figures/',destFolder,'/MvelProf',num2str(intStart),dataYear,'.eps']); %lon
        end


        % Plot Mass Flux
        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');
        
%        profMtAggr = targetDensProf .* profVelAggr;
        switch targetVar
            case 'lat'
                profMtAggr = profLatMassFluxAggr;
            case 'lon'
                profMtAggr = profLonMassFluxAggr;
        end
        
        profMtAggr = mask.*profMtAggr(idx);
        scatterm(profLatAggrSel(idx),profLongAggrSel(idx),10, profMtAggr(idx),'.');
        caxis([quantile(profMtAggr, 0.02), quantile(profMtAggr, 0.98)]);
                
        cb = colorbar;
        cb.Label.String = 'Mass Flux [kg / s]';  
        cb.FontSize = 14;

        cLims = caxis;
%        colormap(darkb2r(cLims(1),cLims(2)));

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        switch targetVar
            case 'lat'
                print('-depsc2',['./Figures/',destFolder,'/ZmtProf',num2str(intStart),dataYear,'.eps']); %lat
            case 'lon'
                print('-depsc2',['./Figures/',destFolder,'/MmtProf',num2str(intStart),dataYear,'.eps']); %lon
        end
        
        % Plot Heat Flux
        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

%        profHtAggr = gsw_cp0 .* targetTempProf .* targetDensProf .* profVelAggr;
        switch targetVar
            case 'lat'
                profHtAggr = gsw_cp0 .* profLatHeatFluxAggr;
            case 'lon'
                profHtAggr = gsw_cp0 .* profLonHeatFluxAggr;
        end
        profHtAggr = mask.*profHtAggr(idx) .* 1e-15;

        scatterm(profLatAggrSel(idx),profLongAggrSel(idx),10, profHtAggr(idx),'.');
        caxis([quantile(profHtAggr, 0.02), quantile(profHtAggr, 0.98)]);

        cb = colorbar;
        cb.Label.String = 'Heat Flux [PW / m^3]';  
        cb.FontSize = 14;

        cLims = caxis;
%        colormap(darkb2r(cLims(1),cLims(2)));

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        switch targetVar
            case 'lat'
                print('-depsc2',['./Figures/',destFolder,'/ZhtProf',num2str(intStart),dataYear,'.eps']); %lat
            case 'lon'
                print('-depsc2',['./Figures/',destFolder,'/MhtProf',num2str(intStart),dataYear,'.eps']); %lon
        end

        % Plot Heat Flux / Mass ratio
        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

        scatterm(profLatAggrSel,profLongAggrSel,10,...
            profHtAggr ./ profMtAggr ./ mean(profHtAggr ./ profMtAggr, 'omitnan'),'.');
        caxis(quantile(profHtAggr ./ profMtAggr ./ mean(profHtAggr ./ profMtAggr, 'omitnan'), [0.02,0.98]));

        cb = colorbar;
        cb.Label.String = 'Ratio';  
        cb.FontSize = 14;

        cLims = caxis;
%        colormap(darkb2r(cLims(1),cLims(2)));

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        title('Heat / Mass')
        switch targetVar
            case 'lat'
                print('-depsc2',['./Figures/',destFolder,'/ZhtMtRatio',num2str(intStart),dataYear,'.eps']); %lat
            case 'lon'
                print('-depsc2',['./Figures/',destFolder,'/MhtMtRatio',num2str(intStart),dataYear,'.eps']); %lon
        end

         % Plot Heat Flux / Vol ratio
        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

        scatterm(profLatAggrSel,profLongAggrSel,10,...
            profHtAggr ./ profVelAggr ./ mean(profHtAggr ./ profVelAggr, 'omitnan'),'.');
        caxis(quantile(profHtAggr ./ profVelAggr ./ mean(profHtAggr ./ profVelAggr, 'omitnan'), [0.02,0.98]));

        cb = colorbar;
        cb.Label.String = 'Ratio';  
        cb.FontSize = 14;

        cLims = caxis;
%        colormap(darkb2r(cLims(1),cLims(2)));

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        title('Heat / Vol')
        switch targetVar
            case 'lat'
                print('-depsc2',['./Figures/',destFolder,'/ZhtVtRatio',num2str(intStart),dataYear,'.eps']); %lat
            case 'lon'
                print('-depsc2',['./Figures/',destFolder,'/MhtVtRatio',num2str(intStart),dataYear,'.eps']); %lon
        end        
    end

end
