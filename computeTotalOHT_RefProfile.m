function computeTotalOHT_RefProfile(responseTag, verticalSelection, targetPres, intStart, dataYear, windowSize, minNumberOfObs, meanTag, kernelType, targetVar, isPlot, isAdjusted, isAbsolute)
%% Compute the in-situ OHT used for is2step
%% For Reference velocity = 0

    refPres = 900
    
    tag = 'PchipPotTemp';
    if isAdjusted
        adjustTag = 'Adjusted';
    else
        adjustTag = [];
    end

    if isAbsolute
        absoluteTag = 'Absolute';
    else
        absoluteTag = [];
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

    [latGrid,longGrid] = meshgrid(linspace(-89.5,89.5,180), linspace(20.5,379.5,360));

    %% Load Fitted Mean and Anomaly
    % MeanField itself is not adjusted: will be added later
%{
    load(['./Results/Step1_',num2str(min(targetPres)),'_',num2str(max(targetPres)),'/meanField',responseTag,meanTag,tag,verticalSelection,dataYear,'_w',num2str(windowSize),'_',num2str(minNumberOfObs),'.mat']);
    % Adjusted pred Anomaly
    load(['./Results/Step1_',num2str(min(targetPres)),'_',num2str(max(targetPres)),'/anomaly',responseTag,verticalSelection,adjustTag,'SeasonSpaceTime',kernelType,targetVar,'Deriv_Profile','.mat'], 'predRes');
%}

%{
    if isAdjusted
        % Fetch pre-adjusted mean Anomaly: meanPredGrid
        srcFolder = ['anomaly_','int',responseTag,'_w',num2str(windowSize),'_',kernelType,'_',verticalSelection];
%        load(['./Results/Step1_',num2str(min(targetPres)),'_',num2str(max(targetPres)),'/',srcFolder,'/MeanAnomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,targetVar,'Deriv','.mat']);
    end
%}


    if isAbsolute
        load(['../Misc/AGVA/','meanZonVelRef',num2str(refPres),'.mat'], 'meanRefZonVel');
        load(['../Misc/AGVA/','meanMerVelRef',num2str(refPres),'.mat'], 'meanRefMerVel');
    end

    % Conversion: Deg / m
%{
    switch targetVar
        case 'lat'
            distGrid = 1 ./ distance(latGrid - 0.5, longGrid, latGrid + 0.5, longGrid,...
                            referenceEllipsoid('GRS80', 'm'));
        case 'lon'
            distGrid = 1 ./ distance(latGrid, longGrid - 0.5, latGrid, longGrid + 0.5,...
                            referenceEllipsoid('GRS80', 'm'));
    end
%}

    nProf = length(profLongAggrSel);
    profLatAggrSelRounded = roundHalf(profLatAggrSel);
    profLongAggrSelRounded = roundHalf(profLongAggrSel);
    profDerivSel = zeros(1, nProf);

    if isAbsolute
        profMeanRefVel = zeros(1, nProf);
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
        %{
        % Predict meanField from the nearest grid
        switch targetVar
            case 'lat'
                betaDeriv = betaGrid(iLong,iLat,14) +...
                            (profLongAggrSel(iProf) - predLong) .* betaGrid(iLong,iLat,16) + ...
                            2 .* (profLatAggrSel(iProf) - predLat) .* betaGrid(iLong,iLat,17);
            case 'lon'
                betaDeriv = betaGrid(iLong,iLat,15) +...
                            (profLatAggrSel(iProf) - predLat) .* betaGrid(iLong,iLat,16) + ...
                            2 .* (profLongAggrSel(iProf) - predLong) .* betaGrid(iLong,iLat,18);
        end
        if isAdjusted
            betaDeriv = betaDeriv + meanPredGrid(iLong, iLat);
        end
        
        profDerivSel(iProf) = (betaDeriv + predRes(iProf)) .* distGrid(iLong, iLat);
%}

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
        fprintf("Drop %d NaN Observations\n", sum(drop))
    end
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
    profVelAggr = zeros(1, numel(profDerivSel));
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
%            profLatMassFluxAggrInt = - profVelAggr; % Mass flux will be integrated from DerivSel

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

    end

end
