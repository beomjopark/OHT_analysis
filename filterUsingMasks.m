function filterUsingMasks(typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, isPlot, fluxType, isAdjusted, isAbsolute)
    tag = 'PchipPotTemp';

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

    if strcmp(windowType, 'spherical')
    % Determine reference distance
    refDist = distance(0, 180, 0+windowSizeMean, 180,...
                            referenceEllipsoid('WGS84', 'm'))
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
    
    % Load Data
    if isnumeric(verticalSelection)
        % Specifically for target / intlatlon
        % Multiple pressure inputs
        targetPres = verticalSelection;
        nTargetPres = numel(targetPres);
        if nTargetPres > 1
            presString = [num2str(min(targetPres)),'_',num2str(max(targetPres))];
        else
            presString = num2str(targetPres);
        end

        if strcmp(typeTag, 'int')
            presString = ['Relative', presString];
        end
        load(['./Data/',typeTag,responseTag,'Prof',tag,presString,dataYear,adjustTag,absoluteTag,'.mat']);
    else
        presString = verticalSelection;
        switch responseTag
            case 'Flux'
                load(['./Data/',typeTag,responseTag,'Prof',verticalSelection,dataYear,adjustTag,absoluteTag,'.mat']);
                %% FILTER NAN count
                drop = isnan(profDerivSel);
                profLatAggrSel = profLatAggrSel(~drop);
                profLongAggrSel = profLongAggrSel(~drop);
                profJulDayAggrSel = profJulDayAggrSel(~drop);  

                profHeatFluxAggrInt = profHeatFluxAggrInt(~drop);
                profMassFluxAggr = profMassFluxAggr(~drop);
                profVelAggr = profVelAggr(~drop);
            case {'Temp', 'Dens'}
                srcName = ['./Data/',typeTag,'TempDens','Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'.mat']
                load(srcName);                
            otherwise
                load(['./Data/',typeTag,responseTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'.mat']);
        end
    end

    % Load Mask
    switch typeTag
        case 'target'
            if isnumeric(verticalSelection)
                load(['./Data/dataMask',typeTag,responseTag,presString,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);
            else
                load(['./Data/dataMask',typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);
            end

        case 'int'
            if isnumeric(verticalSelection)
                load(['./Data/dataMask',presString,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);
            else
                load(['./Data/dataMask',verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);            
            end
        otherwise %'lat' 'lon'
            if isnumeric(verticalSelection) %intlatlon
                load(['./Data/dataMask',typeTag,responseTag,presString,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);
            else
                load(['./Data/dataMask',typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);
            end
    end
    maskJohn = ncread('./RG_climatology/RG_ArgoClim_Temperature_2016.nc','BATHYMETRY_MASK',[1 1 25],[Inf Inf 1]);
    maskJohn(maskJohn == 0) = 1;
    maskJohn = [NaN*ones(360,25) maskJohn NaN*ones(360,25)];

    mask = maskJohn .* dataMask;

    % Initialize
    profLatAggrRounded = roundHalf(profLatAggrSel);
    profLongAggrRounded = roundHalf(profLongAggrSel);
    nProf = length(profLatAggrRounded);
    keep = zeros(1,nProf);

    % retrieve keeping index
    for iProf = 1:nProf
        latIdx = find(abs(latGrid(1,:) - profLatAggrRounded(iProf)) < 1e-03);
        longIdx = find(abs(longGrid(:,1) - profLongAggrRounded(iProf)) < 1e-03);    
        keep(iProf) = ~isnan(mask(longIdx,latIdx));
    end
    fprintf("%d / %d (%d) were kept\n", sum(keep), nProf, sum(keep)/nProf);

    keep = logical(keep);
    profLatAggrSel = profLatAggrSel(keep);
    profLongAggrSel = profLongAggrSel(keep);
    profJulDayAggrSel = profJulDayAggrSel(keep);
    switch typeTag
        case 'target'
            profFloatIDAggrSel = profFloatIDAggrSel(keep);
            switch responseTag
                case 'Temp'
                    targetTempProf = targetTempProf(:,keep);
                    save(['./Data/',typeTag,responseTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat'],...
                        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
                        'targetTempProf','targetPres');
                case 'Sal'
                    targetSalProf = targetSalProf(:,keep);
                    save(['./Data/',typeTag,responseTag,'Prof',tag,verticalSelection,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat'],...
                        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
                        'targetSalProf','targetPres');
            end
        case 'int'
            profFloatIDAggrSel = profFloatIDAggrSel(keep);

            if isnumeric(verticalSelection)
                intTempProfTemp = intTempProf;
                intDensProfTemp = intDensProf;
                targetDynhProfTemp = targetDynhProf;
                for presIdx = 1:nTargetPres
                    intTempProf = intTempProfTemp(presIdx, keep);
                    intDensProf = intDensProfTemp(presIdx, keep);
                    targetDynhProf = targetDynhProfTemp(presIdx, keep);
                    saveName = ['./Data/intTempDensProf',tag,'Relative',num2str(targetPres(presIdx)),dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat'];
                    save(saveName,...
                        'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel',...
                        'intTempProf','intDensProf','targetDynhProf','intStart','intEnd');
                end
            else
                intTempProf = intTempProf(keep);
                intDensProf = intDensProf(keep);
                targetDynhProf = targetDynhProf(keep);
                saveName = ['./Data/intTempDensProf',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat'];
                save(saveName,...
                    'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel',...
                    'intTempProf','intDensProf','targetDynhProf','intStart','intEnd');
            end
        % Flux            
        case {'intlat', 'intlon'} % intlatlon
            % Beware of the choice either Exact(take care of grav) or not
            switch fluxType
                case 'heat'
%                    intFluxProf = intFluxProf(keep);
                    intFluxProf = intHeatFluxExactProf(keep);
                case 'mass'
%                    intFluxProf = intMassFluxProf(keep);
                    intFluxProf = intMassFluxExactProf(keep);
                case 'vol'
%                    intFluxProf = intVolFluxProf(keep);
                    intFluxProf = intVolFluxExactProf(keep);
            end
            save(['./Data/',typeTag,fluxType,responseTag,'Prof',tag,num2str(min(targetPres)),'_',num2str(max(targetPres)),dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat'],...
                'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
                'intFluxProf','intStart','intEnd');
        otherwise %{'lat', 'lon'}
            switch fluxType
                case 'heat'
                    profFluxAggrSel = profHeatFluxAggrInt(keep);
                case 'mass'
                    profFluxAggrSel = profMassFluxAggr(keep);
                case 'vol'
                    profFluxAggrSel = profVelAggr(keep);
            end
            save(['./Data/',typeTag,fluxType,responseTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat'],...
                'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
                'profFluxAggrSel','intStart');
        % intFlux
    end


    if isPlot % Deprecated: See subtractMeanField
        fprintf('Deprecated: See subtractMeanField\n');
        %{
        nProf = length(intTempProf);
        plotIdx = 1:nProf; % randperm(nProf);

        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

        scatterm(profLatAggrSel(plotIdx),profLongAggrSel(plotIdx),10,intTempProf(plotIdx),'.');

        colormap(jet(100));
        colorbar;
        caxis([-1e3,3e4]);

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        print('-depsc2','./Figures/intTempProfFiltered.eps');


        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

        scatterm(profLatAggrSel(plotIdx),profLongAggrSel(plotIdx),10,intDensProf(plotIdx),'.');

        colormap(jet(100));
        colorbar;
        caxis([1.844e06,1.848e06]);

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        print('-depsc2','./Figures/intDensProfFiltered.eps');
        %}
    end
end
