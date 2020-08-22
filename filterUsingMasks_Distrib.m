function filterUsingMasks_Distrib(typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, isPlot, fluxType, isAdjusted, isAbsolute)
%% CURRENT: Just filter out the map mask... Not really necessary procedure anymore.
    tag = 'PchipPotTemp';

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

    if isempty(isAdjusted)
        isAdjusted = false;
    end
    if isempty(isAbsolute)
        isAbsolute = false;
    end

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
        data = load(['./Data/',typeTag,responseTag,'Prof',tag,presString,dataYear,adjustTag,absoluteTag,'.mat']);
    else
        presString = verticalSelection;
        switch responseTag
            case 'Flux'
                load(['./Data/',typeTag,responseTag,'Prof',verticalSelection,dataYear,adjustTag,absoluteTag,'.mat']);
            otherwise
                load(['./Data/',typeTag,responseTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'.mat']);
        end
    end

    intStart = data.intStart;
    intEnd = data.intEnd;

    % Load Mask
%{
    switch typeTag
        case 'target'
            if isnumeric(verticalSelection)
                load(['./Data/dataMask',typeTag,responseTag,presString,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat']);
            else
                load(['./Data/dataMask',typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat']);
            end

        case 'int'
            if isnumeric(verticalSelection)
                load(['./Data/dataMask',presString,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat']);
            else
                load(['./Data/dataMask',verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat']);            
            end
        otherwise %'lat' 'lon'
            if isnumeric(verticalSelection) %intlatlon
                load(['./Data/dataMask',typeTag,responseTag,presString,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat']);
            else
                load(['./Data/dataMask',typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat']);
            end
    end
%}


    maskJohn = ncread('./RG_climatology/RG_ArgoClim_Temperature_2016.nc','BATHYMETRY_MASK',[1 1 25],[Inf Inf 1]);
    maskJohn(maskJohn == 0) = 1;
    maskJohn = [NaN*ones(360,25) maskJohn NaN*ones(360,25)];

    mask = maskJohn; %.* dataMask;
    
    % Initialize
    profLatAggrRounded = roundHalf(data.profLatAggrSel);
    profLongAggrRounded = roundHalf(data.profLongAggrSel);
    nProf = length(profLatAggrRounded);
    keep = zeros(1,nProf);
    [latGrid, longGrid] = meshgrid(linspace(-89.5,89.5,180), linspace(20.5,379.5,360));

    % retrieve keeping index
    for iProf = 1:nProf
%        disp(iProf)
        latIdx = find(abs(latGrid(1,:) - profLatAggrRounded(iProf)) < 1e-03);
        longIdx = find(abs(longGrid(:,1) - profLongAggrRounded(iProf)) < 1e-03);
        keep(iProf) = ~isnan(mask(longIdx,latIdx));
    end
    keep = logical(keep);
    fprintf("%d / %d (%d) were kept\n", sum(keep), nProf, sum(keep)/nProf);
    
%{
    profLatAggrSel = profLatAggrSel(keep);
    profLongAggrSel = profLongAggrSel(keep);
    profJulDayAggrSel = profJulDayAggrSel(keep);
%}
    switch typeTag
        case 'target'
            profFloatIDAggrSel = profFloatIDAggrSel(keep);
            switch responseTag
                case 'Temp'
                    targetTempProf = targetTempProf(:,keep);
                    save(['./Data/',typeTag,responseTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat'],...
                        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
                        'targetTempProf','targetPres');
                case 'Sal'
                    targetSalProf = targetSalProf(:,keep);
                    save(['./Data/',typeTag,responseTag,'Prof',tag,verticalSelection,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat'],...
                        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
                        'targetSalProf','targetPres');
            end
        case 'int'
            if isnumeric(verticalSelection)                
                for presIdx = 1:nTargetPres
                    disp(['Pres: ',  num2str(verticalSelection(presIdx))]);
                    isRetain = ~isnan(data.targetDynhProf(presIdx, :));
                    fprintf("%d / %d (%d) were kept\n", sum(isRetain & keep), sum(isRetain), sum(isRetain & keep)/sum(isRetain));

                    isRetain = isRetain & keep;

                    profLatAggrSel = data.profLatAggrSel(isRetain);
                    profLongAggrSel = data.profLongAggrSel(isRetain);
                    profJulDayAggrSel = data.profJulDayAggrSel(isRetain);
                    profYearAggrSel = data.profYearAggrSel(isRetain);
                    profFloatIDAggrSel = data.profFloatIDAggrSel(isRetain);
                    
                    intTempProf = data.intTempProf(presIdx, isRetain);
                    intDensProf = data.intDensProf(presIdx, isRetain);
                    targetDynhProf = data.targetDynhProf(presIdx, isRetain);
                    saveName = ['./Data/intTempDensProf',tag,'Relative',num2str(targetPres(presIdx)),dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat'];
                    save(saveName,...
                        'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel',...
                        'intTempProf','intDensProf','targetDynhProf','isRetain','intStart','intEnd');
                end
            else
                intTempProf = intTempProf(keep);
                intDensProf = intDensProf(keep);
                targetDynhProf = targetDynhProf(keep);
                saveName = ['./Data/intTempDensProf',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat'];
                save(saveName,...
                    'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel',...
                    'intTempProf','intDensProf','targetDynhProf','intStart','intEnd');
            end
        % Flux
        case 'lat'
            switch fluxType
                case 'heat'
                    profLatFluxAggrSel = profLatFluxAggr(keep);
                case 'mass'
                    profLatFluxAggrSel = profLatMassFluxAggr(keep);
                case 'vol'
                    profLatFluxAggrSel = - profDerivSel(keep) ./ gsw_f(profLatAggrSel);
            end
            save(['./Data/',typeTag,fluxType,responseTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat'],...
                'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
                'profLatFluxAggrSel','intStart');
        case 'lon'
            switch fluxType
                case 'heat'
                    profLonFluxAggrSel = profLonFluxAggr(keep);
                case 'mass'
                    profLonFluxAggrSel = profLonMassFluxAggr(keep);
                case 'vol'
                    profLonFluxAggrSel = profDerivSel(keep) ./ gsw_f(profLatAggrSel);
            end
            save(['./Data/',typeTag,fluxType,responseTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat'],...
                'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
                'profLonFluxAggrSel','intStart');
        % intFlux
        otherwise % intlatlon
            % Beware of the choice either Exact(take care of grav) or not
            switch fluxType
                case 'heat'
%                    intFluxProf = intFluxProf(keep);
                    intFluxProf = intFluxExactProf(keep);
                case 'mass'
%                    intFluxProf = intMassFluxProf(keep);
                    intFluxProf = intMassFluxExactProf(keep);
                case 'vol'
%                    intFluxProf = intVolFluxProf(keep);
                    intFluxProf = intVolFluxExactProf(keep);
            end
            save(['./Data/',typeTag,fluxType,responseTag,'Prof',tag,num2str(min(targetPres)),'_',num2str(max(targetPres)),dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSize),'.mat'],...
                'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
                'intFluxProf','intStart','intEnd');
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