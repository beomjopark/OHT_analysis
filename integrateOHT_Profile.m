function integrateOHT_Profile(intStartList, dataYear, minNumberOfObs, targetVar, isPlot, isAdjusted, isAbsolute)
%% Integrate the Total OHT Profile
%% Plot Random profiles
%% Does not return object but save mat file

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

    %% Load data
    nTargetPres = numel(intStartList);

    % Load first target to grab setting
    intStart = intStartList(1);
    verticalSelection = strcat('Relative', num2str(intStart));
    targetFlux = load(['./Data/',targetVar,'FluxProf',verticalSelection,dataYear,adjustTag,absoluteTag,'.mat']);
    profLatAggrSel    = targetFlux.profLatAggrSel;
    profLongAggrSel   = targetFlux.profLongAggrSel;
    profJulDayAggrSel = targetFlux.profJulDayAggrSel;

    % Get Year for compatability
    nProfile = numel(profLatAggrSel);
    profYearAggrSel = zeros(size(profJulDayAggrSel));
    for iProf = 1:nProfile
        temp = datevec(profJulDayAggrSel(iProf));
        profYearAggrSel(iProf) = temp(1);
    end

    % Initialize Flux for merging
    profScaleIntSel = zeros(nTargetPres, nProfile);    
    profHeatFluxAggrSel = zeros(nTargetPres, nProfile);
    profDerivAggrSel = zeros(nTargetPres, nProfile);
    profVolFluxAggrSel = zeros(nTargetPres, nProfile);
    profHeatFluxAggrSel(1, :) = targetFlux.profHeatFluxAggrInt;
    profVolFluxAggrSel(1, :) = targetFlux.profVolFluxAggrInt;
    profDerivAggrSel(1, :) = targetFlux.profDerivSel;
    profScaleIntSel(1, :) = 1 ./ gsw_grav(targetFlux.profLatAggrSel, intStart);%targetFlux.profScaleInt;

    % Load Leftout targets to merge Flux
    for presIdx = 2:nTargetPres
        intStart = intStartList(presIdx);
        verticalSelection = strcat('Relative', num2str(intStart));
        targetFlux = load(['./Data/',targetVar,'FluxProf',verticalSelection,dataYear,adjustTag,absoluteTag,'.mat']);
        % Sanity Check: Merging the matching float for each pressure level
        if ~(all(profLatAggrSel == targetFlux.profLatAggrSel) && all(profLongAggrSel == targetFlux.profLongAggrSel) && all(profJulDayAggrSel == targetFlux.profJulDayAggrSel))
            error('Unmatching Float exists!');
        end
        
        profHeatFluxAggrSel(presIdx, :) = targetFlux.profHeatFluxAggrInt;
        profVolFluxAggrSel(presIdx, :) = targetFlux.profVolFluxAggrInt;
        profDerivAggrSel(presIdx, :) = targetFlux.profDerivSel;
        profScaleIntSel(presIdx, :) = 1 ./ gsw_grav(targetFlux.profLatAggrSel, intStart);%targetFlux.profScaleInt;
    end
    
    %% Filter NaNs
    isDrop = (sum(isnan(profDerivAggrSel)) > 0);
    profLatAggrSel    = profLatAggrSel(~isDrop);
    profLongAggrSel   = profLongAggrSel(~isDrop);
    profJulDayAggrSel = profJulDayAggrSel(~isDrop);
    profYearAggrSel = profYearAggrSel(~isDrop);
    profHeatFluxAggrSel = profHeatFluxAggrSel(:, ~isDrop);
    profVolFluxAggrSel = profVolFluxAggrSel(:, ~isDrop);
    profDerivAggrSel = profDerivAggrSel(:, ~isDrop);
    profScaleIntSel = profScaleIntSel(:, ~isDrop);
    
    nProfile = numel(profLatAggrSel);

    %% Vertical Integration with PCHIP Interpolation
    intStart = min(intStartList);
    intEnd = max(intStartList);
    gridx = linspace(intStart, intEnd, 5000); % Fine grid, takes about 13 mins to run
    intHeatFluxProf = zeros(1, nProfile); % Integrated Flux postmultiply gravity
    intHeatFluxExactProf = zeros(1, nProfile); % Integrated Flux premultiply gravity

    intDerivProf = zeros(1, nProfile); % Integrated Deriv (for mass Flux)
    intDerivExactProf = zeros(1, nProfile); % Integrated Deriv (for mass Flux)

    intVolFluxProf = zeros(1, nProfile); % Integrated Volume Flux postmultiply gravity
    intVolFluxExactProf = zeros(1, nProfile); % Integrated Volume Flux premultiply gravity

    intScaleProf = zeros(1, nProfile);

    tic;
    parfor iGrid = 1:nProfile
        % Interpolate grids with pchip
        pchip_int_flux=interp1(intStartList, profHeatFluxAggrSel(:,iGrid), gridx, 'pchip');
        pchip_int_deriv=interp1(intStartList, profDerivAggrSel(:,iGrid), gridx, 'pchip');
        pchip_int_volflux=interp1(intStartList, profVolFluxAggrSel(:,iGrid), gridx, 'pchip');

        pchip_int_flux_exact=interp1(intStartList, profHeatFluxAggrSel(:,iGrid) ./ gsw_grav(profLatAggrSel(iGrid), intStartList)', gridx, 'pchip');
        pchip_int_deriv_exact=interp1(intStartList, profDerivAggrSel(:,iGrid)./ gsw_grav(profLatAggrSel(iGrid), intStartList)', gridx, 'pchip');
        pchip_int_volflux_exact=interp1(intStartList, profVolFluxAggrSel(:,iGrid) ./ gsw_grav(profLatAggrSel(iGrid), intStartList)', gridx, 'pchip');

        pchip_int_scale=interp1(intStartList, profScaleIntSel(:,iGrid), gridx, 'pchip');

        % Integrate with trapizoidal rule
        intHeatFluxProf(iGrid)  = trapz(gridx, pchip_int_flux);
        intDerivProf(iGrid)  = trapz(gridx, pchip_int_deriv);
        intVolFluxProf(iGrid)  = trapz(gridx, pchip_int_volflux);

        intHeatFluxExactProf(iGrid)  = trapz(gridx, pchip_int_flux_exact);
        intDerivExactProf(iGrid) = trapz(gridx, pchip_int_deriv_exact); 
        intVolFluxExactProf(iGrid)  = trapz(gridx, pchip_int_volflux_exact);

        intScaleProf(iGrid)  = trapz(gridx, pchip_int_scale);
    end
    toc;

    % Mass Transport
    intMassFluxProf = intDerivProf ./ gsw_f(profLatAggrSel);
    intMassFluxExactProf = intDerivExactProf ./ gsw_f(profLatAggrSel);

    % Analytic solution of integrated reference velocity
    switch targetVar
        case 'lat'
            intMassFluxProf = - intMassFluxProf;
            intMassFluxExactProf = - intMassFluxExactProf;
    end

    if isAbsolute
        intMassFluxProf = intMassFluxProf + (targetFlux.profMeanRefVel(~isDrop) .* intScaleProf);
        intMassFluxExactProf = intMassFluxExactProf + (targetFlux.profMeanRefVel(~isDrop) .* intScaleProf);
    end

    save(['./Data/int',targetVar,'FluxProf',tag,num2str(min(intStartList)),'_',num2str(max(intStartList)),dataYear,adjustTag,absoluteTag,'.mat'], ...
        'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel', ...
        'intHeatFluxProf', 'intHeatFluxExactProf', 'intMassFluxProf', 'intMassFluxExactProf', 'intVolFluxProf', 'intVolFluxExactProf',...
         'intStart','intEnd', 'isAdjusted', 'isAbsolute');


    if isPlot
        %% 1. Plot OHT Profiles
        destFolder = ['int',targetVar,'Flux',num2str(min(intStartList)),'_',num2str(max(intStartList)),dataYear,adjustTag,absoluteTag];
        mkdir(['./Figures/',destFolder])
        mkdir(['./Figures/',destFolder,'/Profile'])

        rng(12345);
        nPlotFloat = 20;
        plotIdxList = randperm(nProfile, nPlotFloat);
        switch targetVar
            case 'lat'
                profMassFluxAggrSelInt = - profDerivAggrSel ./ gsw_f(profLatAggrSel);
                if isAbsolute
                    profMassFluxAggrSelInt = profMassFluxAggrSelInt + targetFlux.profMeanRefVel(~isDrop);
                end
            case 'lon'
                profMassFluxAggrSelInt = profDerivAggrSel ./ gsw_f(profLatAggrSel);
                if isAbsolute
                    profMassFluxAggrSelInt = profMassFluxAggrSelInt + targetFlux.profMeanRefVel(~isDrop);
                end
        end

        for plotIdx = plotIdxList
            plotPchipFlux = interp1(intStartList, profMassFluxAggrSelInt(:,plotIdx), gridx, 'pchip');

            figHandle = figure;
            plot(profMassFluxAggrSelInt(:, plotIdx), intStartList, 'Marker','+',...
                'Color', 'none',...
                'MarkerSize', 10, 'MarkerEdgeColor', 'black')
            hold on
            plot(plotPchipFlux, gridx, '-k',...
                'Color', [0.8 0.8 0.8],...
                'MarkerSize', 10, 'MarkerEdgeColor', 'black')
            hold off

            switch targetVar
                case 'lat'
                    xlabel('Zonal Velocity','FontSize',14);
                case 'lon'
                    xlabel('Meridional Velocity','FontSize',14);
            end
            ylabel('Pressure (dBar)','FontSize',14);
            ylim([intStart - 20,  max(intEnd) + 20]);

            set(gca,'xaxisLocation','top');
            set(gca, 'YDir','reverse');

            str = sprintf('Lat %4.2f%c, Lon %4.2f%c',...
                            profLatAggrSel(plotIdx), char(176),...
                            profLongAggrSel(plotIdx), char(176));
            annotation('textbox',[0.55 0 0.3 0.1],'String',str,'FitBoxToText','on','EdgeColor','none');

            set(gcf,'units','centimeters')
            set(gcf,'pos',[0 0 11 11])
            set(gcf,'paperunits',get(gcf,'units')) 
            set(gcf,'paperpos',get(gcf,'pos'))
            print('-depsc2',['./Figures/',destFolder,'/Profile/',...
                targetVar,'velProfile',num2str(plotIdx),'.eps']);
            close(figHandle);
        end

        rng(12345);
        nPlotFloat = 20;
        plotIdxList = randperm(nProfile, nPlotFloat);

        for plotIdx = plotIdxList
            plotPchipFlux = interp1(intStartList, profHeatFluxAggrSel(:,plotIdx), gridx, 'pchip');

            figHandle = figure;
            plot(profHeatFluxAggrSel(:, plotIdx), intStartList, 'Marker','+',...
                'Color', 'none',...
                'MarkerSize', 10, 'MarkerEdgeColor', 'black')
            hold on
            plot(plotPchipFlux, gridx, '-k',...
                'Color', [0.8 0.8 0.8],...
                'MarkerSize', 10, 'MarkerEdgeColor', 'black')
            hold off

            switch targetVar
                case 'lat'
                    xlabel('Temperature {\times} Zonal Velocity','FontSize',14);
                case 'lon'
                    xlabel('Temperature {\times} Meridional Velocity','FontSize',14);
            end
            ylabel('Pressure (dBar)','FontSize',14);
            ylim([intStart - 20,  max(intEnd) + 20]);

            set(gca,'xaxisLocation','top');
            set(gca, 'YDir','reverse');

            str = sprintf('Lat %4.2f%c, Lon %4.2f%c',...
                            profLatAggrSel(plotIdx), char(176),...
                            profLongAggrSel(plotIdx), char(176));
            annotation('textbox',[0.55 0 0.3 0.1],'String',str,'FitBoxToText','on','EdgeColor','none');

            set(gcf,'units','centimeters')
            set(gcf,'pos',[0 0 11 11])
            set(gcf,'paperunits',get(gcf,'units')) 
            set(gcf,'paperpos',get(gcf,'pos'))
            print('-depsc2',['./Figures/',destFolder,'/Profile/',...
                targetVar,'OHTProfile',num2str(plotIdx),'.eps']);
            close(figHandle);
        end
        
        %% 2. Plot Transport
        % Equitorial Mask
        mask = ~(profLatAggrSel < 1 & profLatAggrSel > - 1) .* 1;
        mask(mask == 0) = NaN;

        % Plot Vol Transport
        figHandle = figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

        profVoltAggr = mask.*intVolFluxExactProf.*1e-6;
%        profVoltAggr = mask.*intMassFluxProf;%./ gsw_grav(profLatAggrSel, mean(intStartList));
        scatterm(profLatAggrSel,profLongAggrSel,10, profVoltAggr,'.');
        caxis(quantile(profVoltAggr, [0.01, 0.99]));

        cb = colorbar;
        cb.Label.String = 'Volume Transport [Sv]';  
        cb.FontSize = 14;

        cLims = caxis;
        colormap(darkb2r(cLims(1),cLims(2)));

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        %set(gcf,'pos',[0 0 22.5 15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        switch targetVar
            case 'lat'
                title(['Zonal Volume Transport ',num2str(intStart), '~', num2str(intEnd),'dBar']);
                print('-depsc2',['./Figures/',destFolder,'/ZvtProf',num2str(intStart),'_',num2str(intEnd),dataYear,'.eps']); %lat
            case 'lon'
                title(['Meridional Volume Transport ',num2str(intStart), '~', num2str(intEnd),'dBar']);
                print('-depsc2',['./Figures/',destFolder,'/MvtProf',num2str(intStart),'_',num2str(intEnd),dataYear,'.eps']); %lon
        end
        close(figHandle);


        % Plot Mass Transport
        figHandle = figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

        profMasstAggr = mask.*intMassFluxExactProf;
%        profMasstAggr = mask.*intMassFluxProf;%./ gsw_grav(profLatAggrSel, mean(intStartList));
        scatterm(profLatAggrSel,profLongAggrSel,10, profMasstAggr,'.');
        caxis(quantile(profMasstAggr, [0.01, 0.99]));

        cb = colorbar;
        cb.Label.String = 'Mass Transport [kg / s]';  
        cb.FontSize = 14;

        cLims = caxis;
        colormap(darkb2r(cLims(1),cLims(2)));

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        %set(gcf,'pos',[0 0 22.5 15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        switch targetVar
            case 'lat'
                title(['Zonal Mass Transport ',num2str(intStart), '~', num2str(intEnd),'dBar']);
                print('-depsc2',['./Figures/',destFolder,'/ZmtProf',num2str(intStart),'_',num2str(intEnd),dataYear,'.eps']); %lat
            case 'lon'
                title(['Meridional Mass Transport ',num2str(intStart), '~', num2str(intEnd),'dBar']);
                print('-depsc2',['./Figures/',destFolder,'/MmtProf',num2str(intStart),'_',num2str(intEnd),dataYear,'.eps']); %lon
        end
        close(figHandle);

        % Plot Heat Transport
        figHandle = figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

        profHtAggr = mask.* gsw_cp0 .* intHeatFluxExactProf.*1e-15;
        scatterm(profLatAggrSel,profLongAggrSel,10, profHtAggr,'.');
        caxis(quantile(profHtAggr, [0.01, 0.99]));

        cb = colorbar;
        cb.Label.String = 'Heat Transport [PW / m]';  
        cb.FontSize = 14;

        cLims = caxis;
        colormap(darkb2r(cLims(1),cLims(2)));

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        %set(gcf,'pos',[0 0 22.5 15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        switch targetVar
            case 'lat'
                title(['Zonal Heat Transport ',num2str(intStart), '~', num2str(intEnd),'dBar']);
                print('-depsc2',['./Figures/',destFolder,'/ZhtProf',num2str(intStart),'_',num2str(intEnd),dataYear,'.eps']); %lat
            case 'lon'
                title(['Meridional Heat Transport ',num2str(intStart), '~', num2str(intEnd),'dBar']);
                print('-depsc2',['./Figures/',destFolder,'/MhtProf',num2str(intStart),'_',num2str(intEnd),dataYear,'.eps']); %lon
        end
        close(figHandle);

        % Plot Heat/Vol Transport
        figHandle = figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

        scatterm(profLatAggrSel,profLongAggrSel,10,...
         profHtAggr ./ profVoltAggr ./ mean(profHtAggr ./ profVoltAggr, 'omitnan'),'.');
        caxis(quantile(profHtAggr ./ profVoltAggr ./ mean(profHtAggr ./ profVoltAggr, 'omitnan'), [0.01, 0.99]));

        cb = colorbar;
        cb.Label.String = 'Ratio';  
        cb.FontSize = 14;

        colormap(parula);

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        %set(gcf,'pos',[0 0 22.5 15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        switch targetVar
            case 'lat'
                title(['Zonal Heat/Vol Ratio ',num2str(intStart), '~', num2str(intEnd),'dBar']);
                print('-depsc2',['./Figures/',destFolder,'/ZhtVtRatio',num2str(intStart),'_',num2str(intEnd),dataYear,'.eps']); %lat
            case 'lon'
                title(['Meridional Heat/Vol Ratio ',num2str(intStart), '~', num2str(intEnd),'dBar']);
                print('-depsc2',['./Figures/',destFolder,'/MhtVtRatio',num2str(intStart),'_',num2str(intEnd),dataYear,'.eps']); %lon
        end
        close(figHandle);


       % Plot Heat/Mass Transport
        figHandle = figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

        scatterm(profLatAggrSel,profLongAggrSel,10,...
         profHtAggr ./ profMasstAggr ./ mean(profHtAggr ./ profMasstAggr, 'omitnan'),'.');
        caxis(quantile(profHtAggr ./ profMasstAggr ./ mean(profHtAggr ./ profMasstAggr, 'omitnan'), [0.01, 0.99]));

        cb = colorbar;
        cb.Label.String = 'Ratio';  
        cb.FontSize = 14;

        cLims = caxis;
        colormap(parula);

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        %set(gcf,'pos',[0 0 22.5 15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        switch targetVar
            case 'lat'
                title(['Zonal Heat/Mass Ratio ',num2str(intStart), '~', num2str(intEnd),'dBar']);
                print('-depsc2',['./Figures/',destFolder,'/ZhtMtRatio',num2str(intStart),'_',num2str(intEnd),dataYear,'.eps']); %lat
            case 'lon'
                title(['Meridional Heat/Mass Ratio ',num2str(intStart), '~', num2str(intEnd),'dBar']);
                print('-depsc2',['./Figures/',destFolder,'/MhtMtRatio',num2str(intStart),'_',num2str(intEnd),dataYear,'.eps']); %lon
        end
        close(figHandle);


    end

end