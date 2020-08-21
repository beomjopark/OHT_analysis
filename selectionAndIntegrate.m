function selectionAndIntegrate(dataYear, verticalSelection, intStart, isPlot, iPart, nPart)
%% Filter profile based on verticalSelection and Integrate the potential Temperature and Density within the pressure region.
%% Does not return object but save mat file
%% Note: Function version of selectionAndVerticalIntegrationPchipTrapzInterpolation_TempDens

    nIntStart = numel(intStart);
    if nIntStart > 1
        fprintf("Multiple intStart is provided.\n")
        isMultiple = true;
    else
        isMultiple = false;
    end
    tag = 'PchipPotTemp';

    % intStart Should be provided if Relative Selection before the code runs
    if strcmp(verticalSelection, 'Relative') && (nargin == 2)
        error('Relative Profile requires intStart input!');
    end

    %% Load data
    S = load(['./Data/Argo_data_aggr',dataYear,'.mat']);

    %% Filter bad observation
    profPresAggrGood = S.profPresAggr;
    profTempAggrGood = S.profTempAggr;
    profPsalAggrGood = S.profPsalAggr;
    profLatAggrGood = S.profLatAggr;
    profLongAggrGood = S.profLongAggr;
    profYearAggrGood = S.profYearAggr;
    profJulDayAggrGood = S.profJulDayAggr;
    profFloatIDAggrGood = S.profFloatIDAggr;
    
    clear S;

    % Data QC
    if strcmp(dataYear, '')
        badIdx = 875223;
        profPresAggrGood(badIdx) = [];
        profTempAggrGood(badIdx) = [];
        profPsalAggrGood(badIdx) = [];

        profLatAggrGood(badIdx) = [];
        profLongAggrGood(badIdx) = [];
        profYearAggrGood(badIdx) = [];
        profJulDayAggrGood(badIdx) = [];
        profFloatIDAggrGood(badIdx) = [];
    end

    % Remove duplicated Floats
    [~,ia,~] = unique([profLatAggrGood',profLongAggrGood',profJulDayAggrGood', profFloatIDAggrGood'], 'rows');
    fprintf('Removed %d duplicated Floats in the same space-time.\n', numel(profLatAggrGood) - numel(ia));
    profPresAggrGood = profPresAggrGood(ia);
    profTempAggrGood = profTempAggrGood(ia);
    profPsalAggrGood = profPsalAggrGood(ia);
    profLatAggrGood = profLatAggrGood(ia);
    profLongAggrGood = profLongAggrGood(ia);
    profYearAggrGood = profYearAggrGood(ia);
    profJulDayAggrGood = profJulDayAggrGood(ia);
    profFloatIDAggrGood = profFloatIDAggrGood(ia);
    

    % Remove duplicated (x,y)
    [~,ia,~] = unique([profLatAggrGood',profLongAggrGood',profJulDayAggrGood'], 'rows');
    dupIdx = find(not(ismember(1:numel(profLatAggrGood), ia)));
    removeIdx = [];
    for idx = dupIdx
        insp = find(profLatAggrGood == profLatAggrGood(idx)...
            & profLongAggrGood == profLongAggrGood(idx)...
            & profJulDayAggrGood == profJulDayAggrGood(idx));
        
        % Check duplicate profile
        if numel(profPresAggrGood{insp(1)}) == numel(profPresAggrGood{insp(2)})
            if ~sum(profTempAggrGood{insp(2)} - profTempAggrGood{insp(1)})
                removeIdx = [removeIdx, idx];
            end
        end    
    end
    fprintf('Removed Additional Duplicated %d Obs\n', numel(removeIdx));
    profPresAggrGood(removeIdx) = [];
    profTempAggrGood(removeIdx) = [];
    profPsalAggrGood(removeIdx) = [];
    profLatAggrGood(removeIdx) = [];
    profLongAggrGood(removeIdx) = [];
    profYearAggrGood(removeIdx) = [];
    profJulDayAggrGood(removeIdx) = [];
    profFloatIDAggrGood(removeIdx) = [];    
    
    % Still duplicated x yet unknown reason
    [~,ia,~] = unique([profLatAggrGood',profLongAggrGood',profJulDayAggrGood'], 'rows');
    fprintf('Removed %d duplicated same space-time obs (Need care).\n', numel(profLatAggrGood) - numel(ia));
    profPresAggrGood = profPresAggrGood(ia);
    profTempAggrGood = profTempAggrGood(ia);
    profPsalAggrGood = profPsalAggrGood(ia);
    profLatAggrGood = profLatAggrGood(ia);
    profLongAggrGood = profLongAggrGood(ia);
    profYearAggrGood = profYearAggrGood(ia);
    profJulDayAggrGood = profJulDayAggrGood(ia);
    profFloatIDAggrGood = profFloatIDAggrGood(ia);

    startPres = cellfun(@min,profPresAggrGood);
    endPres = cellfun(@max,profPresAggrGood);

    %% Profile selection based on start and end pressure
    %%% MidMeso and FullMeso are legacy of Relative selection
    switch verticalSelection
        case 'MidMeso'
            intStart = 600;
            intEnd = 975;        
            verticalSelection = 'Relative';
        case 'FullMeso'
            intStart = 200;
            intEnd = 975;
            verticalSelection = 'Relative';
        case 'Relative'
            % Relative profile
            if min(intStart) > 900 % Deep profile
                intEnd = 1850
            else
                intEnd = 900%975;
            end
    end

    switch verticalSelection
        case 'Anirban'
            selIdx = (startPres >=0 & startPres <= 10 & endPres >= 1800); % Anirban's selection, retains 59.64% of profiles
            intStart = 10;
            intEnd = 1800;
        case 'Mikael'
            selIdx = (startPres >=0 & startPres < 10 & endPres > 1800 & endPres <= 2100); % Mikael's selection, retains 59.14% of profiles
            intStart = 10;
            intEnd = 1800;
        case 'UpperOcean'
            selIdx = (startPres >=0 & startPres <= 15 & endPres >= 975); % Optimized upper ocean selection, retains 88.62% of profiles
            intStart = 15;
            intEnd = 975;
        case 'MidOcean'
            selIdx = (startPres <= 975 & endPres >= 1850); % Optimized mid-ocean selection, retains 62.98% of profiles
            intStart = 975;
            intEnd = 1850;
        case 'Relative'
            if isMultiple
                selIdx = (startPres >=0 & startPres <= max(intStart) & endPres >= intEnd); % Take the most lenient condition
            else
                selIdx = (startPres >=0 & startPres <= intStart & endPres >= intEnd);
            end
    end

    if nargin > 4 % iPart exists
        fprintf("Part %d/%d\n", iPart, nPart);

        nOriginal = length(startPres);
        nInterval = ceil(nOriginal / nPart);

        partIdx = false(size(startPres));
        if iPart == 1
            partIdx(iPart : (iPart * nInterval)) = true;
        elseif iPart == nPart
            partIdx((1+(iPart-1) * nInterval) : nOriginal) = true;
        else
            partIdx((1+(iPart-1)*nInterval) : iPart * nInterval) = true;
        end
        
        % Only consider the part
        selIdx = selIdx & partIdx;
    end

%{
                % Note that this breaks down the backward compatibility
                selIdx = struct;
                for intPres = intStart
                    selIdx.(intPres) = (startPres >=0 & startPres <= intPres & endPres >= intEnd);
                end
%}

    profPresAggrSel = profPresAggrGood(selIdx);
    profTempAggrSel = profTempAggrGood(selIdx);
    profPsalAggrSel = profPsalAggrGood(selIdx);
    profLatAggrSel = profLatAggrGood(selIdx);
    profLongAggrSel = profLongAggrGood(selIdx);
    profYearAggrSel = profYearAggrGood(selIdx);
    profJulDayAggrSel = profJulDayAggrGood(selIdx);
    profFloatIDAggrSel = profFloatIDAggrGood(selIdx);

    startPres = startPres(selIdx); % for additional selection later.

    %% Compute absolute salinity and conservative and potential temperature, you'll need to have the GSW toolbox in Matlab path to run this section, see http://www.teos-10.org/software.htm
    % Convert longitude from 20-380 range to 0-360 range
    profLongAggrSelTemp = (profLongAggrSel > 360).*(profLongAggrSel - 360) + (profLongAggrSel <= 360).*profLongAggrSel;

    % Compute absolute salinity, line below takes ~17 mins to run
    %profAbsSalAggrSel = cellfun(@gsw_SA_from_SP,profPsalAggrSel,profPresAggrSel,num2cell(profLongAggrSelTemp),num2cell(profLatAggrSel),'UniformOutput',0);
    nProfile = length(profPresAggrSel);
    profAbsSalAggrSel = cell(1,nProfile);
    profConsTempAggrSel = cell(1,nProfile);
    profPotTempAggrSel = cell(1,nProfile);
    profDensAggrSel = cell(1,nProfile);
    
    fprintf('Sanity Check for SA_from_SP\n');
    tic;
    for iGrid = 1:nProfile
        if(~mod(iGrid, nProfile/20))
            fprintf('%d / %d\n',iGrid, numel(profAbsSalAggrSel));
        end
        % Check sanity for error('gsw_SA_from_SP: pressure is out of range')
        badPresIdx = find(profPresAggrSel{iGrid} < -1.5 | profPresAggrSel{iGrid} > 12000);
        if badPresIdx
            profPsalAggrSel{iGrid}(badPresIdx) = [];
            profPresAggrSel{iGrid}(badPresIdx) = [];
            profTempAggrSel{iGrid}(badPresIdx) = [];
        end
    end
    toc;

    fprintf('Absolute Salinity, Conservative Temperature and in-situ density Computation\n');
    tic;
    parfor iGrid = 1:numel(profAbsSalAggrSel)
        % Compute absolute salinity, line below takes ~17 mins to run        
        profAbsSalAggrSel{iGrid} = gsw_SA_from_SP(profPsalAggrSel{iGrid}, profPresAggrSel{iGrid}, profLongAggrSelTemp(iGrid), profLatAggrSel(iGrid));
        % Compute conservative temperature
        profConsTempAggrSel{iGrid} = gsw_CT_from_t(profAbsSalAggrSel{iGrid}, profTempAggrSel{iGrid}, profPresAggrSel{iGrid});
        % Compute potential temperature (added by Anirban)
        profPotTempAggrSel{iGrid} = gsw_pt_from_t(profAbsSalAggrSel{iGrid}, profTempAggrSel{iGrid}, profPresAggrSel{iGrid});
        % Compute in-situ density
        profDensAggrSel{iGrid} = gsw_rho(profAbsSalAggrSel{iGrid}, profConsTempAggrSel{iGrid}, profPresAggrSel{iGrid});
    end
    toc;
    %{
    % Compute conservative temperature
    tic;
    profConsTempAggrSel = cellfun(@gsw_CT_from_t,profAbsSalAggrSel,profTempAggrSel,profPresAggrSel,'UniformOutput',0);
    toc;

    % Compute potential temperature (added by Anirban)
    tic;
    profPotTempAggrSel = cellfun(@gsw_pt_from_t,profAbsSalAggrSel,profTempAggrSel,profPresAggrSel,'UniformOutput',0);
    toc;

    % Compute in-situ density
    tic;
    profDensAggrSel = cellfun(@gsw_rho, profAbsSalAggrSel, profConsTempAggrSel, profPresAggrSel, 'UniformOutput', 0);
    toc;
    %}
%{    
%    Check infunnel to see inside oceanographic funnel: Approximation is accurate
    isInFunnel = zeros(1,numel(profAbsSalAggrSel));
    for iGrid = 1:numel(profAbsSalAggrSel)
        if(~mod(iGrid, numel(profAbsSalAggrSel)/20))
            disp([iGrid, '/',numel(profAbsSalAggrSel) ]);
        end
        isInFunnel(iGrid) = sum(gsw_infunnel(profAbsSalAggrSel{iGrid}, profConsTempAggrSel{iGrid}, profPresAggrSel{iGrid})) == numel(profAbsSalAggrSel{iGrid});
    end
%}

    fprintf('Dynamic Height Computation\n');
    targetDynhProf = NaN(nIntStart, nProfile);

    % Check the top pressure bottle Idx for each profiles
    startPresIdx = sum( cumprod((startPres > intStart'), 1), 1) + 1; %sum( cumprod((startPres <= intStart') == 0, 1), 1) + 1;

    tic;
    parfor iGrid = 1:nProfile
        startPresIgrid = startPresIdx(iGrid);
        intQuery = intStart(startPresIgrid:end); %feasible bottles

        % Check the overlapping intQuery with pressure bottles and include the non-existing query bottles as a interpolation points
        [commonInt, targetIdx, ~] = intersect(profPresAggrSel{iGrid}, intQuery);
        if numel(commonInt) == numel(intQuery) %all(profPresAggrSel{iGrid} == intQuery)
            profDynhAggrSel = gsw_geo_strf_dyn_height(profAbsSalAggrSel{iGrid}, profConsTempAggrSel{iGrid}, profPresAggrSel{iGrid}, intEnd);
        else
            % First Interpolate the target intQuery
           pres_i = [profPresAggrSel{iGrid}; setdiff(intQuery, commonInt)'];
           pres_i = sort(pres_i);
           [~, targetIdx, ~] = intersect(pres_i, intQuery);
%           targetIdx = (pres_i == intQuery);
           
           [SA_i, CT_i] =  gsw_SA_CT_interp(profAbsSalAggrSel{iGrid},profConsTempAggrSel{iGrid},profPresAggrSel{iGrid},pres_i);
            if any(isnan(SA_i))
                [Inan] = find(isnan(SA_i));
                [SA_i(Inan), CT_i(Inan)] = gsw_linear_interp_SA_CT(profAbsSalAggrSel{iGrid},profConsTempAggrSel{iGrid},profPresAggrSel{iGrid},pres_i(Inan));
            end 
           profDynhAggrSel = gsw_geo_strf_dyn_height(SA_i, CT_i, pres_i, intEnd);
        end
    %{
    % Check interpolation makes spurious term
    sum(SA_i([1:(targetIdx-1), (targetIdx+1):end]) -  profAbsSalAggrSel{iGrid})
    sum(CT_i([1:(targetIdx-1), (targetIdx+1):end]) -  profConsTempAggrSel{iGrid})

    figure;
    hold on;
    plot(profAbsSalAggrSel{iGrid},profConsTempAggrSel{iGrid}, '-ok')
    plot(SA_i, CT_i, '-+r')
    hold off;
    [SA_i(targetIdx), CT_i(targetIdx)]
    %}
        targetDynhProf(:, iGrid) = [NaN(startPresIgrid-1, 1); profDynhAggrSel(targetIdx)];
    end

%   profDynhAggrSel = cellfun(@gsw_geo_strf_dyn_height,...
%       profAbsSalAggrSel, profConsTempAggrSel, profPresAggrSel, repmat({intEnd}, size(profLatAggrSel)), 'UniformOutput', false);   
%   profDynhAggrSel{iGrid} = gsw_geo_strf_dyn_height(profAbsSalAggrSel{iGrid}, profConsTempAggrSel{iGrid}, profPresAggrSel{iGrid}, intEnd);
   toc;
    
    %% Vertical integration of the potential temperature and conservative temperature profiles
    fprintf('Vertical integration of Potential Temp and in-situ Dens\n');
    nGridx = 5000;
    gridx = linspace_vec(intStart', intEnd.*ones(size(intStart')), nGridx); % Fine grid, takes about 13 mins to run

    nProfile = length(profPresAggrSel);
    intTempProf = NaN(nIntStart, nProfile); % Potential Temperature
    intDensProf = NaN(nIntStart, nProfile);

    tic;
    parfor iGrid = 1:nProfile
        % Select the profile
        press=profPresAggrSel{iGrid};
        pottemp=profPotTempAggrSel{iGrid};
        dens = profDensAggrSel{iGrid};
                
        % Interpolate grids with pchip
        intTempIgrid = NaN(nIntStart, 1);
        intDensIgrid = NaN(nIntStart, 1);
        for idxIntStart = startPresIdx(iGrid):nIntStart
            iGridx = gridx(idxIntStart,:);
            pchipIntPot=interp1(press, pottemp, iGridx, 'pchip');
            pchipIntDens=interp1(press, dens, iGridx, 'pchip');
            
            % Integrate with trapizoidal rule
            intTempIgrid(idxIntStart)  = trapz(iGridx, pchipIntPot);
            intDensIgrid(idxIntStart)  = trapz(iGridx, pchipIntDens); 
        end
        intTempProf(:,iGrid) = intTempIgrid;
        intDensProf(:,iGrid) = intDensIgrid;
    end
    toc;


    if isMultiple
        intStartString = [num2str(min(intStart)),'_',num2str(max(intStart))];
    else
        intStartString = num2str(intStart);
    end
    if nargin > 4 % part
        if strcmp(verticalSelection, 'Relative')
            saveName = ['./Data/intTempDensProf',tag,verticalSelection,intStartString,dataYear,'Part',num2str(iPart),'_',num2str(nPart),'.mat'];
        else
            saveName = ['./Data/intTempDensProf',tag,verticalSelection,dataYear,'Part',num2str(iPart),'_',num2str(nPart),'.mat'];
        end
        save(saveName,...
            'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel', ...
            'intTempProf', 'intDensProf','targetDynhProf','intStart','intEnd', 'startPresIdx');
    else
        if strcmp(verticalSelection, 'Relative')
            saveName = ['./Data/intTempDensProf',tag,verticalSelection,intStartString,dataYear,'.mat'];
        else
            saveName = ['./Data/intTempDensProf',tag,verticalSelection,dataYear,'.mat'];
        end
        save(saveName,...
            'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel', ...
            'intTempProf', 'intDensProf','targetDynhProf','intStart','intEnd', 'startPresIdx');
    end

    %% Plot the profile
    if isPlot
        idx = 1:nProfile;

        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

%        scatterm(profLatAggrSel(idx),profLongAggrSel(idx),10,intTempProf(idx),'.');
        scatterm(profLatAggrSel,profLongAggrSel,10,targetDynhProf,'.');

%        colormap(jet(100));
        colormap(parula);
        colorbar;
        caxis(quantile(targetDynhProf, [0.01, 0.99]))
        %caxis([-1e3,3e4]);

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        print('-depsc2','./Figures/intTempProf.eps');


        figure;
        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');

        load coast;
        plotm(lat,long,'k');

        scatterm(profLatAggrSel(idx),profLongAggrSel(idx),10,intDensProf(idx),'.');

        colormap(jet(100));
        colorbar;
        %caxis([1.844e06,1.848e06]);

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        print('-depsc2','./Figures/intDensProf.eps');
    end

end