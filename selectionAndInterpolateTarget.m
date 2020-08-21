function selectionAndInterpolateTarget(dataYear, responseTag, intStart, isPlot)
%% Interpolate the target response at intStart from +-10 pressure grid.
%% If Multiple intStart is given, select the profile whose obs is in all range
%% and return interpolated Target for each intStart as a matrix
%%
%% Does not return object but save mat file
%% Note: Function version of selectionAndInterpolateTemp_target

    nIntStart = numel(intStart);
    if nIntStart > 1
        fprintf("Multiple intStart is provided.\n")
        isMultiple = true;
    else
        isMultiple = false;
    end
    tag = 'PchipPotTemp';

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

    startPresIdx = sum( cumprod((startPres > intStart'), 1), 1) + 1; %sum( cumprod((startPres <= intStart') == 0, 1), 1) + 1;
    %% Profile selection based on start(target) pressure
    if isMultiple
        selIdx = (startPres >=0 & startPres <= min(intStart) & endPres >= max(intStart)); 
    else
        selIdx = (startPres >=0 & startPres <= (intStart-10) & endPres >= (intStart+10)); 
    end


    profPresAggrSel = profPresAggrGood(selIdx);
    profTempAggrSel = profTempAggrGood(selIdx);
    profPsalAggrSel = profPsalAggrGood(selIdx);
    profLatAggrSel = profLatAggrGood(selIdx);
    profLongAggrSel = profLongAggrGood(selIdx);
    profYearAggrSel = profYearAggrGood(selIdx);
    profJulDayAggrSel = profJulDayAggrGood(selIdx);
    profFloatIDAggrSel = profFloatIDAggrGood(selIdx);

    %% Compute absolute salinity and conservative and potential temperature, you'll need to have the GSW toolbox in Matlab path to run this section, see http://www.teos-10.org/software.htm
    % Convert longitude from 20-380 range to 0-360 range
    profLongAggrSelTemp = (profLongAggrSel > 360).*(profLongAggrSel - 360) + (profLongAggrSel <= 360).*profLongAggrSel;

    % Compute absolute salinity, line below takes ~17 mins to run
    tic;
    %profAbsSalAggrSel = cellfun(@gsw_SA_from_SP,profPsalAggrSel,profPresAggrSel,num2cell(profLongAggrSelTemp),num2cell(profLatAggrSel),'UniformOutput',0);
    profAbsSalAggrSel = cell(1,numel(profPsalAggrSel));
    for iGrid = 1:numel(profAbsSalAggrSel)
        % Check sanity for error('gsw_SA_from_SP: pressure is out of range')
        badPresIdx = find(profPresAggrSel{iGrid} < -1.5 | profPresAggrSel{iGrid} > 12000);
        if badPresIdx
            profPsalAggrSel{iGrid}(badPresIdx) = [];
            profPresAggrSel{iGrid}(badPresIdx) = [];
            profTempAggrSel{iGrid}(badPresIdx) = [];
        end
        profAbsSalAggrSel{iGrid} = gsw_SA_from_SP(profPsalAggrSel{iGrid}, profPresAggrSel{iGrid}, profLongAggrSelTemp(iGrid), profLatAggrSel(iGrid));
    end
    toc;

    % Compute potential temperature (added by Anirban)
    switch responseTag
        case 'Temp'
            tic;
            profPotTempAggrSel = cellfun(@gsw_pt_from_t,profAbsSalAggrSel,profTempAggrSel,profPresAggrSel,'UniformOutput',0);
            toc;
        case 'Dens'
            % Compute conservative temperature
            tic;
            profConsTempAggrSel = cellfun(@gsw_CT_from_t,profAbsSalAggrSel,profTempAggrSel,profPresAggrSel,'UniformOutput',0);
            toc;

            % Compute in-situ density
            tic;
            profDensAggrSel = cellfun(@gsw_rho, profAbsSalAggrSel, profConsTempAggrSel, profPresAggrSel, 'UniformOutput', 0);
            toc;
    end

    %% Temperature interpolation at intStart with pchip
    nProfile = length(profPresAggrSel);
    tic;
    switch responseTag
        case 'Temp'
            targetTempProf = zeros(nIntStart, nProfile); % Potential Temperature
            parfor iGrid = 1:nProfile
            %for iGrid = 1:nProfile
                targetTempProf(:,iGrid)  = interp1(profPresAggrSel{iGrid}, profPotTempAggrSel{iGrid}, intStart, 'pchip');
            end

            if isMultiple
                save(['./Data/target',responseTag,'Prof',tag,num2str(min(intStart)),'_',num2str(max(intStart)),dataYear,'.mat'],...
                    'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel',...
                    'targetTempProf','startPresIdx','intStart');
            else
                save(['./Data/target',responseTag,'Prof',tag,num2str(intStart),dataYear,'.mat'],...
                    'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel',...
                    'targetTempProf','startPresIdx','intStart');
            end
        case 'Dens'
            targetDensProf = zeros(nIntStart, nProfile); % Potential Temperature
            %parfor iGrid = 1:nProfile
            parfor iGrid = 1:nProfile
                targetDensProf(:,iGrid)  = interp1(profPresAggrSel{iGrid}, profDensAggrSel{iGrid}, intStart, 'pchip');
            end

            if isMultiple
                save(['./Data/target',responseTag,'Prof',tag,num2str(min(intStart)),'_',num2str(max(intStart)),dataYear,'.mat'],...
                    'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel',...
                    'targetDensProf','startPresIdx','intStart');
            else
                save(['./Data/target',responseTag,'Prof',tag,num2str(intStart),dataYear,'.mat'],...
                    'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel',...
                    'targetDensProf','startPresIdx','intStart');
            end
            
        case 'Sal'
            targetSalProf = zeros(nIntStart, nProfile); % Potential Temperature
            %parfor iGrid = 1:nProfile
            for iGrid = 1:nProfile
                targetSalProf(:,iGrid)  = interp1(profPresAggrSel{iGrid}, profAbsSalAggrSel{iGrid}, intStart, 'pchip');
            end
            save(['./Data/target',responseTag,'Prof',tag,num2str(intStart),dataYear,'.mat'],...
                'profLatAggrSel','profLongAggrSel','profYearAggrSel','profJulDayAggrSel','profFloatIDAggrSel',...
                'targetSalProf','startPresIdx','intStart');
    end
    toc;


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

        switch responseTag
            case 'Temp'
                scatterm(profLatAggrSel(idx),profLongAggrSel(idx),10,targetTempProf(idx),'.');
            case 'Dens'
                scatterm(profLatAggrSel(idx),profLongAggrSel(idx),10,targetDensProf(idx),'.');
            case 'Sal'
                scatterm(profLatAggrSel(idx),profLongAggrSel(idx),10,targetSalProf(idx),'.');
        end

        colormap(jet(100));
        colorbar;
        %caxis([-1e3,3e4]);

        set(gcf,'units','centimeters')
        set(gcf,'pos',[0 0 1.5*22.5 1.5*15])
        set(gcf,'paperunits',get(gcf,'units')) 
        set(gcf,'paperpos',get(gcf,'pos'))
        print('-depsc2',['./Figures/target',num2str(intStart),responseTag,'Prof.eps']);
    end

end