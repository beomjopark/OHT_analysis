function createDataMask(typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, isAdjusted, isAbsolute)
    %% Create DataMask which filters whose number of observation in each grid is larger than minNumberOfObs
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

    yearRange = strsplit(dataYear, '_');
    years = str2double(yearRange{2}) : str2double(yearRange{3});
    nYears = numel(years);

    %% Load Data
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
        srcName = ['./Data/',typeTag,responseTag,'Prof',tag,presString,dataYear,adjustTag,absoluteTag,'.mat']
        load(srcName);
    else
        switch responseTag
            case 'Flux'
                srcName = ['./Data/',typeTag,responseTag,'Prof',verticalSelection,dataYear,adjustTag,absoluteTag,'.mat']
                load(srcName);
                %% FILTER NAN count
                drop = isnan(profDerivSel);
                profLatAggrSel = profLatAggrSel(~drop);
                profLongAggrSel = profLongAggrSel(~drop);
                profJulDayAggrSel = profJulDayAggrSel(~drop);
            otherwise
                srcName = ['./Data/',typeTag,responseTag,'Prof',tag,verticalSelection,dataYear,adjustTag,absoluteTag,'.mat']
                load(srcName);
        end
    end

    
    %% Count number of data points for each month and year within each window
    [latGrid, longGrid] = meshgrid(linspace(-89.5,89.5,180), linspace(20.5,379.5,360));
    nGrid = numel(latGrid);

    % Enable wrap around by duplicating boundary data
    leftBoundaryIdx = find(profLongAggrSel <= 20 + windowSizeMean);
    rightBoundaryIdx = find(profLongAggrSel >= 380 - windowSizeMean);
    profLongAggrSel = [profLongAggrSel profLongAggrSel(leftBoundaryIdx) + 360 profLongAggrSel(rightBoundaryIdx) - 360];
    profLatAggrSel = [profLatAggrSel profLatAggrSel(leftBoundaryIdx) profLatAggrSel(rightBoundaryIdx)];
    profJulDayAggrSel = [profJulDayAggrSel profJulDayAggrSel(leftBoundaryIdx) profJulDayAggrSel(rightBoundaryIdx)];

    % Count observation for each Grid, Year, Month
    yearMonthCounts = zeros([size(latGrid), nYears, 12]);
    for iGrid = 1:nGrid
        if ~mod(iGrid, floor(nGrid/20))
            disp([int2str(iGrid), '/', int2str(nGrid)]);
        end
        
        latSel = latGrid(iGrid);
        longSel = longGrid(iGrid);

        latMin = latSel - windowSizeMargined;
        latMax = latSel + windowSizeMargined;
        longMin = longSel - windowSizeMargined;
        longMax = longSel + windowSizeMargined;
        
        idx = find(profLatAggrSel > latMin & profLatAggrSel < latMax & profLongAggrSel > longMin & profLongAggrSel < longMax);
        switch windowType
            case 'spherical'
              is_in_circle = distance(latSel, longSel, profLatAggrSel(idx), profLongAggrSel(idx),...
                                referenceEllipsoid('WGS84', 'm')) < refDist;
              idx = idx(is_in_circle);
            case 'box_var'
              idx = find(profLatAggrSel > latMin & profLatAggrSel < latMax & profLongAggrSel > longMin & profLongAggrSel < longMax);          
        end        
        [iGridSub1,iGridSub2] = ind2sub(size(latGrid), iGrid);
        nProf = length(idx);

        % Split Year-Month
        profJulDayAggrWindow = profJulDayAggrSel(idx)';
        profYearAggrWindow = zeros(1, nProf);
        profMonthAggrWindow = zeros(1, nProf);
        for iProf = 1:nProf
            temp = datevec(profJulDayAggrWindow(iProf));
            profYearAggrWindow(iProf) = temp(1);
            profMonthAggrWindow(iProf) = temp(2);
        end

        yearMonthCountsWindow = zeros(nYears, 12);
        for iYear = 1:nYears
            yearMonthCountsWindow(iYear, :) = histcounts(profMonthAggrWindow(profYearAggrWindow == years(iYear)), 0.5:12.5);
        end
        yearMonthCounts(iGridSub1, iGridSub2, :, :) = yearMonthCountsWindow;   
    end

    %% Form data-driven landmask
    dataMask = zeros(size(latGrid));
    for iGrid = 1:nGrid
        [iGridSub1, iGridSub2] = ind2sub(size(latGrid), iGrid);
        %{
        predLat = latGrid(iGrid);
        predLong = longGrid(iGrid);
        %}
        yearMonthCountsWindow = squeeze(yearMonthCounts(iGridSub1, iGridSub2, :, :));
        dataMask(iGrid) = all(sum(yearMonthCountsWindow,1) >= minNumberOfObs);        
    end
    dataMask(dataMask == 0) = NaN;

    switch typeTag
        case 'target'
            if isnumeric(verticalSelection)
                saveName = ['./Data/dataMask',typeTag,responseTag,presString,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']
            else
                saveName = ['./Data/dataMask',typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']
            end
        case 'int'
            if isnumeric(verticalSelection)
                saveName = ['./Data/dataMask',presString,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']
            else
                saveName = ['./Data/dataMask',verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']
            end
        otherwise %'lat' or 'lon'
            if isnumeric(verticalSelection) % intlatlon
                saveName = ['./Data/dataMask',typeTag,responseTag,presString,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']
            else
                saveName = ['./Data/dataMask',typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']
            end
    end
    save(saveName, 'dataMask', 'latGrid', 'longGrid');

end
