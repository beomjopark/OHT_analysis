%% Wrapper function for createDataMask / filterUsingMasks
%% Run after selection and Merge

addpath(genpath('./cbrewer'));
addpath(genpath('../gsw_matlab'));
addpath(genpath('./Util'));

%% Parameters
%{
dataYear = '_2007_2018' % Use Merged Data
windowSize = 5
minNumberOfObs = 20
typeTag = 'target' %'int' % target;
%}


%% Run Selection
nCore = feature('numcores')
%{
refPres = 900
intStartList = [5, 10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800] %
 Full upper ocean
intStartList = [1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1850, 1900] % Mid Ocean
%}

isMid = (min(intStartList) > refPres)

%% QC by distance: BEFORE EVERYTHING ELSE
if ~is2step
    if isMid
        filterUsingDist(typeTag, responseTag, [1000, 1900]);
    else
        filterUsingDist(typeTag, responseTag, [5, 800]);
    end
end

%% Data Mask and Filter
if is2step && strcmp(typeTag, 'int') % This case is for intlatflux/intlonflux
    if isMid
        intStartList = [refPres, intStartList]
    else
%        intStartList = [10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800] % Full upper ocean from 10 with reflv
        if intStartList(1) == 10 && intStartList(end) == 800
            intStartList = [intStartList, refPres]
        end
    end
    fprintf('Target pressure: %d to %d\n', min(intStartList), max(intStartList));
    typeTag = strcat(typeTag, targetVar); %'intlat / intlon'
    createDataMask(typeTag, responseTag, intStartList, dataYear, windowType, windowSize, minNumberOfObs, isAdjusted, isAbsolute);
%    filterUsingMasks(typeTag, responseTag, intStartList, dataYear, windowSize, minNumberOfObs, false, fluxType, isAdjusted, isAbsolute);
    filterUsingMasks_Distrib(typeTag, responseTag, intStartList, dataYear, windowType, windowSize, minNumberOfObs, false, fluxType, isAdjusted, isAbsolute)    
else
    switch typeTag
        case 'target'
            fprintf('Target pressure: %d to %d\n', min(intStartList), max(intStartList));
            createDataMask(typeTag, responseTag, intStartList, dataYear, windowSize, minNumberOfObs, false, false);
            filterUsingMasks(typeTag, responseTag, intStartList, dataYear, windowSize, minNumberOfObs, false, [], isAdjusted, isAbsolute);        

        case 'int'
            %{
            poolobj = parpool(nCore, 'IdleTimeout', 1200);
            parfor intStartIdx = 1:numel(intStartList)
                intStart = intStartList(intStartIdx); % Parfor requires increasing 1 index
                switch typeTag
                    case 'int'
                        verticalSelection = strcat('Relative', num2str(intStart)); %
                    case 'target'
                        verticalSelection = num2str(intStart); % Deprecated
                end

                fprintf('Target pressure: %d\n', intStart);
                createDataMask(typeTag, responseTag, verticalSelection, dataYear, windowSize, minNumberOfObs);
                filterUsingMasks(typeTag, responseTag, verticalSelection, dataYear, windowSize, minNumberOfObs, false);
            end
            delete(poolobj)
            %}
           createDataMask_Distrib(typeTag, responseTag, intStartList, dataYear, windowType, windowSize, minNumberOfObs, isAdjusted, isAbsolute);
           filterUsingMasks_Distrib(typeTag, responseTag, intStartList, dataYear, windowType, windowSize, minNumberOfObs, false, [], isAdjusted, isAbsolute);
%{
            % QC2
            filterUsingThres_Distrib(typeTag, responseTag, intStartList, dataYear, windowSize, minNumberOfObs, false, [], isAdjusted, isAbsolute, kernelType, 10);
%}

        otherwise % lat lon
%            poolobj = parpool(nCore, 'IdleTimeout', 1200);
            for intStartIdx = 1:numel(intStartList)
                intStart = intStartList(intStartIdx); % Parfor requires increasing 1 index
                verticalSelection = strcat('Relative', num2str(intStart)); 

                fprintf('Target pressure: %d\n', intStart);
                createDataMask(typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, isAdjusted, isAbsolute);
                filterUsingMasks_Distrib(typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, false, fluxType, isAdjusted, isAbsolute, nAdjust);
            end
 %           delete(poolobj)
    end        
end
