%% Wrapper function for estimate MeanField and Prepare for localMLE data
%% Run after selection and Merge
addpath(genpath('../gsw_matlab'));
addpath(genpath('./cbrewer'));
addpath(genpath('./Util'));

%% Run Selection
nCore = feature('numcores')
refPres = 900

isMid = (min(intStartList) > refPres)

if is2step && strcmp(typeTag, 'int') % This case is for intlatflux/intlonflux
    if isMid
        intStartList = [refPres, intStartList]
    else
        if intStartList(1) == 10 && intStartList(end) == 800
            intStartList = [intStartList, refPres]
        end
    end
    fprintf('Target pressure: %d to %d\n', min(intStartList), max(intStartList));
    typeTag = strcat(typeTag, targetVar); %'intlat/intlon'

    verticalSelection = strcat(num2str(min(intStartList)),'_',num2str(max(intStartList)));
    poolobj = parpool(nCore, 'IdleTimeout', 1200);
    for iterEM = 0:nIterEM
        estimateMeanField_EM(kernelType, month, meanTag, typeTag, responseTag, intStartList, dataYear, windowType, windowSize, minNumberOfObs, is2step, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM);

        subtractMeanSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, false, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth);
        divideDataToMonthsSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth)
        extendedDataSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, minNumberOfObs, is2step, isAdjusted, isAbsolute, nAdjust);
        close all;

        localMLESpaceTimeSeason(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType,windowSize, minNumberOfObs, is2step, false, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth);

        % Reoptimize parfor bugs
        localMLESpaceTimeSeason_Reestimate(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, false, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth);
    end
    delete(poolobj);
else
    %poolobj = parpool(nCore, 'IdleTimeout', 1200);
    for intStartIdx = 1:numel(intStartList)
        intStart = intStartList(intStartIdx); % Parfor requires increasing 1 index

        fprintf('Target pressure: %d\n', intStart);
        verticalSelection = strcat('Relative', num2str(intStart));
        poolobj = parpool(nCore, 'IdleTimeout', 1200);
        for iterEM = 0:nIterEM  %3
            if is2step
                estimateMeanField_EM(kernelType, month, meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM);
                subtractMeanSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, false, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth);
                divideDataToMonthsSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth)
            else
                estimateMeanField_EM(kernelType, month, meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, [], [], false, false, 0, iterEM);
    %           estimateMeanField_YFNN(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowSize, minNumberOfObs, is2step, [], [], h_YFParam);
               subtractMeanSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, false, isStandardize, [], [], isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth);
               divideDataToMonthsSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, [], [], isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth);
            end
            extendedDataSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, minNumberOfObs, is2step, isAdjusted, isAbsolute, nAdjust);
            close all;

            % MLE
            localMLESpaceTimeSeason(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType,windowSize, minNumberOfObs, is2step, false, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth);

            % Reoptimize parfor bugs
            localMLESpaceTimeSeason_Reestimate(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, false, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth);
        end
        delete(poolobj);
    end
end
