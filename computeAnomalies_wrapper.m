%% Wrapper function for localMLE and Anomaly Prediction
%% Run after selection / Merge and MeanField
addpath(genpath('../gsw_matlab'));
addpath(genpath('./Util'));
addpath(genpath('../OHC_dynamics'));

%% Run Selection
nCore = feature('numcores')
refPres = 900

isMid = (min(intStartList) > 900)

if is2step && strcmp(typeTag, 'int') % This case is for intlatflux/intlonflux
    if isMid
        intStartList = [refPres, intStartList]
    else
        if intStartList(1) == 10 && intStartList(end) == 800
            intStartList = [intStartList, refPres]
        end
    end
    fprintf('Target pressure: %d to %d\n', min(intStartList), max(intStartList));
    typeTag = strcat(typeTag, targetVar); %'intlatlon'
    verticalSelection = strcat(num2str(min(intStartList)),'_',num2str(max(intStartList)));

    poolobj = parpool(nCore-1, 'IdleTimeout', 1200);
    computeAnomaliesSeasonSpaceTime(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM);
    computeMeanAnomalies(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM);    
    delete(poolobj);
else
    for intStartIdx = 1:numel(intStartList)
        intStart = intStartList(intStartIdx); % Parfor requires increasing 1 index

        fprintf('Target pressure: %d\n', intStart);
        verticalSelection = strcat('Relative', num2str(intStart)); %'MidMeso';%'MidMeso';%'UpperOcean';%'Mikael';

        if isResetRes
            fprintf('Reset Residual')
            if is2step
                divideDataToMonthsSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth);
            else
               divideDataToMonthsSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, [], [], isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth);
            end
            extendedDataSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, minNumberOfObs, is2step, isAdjusted, isAbsolute, nAdjust);
            close all;
        end

        fprintf('Anomaly Computation')
        if is2step  % This case is for each latflux / lonflux
            poolobj = parpool(nCore, 'IdleTimeout', 1200);

            computeAnomaliesSeasonSpaceTime(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM);
            computeAnomaliesLinearCovariance(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM);            
            computeMeanAnomalies(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM);
        else
            % CHECK TO INCLUDE REFPRESS
            if isProfile
                poolobj = parpool(18, 'IdleTimeout', 1200);
                if isMid
                    computeAnomaliesSeasonSpaceTime_Profile(kernelType, month, typeTag, responseTag, verticalSelection, [900, intStartList], dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, isStandardize, [], [], isAdjusted, false, nAdjust, iterEM);
                else
                    computeAnomaliesSeasonSpaceTime_Profile(kernelType, month, typeTag, responseTag, verticalSelection, [10, 900], dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, isStandardize, [], [], isAdjusted, false, nAdjust, iterEM); % isAbsolute is ineffective here
    %                checkAnomaliesSeasonSpaceTime_Profile(kernelType, month, typeTag, responseTag, verticalSelection, [intStartList, 900], dataYear, windowSize, minNumberOfObs, isDeriv, targetVar, isStandardize, isAdjusted, isAbsolute, nAdjust, false);
                end
            else
                poolobj = parpool(24, 'IdleTimeout', 1200);
                computeAnomaliesSeasonSpaceTime(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, isStandardize, [], [], isAdjusted, isAbsolute, nAdjust, iterEM);
                computeMeanAnomalies(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, [], [], isAdjusted, isAbsolute, nAdjust, iterEM);
            end

        end
       delete(poolobj);
    end
end