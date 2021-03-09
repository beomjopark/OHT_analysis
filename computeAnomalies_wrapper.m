%% Wrapper function for localMLE and Anomaly Prediction
%% Run after selection / Merge and MeanField
addpath(genpath('../gsw_matlab'));
addpath(genpath('./Util'));
addpath(genpath('../OHC_dynamics'));

%% Parameters
%{
kernelType = 'Matern' %'Matern'; %'ExpGeom'
month = 2

dataYear = '_2007_2018' % Use Merged Data
windowSize = 5
minNumberOfObs = 20

typeTag = 'int' % target; % 'lat' 'lon'
responseTag = 'Dens' %'Sal'; 'Temp'; %'Flux'; if is2step

% If gridding OHT, set flag and responseTag, typeTag
is2step = false %true;

% Derivative Kernel: For Anomaly prediction of velocity
isDeriv = true
targetVar = 'lat' % 'lon'
%}


%% Sanity Check: Deprecated
%{
if is2step & ~(strcmp(typeTag, 'lat') | strcmp(typeTag, 'lon') | strcmp(responseTag, 'Flux'))
    error('typeTag should be set to lat or lon!\n');
end
%}
%iterEM = 0 : control from PBS

%% Run Selection
%% CHECK isProfile
nCore = feature('numcores')
refPres = 900

%intStartList = [10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800] % Full upper ocean
%{
intStartList = [1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1850, 1900] % Mid Ocean
intStartList = [10]
%}

isMid = (min(intStartList) > 900)

if is2step && strcmp(typeTag, 'int') % This case is for intlatflux/intlonflux
    if isMid
        intStartList = [refPres, intStartList]
    else
        intStartList = [10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800] % Full upper ocean from 10 with reflv
        intStartList = [intStartList, refPres]
    end
    fprintf('Target pressure: %d to %d\n', min(intStartList), max(intStartList));
    typeTag = strcat(typeTag, targetVar); %'intlatlon'
    verticalSelection = strcat(num2str(min(intStartList)),'_',num2str(max(intStartList)));

    poolobj = parpool(nCore-1, 'IdleTimeout', 1200);
    computeAnomaliesSeasonSpaceTime(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust);
    computeMeanAnomalies(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust);
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


        if is2step  % This case is for each latflux / lonflux
            poolobj = parpool(20, 'IdleTimeout', 1200);
            computeAnomaliesSeasonSpaceTime(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM);
            computeMeanAnomalies(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM);
        else
            % CHECK TO INCLUDE REFPRESS
            fprintf('Anomaly Computation')
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