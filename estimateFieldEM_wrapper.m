%% Wrapper function for estimate MeanField and Prepare for localMLE data
%% Run after selection and Merge
addpath(genpath('../gsw_matlab'));
addpath(genpath('./cbrewer'));
addpath(genpath('./Util'));

%% Parameters
%{
meanTag = 'NoTrend' %'TrendDens';'TrendDens2';

dataYear = '_2007_2018' % Use Merged Data
windowSize = 5
minNumberOfObs = 20

typeTag = 'int' % target; % 'lat' 'lon'
responseTag = 'Dens' %'Sal'; 'Temp'; %'Flux'; if is2step

% If gridding OHT, set flag and typeTag
is2step = false %true;
%}
%{ 
Depricated: intlatlon
if is2step & ~(strcmp(typeTag, 'lat') | strcmp(typeTag, 'lon'))
    error('typeTag should be set to lat or lon\n');
end
%}

%h_YFParam = 3

%% Run Selection
nCore = feature('numcores')
refPres = 900
%{
intStartList = [10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800] % Full upper ocean
%intStartList = [1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1850, 1900] % Mid Ocean
nIterEM=3
%}

isMid = (min(intStartList) > refPres)

if is2step && strcmp(typeTag, 'int') % This case is for intlatflux/intlonflux
    if isMid
        intStartList = [refPres, intStartList]
    else
        intStartList = [10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800] % Full upper ocean from 10 with reflv
        intStartList = [intStartList, refPres]
    end
    fprintf('Target pressure: %d to %d\n', min(intStartList), max(intStartList));
    typeTag = strcat(typeTag, targetVar); %'intlat/intlon'

    estimateMeanField(meanTag, typeTag, responseTag, intStartList, dataYear, windowType, windowSize, minNumberOfObs, is2step, fluxType, eqBorder, isAdjusted, isAbsolute);

    verticalSelection = strcat(num2str(min(intStartList)),'_',num2str(max(intStartList)));
    subtractMeanSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, false, isStandardize, fluxType, eqBorder, isAdjusted, isAbsolute);
    divideDataToMonthsSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust);
    extendedDataSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, minNumberOfObs, is2step, isAdjusted, isAbsolute, nAdjust);
else
    %poolobj = parpool(nCore, 'IdleTimeout', 1200);
    for intStartIdx = 1:numel(intStartList)
        intStart = intStartList(intStartIdx); % Parfor requires increasing 1 index

        fprintf('Target pressure: %d\n', intStart);
        verticalSelection = strcat('Relative', num2str(intStart)); %'MidMeso';%'MidMeso';%'UpperOcean';%'Mikael';
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
