%% Wrapper function for computeTotalOHT_Profile / integrateOHT_Profile
%% Run after selection and Merge

%$ Load GSW toolbox
addpath(genpath('../gsw_matlab'));
addpath(genpath('../cbrewer'));
addpath(genpath('./Util'));

%% Parameters
%Params_LonFlux_Step1
%{
dataYear = '_2007_2018' % Use Me\rged Data
windowSize = 5
minNumberOfObs = 20

responseTag = 'Dens' %'Sal'; 'Temp'; %'Flux'; if is2step

% Field Parameters
meanTag = 'NoTrend';
kernelType = 'Matern' %'Matern'; %'ExpGeom'
targetVar = 'lat' % 'lon'
%}
isAdjusted = true
isAbsolute = true
iterEM = 3

isPlot = false
%% Run Selection
%poolobj = parpool(36, 'IdleTimeout', 1200);
refPres = 900;
intStartList = [10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800] % Full upper ocean - 5 % 50:50:500;
intStartList = [intStartList, refPres]
for intStartIdx = 1:numel(intStartList)
    intStart = intStartList(intStartIdx); % Parfor requires increasing 1 index
    fprintf('Target pressure: %d\n', intStart);
    verticalSelection = strcat('Relative', num2str(intStart)); %'MidMeso';%'MidMeso';%'UpperOcean';%'Mikael';

    if intStart == refPres       
        computeTotalOHT_RefProfile(responseTag, verticalSelection, intStartList, intStart, dataYear, windowSize, minNumberOfObs, meanTag, kernelType, targetVar, isPlot, isAdjusted, isAbsolute)
        close all;
        continue;
    end
    computeTotalOHT_Profile(responseTag, verticalSelection, intStartList, intStart, dataYear, windowType, windowSize, minNumberOfObs, meanTag, kernelType, month, targetVar, isPlot, isAdjusted, isAbsolute, nAdjust, iterEM)
%    computeTotalOHT_ProfileDUACS(responseTag, verticalSelection, intStartList, intStart, dataYear, windowType, windowSize, minNumberOfObs, meanTag, kernelType, month, targetVar, isPlot, isAdjusted, isAbsolute, nAdjust, iterEM)
    close all;
end

%delete(poolobj);
%poolobj = parpool(36, 'IdleTimeout', 1200);
integrateOHT_Profile(intStartList, dataYear, minNumberOfObs, targetVar, isPlot, isAdjusted, isAbsolute);
%delete(poolobj);

%% For ENSO
integrateOHT_Profile([10, 15, 20, 30, 50, 75, 100], dataYear, minNumberOfObs, targetVar, isPlot, isAdjusted, isAbsolute);
integrateOHT_Profile([100, 125, 150, 200, 250, 300], dataYear, minNumberOfObs, targetVar, isPlot, isAdjusted, isAbsolute);
integrateOHT_Profile([300, 400, 500, 600, 700, 800, 900], dataYear, minNumberOfObs, targetVar, isPlot, isAdjusted, isAbsolute);
