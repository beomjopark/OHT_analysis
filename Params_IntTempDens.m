%% Setup for Step 1: Selection Int 
% Data Type
typeTag = 'int' % target; % 'lat' 'lon'
switch typeTag
    case 'int'
        responseTag = 'TempDens' %'Dens';%'Sal';%'Dens';
    case 'target'
        responseTag = 'Temp'
    otherwise %'lat' 'lon'
        responseTag = 'Flux'    
end
%dataYear = '_2007_2018' % Use Merged Data

% Filter
windowType = 'spherical' % 'box'
windowSize = 4
minNumberOfObs = 20

% MeanField
meanTag = 'NoTrend';

% AnomalyField
kernelType = 'Matern' %'Matern'; %'ExpGeom'
month = 2
% If gridding OHT, set flag and responseTag, typeTag
is2step = false %true;

% Anomaly Prediction of velocity
isDeriv = false
targetVar = 'lat' % 'lon'

verticalSelection = 'Relative'     %'MidMeso';%'Mikael'; 'Anirban'; 'UpperOcean'; 'MidOcean';

isAdjusted = false
isAbsolute = false

nAdjust = 0 %1%0%2 % How many adjustment has been made so far?
