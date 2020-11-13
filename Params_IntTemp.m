%% Setup for Step 1: Anomaly lat Flux
% Data Type
typeTag = 'int' % target; % 'lat' 'lon'
responseTag = 'ESA' %'Temp' %'Sal'; 'Temp'; %'Flux'; if is2step
dataYear = '_2007_2018' % Use Merged Data

% Filter
windowType = 'spherical' % 'box'
windowSize = [4,4,4]
minNumberOfObs = 20

% MeanField
meanTag = 'NoTrend' % 'NoTrendVelSeas3'%'NoTrend'

% AnomalyField
kernelType = 'Matern' %'Matern'; %'ExpGeom'
month = 2
% If gridding OHT, set flag and responseTag, typeTag
is2step = false %true;
isProfile = false;

% Anomaly Prediction of velocity
isDeriv = false
targetVar = 'lat' % 'lon'

isStandardize = true
fluxType = []
eqBorder = []

% Adjustment
isAdjusted = false%true
isAbsolute = false%true from totatl OHT % false for firststep adjustment

nAdjust = 0 % How many adjustment has been made so far?
isFullMonth = false
