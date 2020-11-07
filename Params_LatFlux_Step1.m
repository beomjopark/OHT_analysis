%% Setup for Step 1
% Data Type
typeTag = 'int' % target; % 'lat' 'lon'
responseTag = 'Dens' %'Sal'; 'Temp'; %'Flux'; if is2step
dataYear = '_2007_2018' % Use Merged Data

% Window Specification
windowType = 'spherical' % 'box'
windowSize = [4,4,4] % Mean/Cov/Krig
minNumberOfObs = 20

% MeanField
meanTag = 'NoTrend' % 'NoTrendVelSeas3'

% AnomalyField
kernelType = 'Matern' %'Matern'; %'ExpGeom'
month = 11

% If gridding OHT, set flag and responseTag, typeTag
is2step = false %true;

% Kriging: Prediction on derivative or original
isDeriv = false
targetVar = 'lat' % 'lon'
isProfile = false;

% 2 step specification
isStandardize = true
fluxType = []
eqBorder = []

% Bias Adjustment
isAdjusted = false%true
nAdjust = 0 % How many adjustment has been made so far?

% Extra Option
isAbsolute = false%true from totatl OHT % false for firststep adjustment
isFullMonth = false
