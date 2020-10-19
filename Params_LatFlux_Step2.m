%% Setup for Step 2
% Data Type
typeTag = 'lat'
responseTag = 'Flux'
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
is2step = true

% Kriging: Prediction on derivative or original
isDeriv = false
targetVar = typeTag % 'lon'

% 2 step specification
isStandardize = true
fluxType = 'heat' %'heat' 'mass' 'vol'
eqBorder = 2

% Bias Adjustment
isAdjusted = false%true
nAdjust = 0 % How many adjustment has been made so far?

% Extra Option
isAbsolute = true%true from totatl OHT % false for firststep adjustment
isFullMonth = false
