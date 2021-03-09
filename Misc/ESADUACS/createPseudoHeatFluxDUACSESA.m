%% Create latPseudoFlux / lonPseudoFlux

Params_LatFlux_Step1
targetVar = 'lon'

switch numel(windowSize)
  case 3
    windowSizeMean = windowSize(1);
    windowSizeCov = windowSize(2);
    windowSizeKrig = windowSize(3);
  case 2
    windowSizeMean = windowSize(1);
    windowSizeCov = windowSize(2);
    windowSizeKrig = windowSize(1);
   otherwise
    windowSizeMean = windowSize;
    windowSizeCov = windowSize;
    windowSizeKrig = windowSize;
end

switch windowType
    case 'spherical'
      windowTypeTag = 'spherical'
      windowSizeMargined = windowSizeMean * 2;
    otherwise
      windowTypeTag = []
      windowSizeMargined = windowSizeMean;
end

if isAdjusted
  adjustTag = 'Adjusted';
  adjustNumTag = ['Adjusted', num2str(nAdjust)];
  windowSizeTag = windowSizeFullTag;
else
  adjustTag = [];
  adjustNumTag = [];
end

if isAbsolute
    absoluteTag = 'Absolute';
else
    absoluteTag = [];
end


verticalSelection = 'Relative10'



DUACS = load(['./Data/int',targetVar,'Vel','DUACSProf','PchipPotTemp','Relative10',dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);

ESA = load(['./Data/intESAProf','PchipPotTemp','Relative10',dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);

% Check merger
if all(DUACS.profLatAggrSel == ESA.profLatAggrSel) & all(DUACS.profLongAggrSel == ESA.profLongAggrSel) & all(DUACS.profJulDayAggrSel == ESA.profJulDayAggrSel)
    disp('Profile Match')
else
    error('Profile Not matching!')
end

switch targetVar
    case 'lat'
        profFluxAggrSel = DUACS.targetlatVelProf .* ESA.targetSSTProf;
    case 'lon'
        profFluxAggrSel = DUACS.targetlonVelProf .* ESA.targetSSTProf;
    otherwise
        
end
% Filter NaN
nanCheck = isnan(profFluxAggrSel);
profFluxAggrSel = profFluxAggrSel(~nanCheck);

profLatAggrSel = DUACS.profLatAggrSel(~nanCheck);
profLongAggrSel = DUACS.profLongAggrSel(~nanCheck);
profJulDayAggrSel = DUACS.profJulDayAggrSel(~nanCheck);


DynhMaskName = ['./Data/dataMask',verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']


Params_LatFlux_Step2
typeTag = 'lon'
typeTag = ['DUACSESA', typeTag]
fluxType = 'heat'
tag = 'PchipPotTemp';
save(['./Data/',typeTag,fluxType,responseTag,'Prof',tag,verticalSelection,dataYear,[],[],'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat'],...
    'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
    'profFluxAggrSel');


% Duplicate dataMask
tgtMaskName = ['./Data/dataMask',typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']
copyfile(DynhMaskName, tgtMaskName, 'f')
