%% Fetch closest grid point of DUACS matching the insitu Argo

varType = 'latVel' % 'Dynh'; 'lonVel'
%addpath(genpath('../OHC_dynamics'));


switch varType
  case 'Dynh'
    typeName = 'adt';
  case 'latVel' % zonal 
    typeName = 'ugos';
  case 'lonVel' % meridional
    typeName = 'vgos';
end


%srcFolder = 'D:/dataset-duacs-rep-global-merged-allsat-phy-l4/' % LOCAL
srcFolder = './Misc/CMEMS/dataset-duacs-rep-global-merged-allsat-phy-l4/' % SERVER

% Read DUACS
yearList = 2007:2018;   nYear = numel(yearList);
monthList = 1:12;


% Collect the time stamp
timeList = [];
for iYear = yearList
    for iMonth = 1:12
        targetDateStr = [num2str(iYear),num2str(iMonth, '%02d'),'15']

        curMonDir = dir([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d')]);

        dayNames = {curMonDir.name};
        dayNames(1:2) = []; % Remove ., ..

        isNC = cellfun(@(name) strcmp(name(end-1:end), 'nc'), dayNames);
        dayNames = dayNames(isNC);

        for iDay = 1:length(dayNames)
            timeList = [timeList, ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'time')];
        end
    end
end

dateList = datenum(datetime(1950,1,1) + days(timeList));

% Lat/Lon Grid
latListDUACS = ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'latitude');
longListDUACS = ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'longitude');

[latGridDUACS, longGridDUACS] = meshgrid(latListDUACS, longListDUACS);
spatialGridDUACS = [latGridDUACS(:), longGridDUACS(:)];

% Read Argo data to matching timestamp
Params_LatFlux_Step1

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

ARGO = load(['./Data/',typeTag,'TempDens','Prof','PchipPotTemp',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);


%% For each day, check the closest grid
%% at most resolution: 1Day, 1/4 * 1/4 resolution
minDateGridIdx = NaN(size(ARGO.profJulDayAggrSel));
for idx = 1:numel(minDateGridIdx)
  [~, minDateGridIdx(idx)] = min(abs(ARGO.profJulDayAggrSel(idx) - dateList));
end
profJulDayAggrSel_Interp = dateList(minDateGridIdx);

cnt = 1
NNgridIdx = NaN(size(ARGO.profLatAggrSel));
distGap = NaN(size(ARGO.profLatAggrSel));



% Convert Lon to 0~360 range
profLonConvSel = ARGO.profLongAggrSel;
profLonConvSel(profLonConvSel > 360) = profLonConvSel(profLonConvSel > 360) - 360;
ARGO.profSpatialAggrSel = [ARGO.profLatAggrSel; profLonConvSel];


targetProf = NaN(size(ARGO.targetDynhProf));
for iYear = yearList
    disp(iYear);
    for iMonth = 1:12
        curMonDir = dir([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d')]);

        dayNames = {curMonDir.name};
        dayNames(1:2) = []; % Remove ., ..
        for iDay = 1:length(dayNames)
            isTargetProf = (minDateGridIdx == cnt);
            %% Find closeset Spatial Grid
            [NNgridIdx(isTargetProf), distGap(isTargetProf)] = dsearchn(spatialGridDUACS, ARGO.profSpatialAggrSel(:,isTargetProf)');


            %% Fetch ADT
            res = ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], typeName);

            targetProf(isTargetProf) = res(NNgridIdx(isTargetProf));
            cnt = cnt + 1;
        end
    end
end

switch varType
  case 'ADT'
    targetADTProf = targetProf;
  case 'latVel'
    targetlatVelProf = targetProf;
  case 'lonVel'
    targetlonVelProf = targetProf;
end

profLatAggrSel = ARGO.profLatAggrSel;
profLongAggrSel = ARGO.profLongAggrSel;
profJulDayAggrSel = ARGO.profJulDayAggrSel;
profspatialAggrSel_Interp = spatialGridDUACS(NNgridIdx, :);

saveName = ['./Data/int',varType,'DUACSProf','PchipPotTemp','Relative10',dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat'];

switch varType
  case 'ADT'
    save(saveName,...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
        'profJulDayAggrSel_Interp', 'profspatialAggrSel_Interp',...
        'targetADTProf');
  case 'latVel'
    save(saveName,...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
        'profJulDayAggrSel_Interp', 'profspatialAggrSel_Interp',...
        'targetlatVelProf');
  case 'lonVel'
    save(saveName,...
        'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
        'profJulDayAggrSel_Interp', 'profspatialAggrSel_Interp',...
        'targetlonVelProf');
end
