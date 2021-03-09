%% Fetch closest grid point of ESA SST matching the insitu Argo
%% Due to data size and computation limit, divide for each year:
%% Supply targetYear = 2007; ... 2018
%% isTarget = True to check for targetTemp

%addpath(genpath('../OHC_dynamics'));
%{
nCore = feature('numcores')
poolobj = parpool(nCore, 'IdleTimeout', 1200);
%}

%srcFolder = 'D:/SST_GLO_SST_L4_REP_OBSERVATIONS_010_024/' % LOCAL
srcFolder = './Misc/CMEMS/SST_GLO_SST_L4_REP_OBSERVATIONS_010_024/' % SERVER

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
        for iDay = 1:length(dayNames)
            timeList = [timeList, ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'time')];
        end
    end
end

dateList = datenum(datetime(1981,1,1) + seconds(timeList)); % For ESA2017

% Lat/Lon Grid
latListDUACS = ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'lat');
longListDUACS = ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'lon');

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

if isTarget
  ARGO = load(['./Data/targetTempProf','PchipPotTemp',num2str(10),'_',num2str(900),'_2007_2018','.mat'], 'profLongAggrSel', 'profLatAggrSel','profJulDayAggrSel');
else
  ARGO = load(['./Data/',typeTag,'TempDens','Prof','PchipPotTemp',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat']);
end

%% For each day, check the closest grid
%% at most resolution: 1Day, 1/4 * 1/4 resolution
minDateGridIdx = NaN(size(ARGO.profJulDayAggrSel));
for idx = 1:numel(minDateGridIdx)
  [~, minDateGridIdx(idx)] = min(abs(ARGO.profJulDayAggrSel(idx) - dateList));
end
profJulDayAggrSel_Interp = dateList(minDateGridIdx);

NNgridIdx = NaN(size(ARGO.profLatAggrSel));
distGap = NaN(size(ARGO.profLatAggrSel));



% Convert Lon to -180~180 range
profLonConvSel = ARGO.profLongAggrSel;
profLonConvSel(profLonConvSel > 180) = profLonConvSel(profLonConvSel > 180) - 360;
ARGO.profSpatialAggrSel = [ARGO.profLatAggrSel; profLonConvSel];

parsedDateList = datevec(dateList);
parsedDateList = parsedDateList(:, 1:3);
profSpatialAggrSel = ARGO.profSpatialAggrSel;

%{
NNgridIdx = cell(1, max(minDateGridIdx));
distGap = cell(1, max(minDateGridIdx));
%% Find closeset Spatial Grid

tic;
for cnt = find(parsedDateList(:,1) == targetYear)'%1:max(minDateGridIdx)
  isTargetProf = (minDateGridIdx == cnt);
  [NNgridIdx_cur, distGap_cur] = knnsearch(spatialGridDUACS, profSpatialAggrSel(:,isTargetProf)');
  NNgridIdx{cnt} = NNgridIdx_cur;
  distGap{cnt} = distGap_cur;
end
toc;
%}

%% Fetch ADT
disp('Fetch ADT')
cnt = 1;
targetSSTProf = NaN(size(ARGO.profJulDayAggrSel));
profspatialAggrSel_Interp = NaN(length(ARGO.profJulDayAggrSel), 2);
profspatialDistGap_Interp = NaN(size(ARGO.profJulDayAggrSel));

for iYear = yearList
    if isTarget
      saveName = ['./Data/targetTempProf','ESA','PchipPotTemp',num2str(10),'_',num2str(900),'_',num2str(iYear),'.mat'];
    else
      saveName = ['./Data/intESAProf','PchipPotTemp','Relative10',dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'Year',num2str(iYear),'.mat'];
    end
    NNpart = load(saveName,'NNgridIdx', 'distGap');

    for iMonth = 1:12
        disp([iYear, '/', iMonth]);
        curMonDir = dir([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d')]);

        dayNames = {curMonDir.name};
        dayNames(1:2) = []; % Remove ., ..
        for iDay = 1:length(dayNames)
            isTargetProf = (minDateGridIdx == cnt);

            %% Fetch ADT
            SST = ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'analysed_sst');
            targetSSTProf(isTargetProf) = SST(NNpart.NNgridIdx{cnt});
            profspatialAggrSel_Interp(isTargetProf,:) = spatialGridDUACS(NNpart.NNgridIdx{cnt}, :);
            profspatialDistGap_Interp(isTargetProf) = NNpart.distGap{cnt};
            cnt = cnt + 1;
        end
    end
end

profLatAggrSel = ARGO.profLatAggrSel;
profLongAggrSel = ARGO.profLongAggrSel;
profJulDayAggrSel = ARGO.profJulDayAggrSel;


if isTarget
  saveName = ['./Data/targetTempProf','ESA','PchipPotTemp',num2str(10),'_',num2str(900),'_2007_2018','.mat'];
else
  saveName = ['./Data/intESAProf','PchipPotTemp','Relative10',dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',num2str(windowSizeMean),'.mat'];
end
save(saveName,...
    'profLatAggrSel','profLongAggrSel','profJulDayAggrSel',...
    'profJulDayAggrSel_Interp', 'profspatialAggrSel_Interp', 'profspatialDistGap_Interp',...
    'targetSSTProf');
