addpath(genpath('./Util'));

%% Fetch DUACS
srcFolder = './Misc/CMEMS/dataset-duacs-rep-global-merged-allsat-phy-l4/' % SERVER

% Read DUACS
%%yearList = 2007:2018; % Now from the -v year=2007:2018

nYear = numel(yearList);
monthList = 1:12;


% Collect the time stamp
iYear = yearList(1);
iMonth = 1;
iDay = 1;
curMonDir = dir([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d')]);
dayNames = {curMonDir.name};
dayNames(1:2) = []; % Remove ., ..

latListDUACS = ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'latitude');
longListDUACS = ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'longitude');

%adt = ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'adt');

[latGridDUACS, longGridDUACS] = meshgrid(double(latListDUACS), double(longListDUACS));
spatialGridDUACS = [latGridDUACS(:), longGridDUACS(:)];

%ncdisp([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))])

%% Fetch ESA and Regrid
%srcFolder = 'D:/SST_GLO_SST_L4_REP_OBSERVATIONS_010_024/' % LOCAL
srcFolder = './Misc/CMEMS/SST_GLO_SST_L4_REP_OBSERVATIONS_010_024/' % SERVER

curMonDir = dir([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d')]);
dayNames = {curMonDir.name};
dayNames(1:2) = []; % Remove ., ..

latListESA = ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'lat');
longListESA = ncread([srcFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'lon');

longListESA = longListESA + 180; % MATCH DUACS

[latGridESA, longGridESA] = meshgrid(double(latListESA), double(longListESA));
spatialGridESA = [latGridESA(:), longGridESA(:)];


% Collect the time stamp
tic;
timeList = [];
dayNamesCell = {};
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
            dayNamesCell{end+1} = cell2mat(dayNames(iDay));
        end
    end
end
toc;

dateList = datenum(datetime(1981,1,1) + seconds(timeList)); % For ESA2017
parsedDateList = datevec(dateList);
parsedDateList = parsedDateList(:, 1:3);

poolobj = parpool(6, 'IdleTimeout', 1200); % Constrained to 
tic
parfor cnt = 1:size(parsedDateList,1)
    SST = ncread([srcFolder, '/', num2str(parsedDateList(cnt,1)), '/', num2str(parsedDateList(cnt,2), '%02d'), '/',dayNamesCell{cnt}], 'analysed_sst');
    goodIdx = ~isnan(SST);

%    SSTInterp = griddata(latGridESA(goodIdx), longGridESA(goodIdx), SST(goodIdx),...
%        latGridDUACS, longGridDUACS, 'natural');
    F = scatteredInterpolant(latGridESA(goodIdx), longGridESA(goodIdx), SST(goodIdx), 'natural', 'none');
    SSTInterp = F(latGridDUACS, longGridDUACS);
    parsave_interp([srcFolder,'/', num2str(parsedDateList(cnt,1)), '/', num2str(parsedDateList(cnt,2), '%02d'), '/',dayNamesCell{cnt}(1:(end-3)),'_','interp','.mat'],...
        latGridDUACS, longGridDUACS, SSTInterp)
end
toc
