addpath(genpath('./Util'));


%varType = 'lat' % 'Dynh'; 'lon'
%addpath(genpath('../OHC_dynamics'));

switch varType
  case 'Dynh'
    typeName = 'adt';
  case 'lat' % zonal 
    typeName = 'ugos';
  case 'lon' % meridional
    typeName = 'vgos';
end

%% Fetch DUACS
DuacsFolder = './Misc/CMEMS/dataset-duacs-rep-global-merged-allsat-phy-l4/' % SERVER
EsaFolder = './Misc/CMEMS/SST_GLO_SST_L4_REP_OBSERVATIONS_010_024/' % SERVER

outFolder = ['./Misc/CMEMS/DUACS_ESA_',varType,'Flux','/'] % SERVER
mkdir(outFolder)

% Read DUACS
%%yearList = 2007:2018; % Now from the -v year=2007:2018

nYear = numel(yearList);
monthList = 1:12;


% Collect the time stamp
iYear = yearList(1);
iMonth = 1;
iDay = 1;
curMonDir = dir([DuacsFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d')]);
dayNames = {curMonDir.name};
dayNames(1:2) = []; % Remove ., ..

latListDUACS = double(ncread([DuacsFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'latitude'));
longListDUACS = double(ncread([DuacsFolder, '/', num2str(iYear), '/', num2str(iMonth, '%02d'), '/',cell2mat(dayNames(iDay))], 'longitude'));

% Target ARGO
latListARGO = linspace(-89.5,89.5,180);
longListARGO = linspace(20.5,379.5,360);
[latGridARGO, longGridARGO] = meshgrid(latListARGO, longListARGO);

% Match ARGO long
longListDUACS(longListDUACS < 20) = longListDUACS(longListDUACS < 20) + 360;
[lonSortDUACS, lonSortDUACSIdx] = sort(longListDUACS);
[latGridDUACS, longGridDUACS] = meshgrid(double(latListDUACS), double(lonSortDUACS));
spatialGridDUACS = [latGridDUACS(:), longGridDUACS(:)];

% Collect the time stamp
tic;
timeESAList = [];
dayNamesESACell = {};
timeDUACSList = [];
dayNamesDUACSCell = {};
for iYear = yearList
    for iMonth = 1:12
        targetYMStr = [num2str(iYear),'/',num2str(iMonth, '%02d')]

        mkdir([outFolder, targetYMStr])

        curMonESADir = dir([EsaFolder, '/', targetYMStr]);
        dayNamesESA = {curMonESADir.name};
        dayNamesESA(1:2) = []; % Remove ., ..
        isNC = cellfun(@(name) strcmp(name(end-1:end), 'nc'), dayNamesESA);
        dayNamesESA = dayNamesESA(isNC);


        curMonDUACSDir = dir([DuacsFolder, '/', targetYMStr]);
        dayNamesDUACS = {curMonDUACSDir.name};
        dayNamesDUACS(1:2) = []; % Remove ., ..
        isNC = cellfun(@(name) strcmp(name(end-1:end), 'nc'), dayNamesDUACS);
        dayNamesDUACS = dayNamesDUACS(isNC);


        for iDay = 1:length(dayNamesESA)
            timeESAList = [timeESAList, ncread([EsaFolder, '/', targetYMStr, '/',cell2mat(dayNamesESA(iDay))], 'time')];
            dayNamesESACell{end+1} = cell2mat(dayNamesESA(iDay));

            timeDUACSList = [timeDUACSList, ncread([DuacsFolder, '/', targetYMStr, '/',cell2mat(dayNamesDUACS(iDay))], 'time')];
            dayNamesDUACSCell{end+1} = cell2mat(dayNamesDUACS(iDay));
        end
    end
end
toc;

dateList = datenum(datetime(1981,1,1) + seconds(timeESAList)); % For ESA2017
parsedDateList = datevec(dateList);
parsedDateList = parsedDateList(:, 1:3);


poolobj = parpool(10, 'IdleTimeout', 1200); % Constrained to 
tic
parfor cnt = 1:size(parsedDateList,1)
    iYear = num2str(parsedDateList(cnt,1));
    iMonth = num2str(parsedDateList(cnt,2), '%02d');

    ESA = load([EsaFolder,'/', iYear, '/', iMonth, '/',dayNamesESACell{cnt}(1:(end-3)),'_','interp','.mat']);

    DUACS = ncread([DuacsFolder, '/', iYear, '/', iMonth, '/',dayNamesDUACSCell{cnt}], typeName);

    Flux = DUACS .* ESA.SSTInterp;
    Flux = Flux(lonSortDUACSIdx, :);
    goodIdx = ~isnan(Flux);


    F = scatteredInterpolant(latGridDUACS(goodIdx), longGridDUACS(goodIdx), Flux(goodIdx), 'natural', 'none');
    FluxInterp = F(latGridARGO, longGridARGO);

    parsave_upscale([outFolder,'/', iYear, '/', iMonth, '/',dayNamesESACell{cnt}(1:14),'.mat'],...
        latGridARGO, longGridARGO, FluxInterp)
end
toc
