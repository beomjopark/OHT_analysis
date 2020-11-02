addpath(genpath('./Util'));

yearList = 2007:2018
varType = 'lat' % 'Dynh'; 'lon'
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

% Read DUACS
%%yearList = 2007:2018; % Now from the -v year=2007:2018

nYear = numel(yearList);
monthList = 1:12;


% Target ARGO
latListARGO = linspace(-89.5,89.5,180);
longListARGO = linspace(20.5,379.5,360);
[latGridARGO, longGridARGO] = meshgrid(latListARGO, longListARGO);

% Collect the time stamp
tic;
timeESAList = [];
dayNamesESACell = {};
for iYear = yearList
    for iMonth = 1:12
        targetYMStr = [num2str(iYear),'/',num2str(iMonth, '%02d')]

        curMonDir = dir([EsaFolder, '/', targetYMStr]);
        dayNames = {curMonDir.name};
        dayNames(1:2) = []; % Remove ., ..
        isNC = cellfun(@(name) strcmp(name(end-1:end), 'nc'), dayNames);
        dayNames = dayNames(isNC);

        for iDay = 1:length(dayNames)
            timeESAList = [timeESAList, ncread([EsaFolder, '/', targetYMStr, '/',cell2mat(dayNames(iDay))], 'time')];
            dayNamesESACell{end+1} = cell2mat(dayNames(iDay));
        end
    end
end
toc;

dateList = datenum(datetime(1981,1,1) + seconds(timeESAList)); % For ESA2017
parsedDateList = datevec(dateList);
parsedDateList = parsedDateList(:, 1:3);


tic
meanFlux = zeros(size(latGridARGO));
for cnt = 1:size(parsedDateList,1)
    iYear = num2str(parsedDateList(cnt,1));
    iMonth = num2str(parsedDateList(cnt,2), '%02d');

    FluxMon = load([outFolder,'/', iYear, '/', iMonth, '/',dayNamesESACell{cnt}(1:14),'.mat']);
    meanFlux = meanFlux + FluxMon.FluxInterp;
end
meanFlux = meanFlux ./ size(parsedDateList,1);
toc


save([outFolder,'mean',varType,'DUACS_ESA','.mat'],...
    'latGridARGO', 'longGridARGO', 'meanFlux')
