%% Wrapper function for SelectionAndIntegrate

% Data Spliting for Walltime Limit
%{
if strcmp(dataYear, '')
    nPart = 4
    iPart = 4
else
    nPart = 2
    iPart = 2
end
%}


%$ Load GSW toolbox
addpath(genpath('./Util'));
addpath(genpath('../gsw_matlab'));


%% Extract the full dataset
if ~exist(['Data/Argo_data_aggr',dataYear,'.mat'], 'file')
    error('Data Does not exists!')
end
% cd ./Data/
% system('tar -xvzf Argo_data_aggr.tar.gz');
% cd ..


%% Run Selection : Inner parfor Version
nCore = feature('numcores')
%{
intStartList = [5, 10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800] % Full upper ocean
intStartList = [1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1850, 1900] % Mid Ocean
intStartList = [1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1850]%, 1900] % Mid Ocean
%}

isMid = (min(intStartList) > 900)

if strcmp(typeTag, 'target')
    if isMid
        selectionAndInterpolateTarget(dataYear, responseTag, [900, intStartList], false);
    else
        selectionAndInterpolateTarget(dataYear, responseTag, [intStartList, 900], false);
    end
else
    %{
    for intStart = intStartList
        poolobj = parpool(36, 'IdleTimeout', 1200);

        fprintf('Target pressure: %d\n', intStart);
%        if strcmp(typeTag, 'target')
%            selectionAndInterpolateTarget(dataYear, 'Temp', intStart, false);
%        else
        if exist('iPart', 'var')
            selectionAndIntegrate(dataYear, verticalSelection, intStart, false, iPart, nPart);
        else
            selectionAndIntegrate(dataYear, verticalSelection, intStart, false);
        end
%        end
        delete(poolobj)
    end
    %}
    % multiple in once in CT-SA interpolation: May incur bias!
    poolobj = parpool(nCore, 'IdleTimeout', 1200);
    if exist('iPart', 'var')
        if isMid
            selectionAndIntegrate_Mid(dataYear, verticalSelection, intStartList, false, iPart, nPart);
        else
            selectionAndIntegrate(dataYear, verticalSelection, intStartList, false, iPart, nPart);
        end
    else
        if isMid
            selectionAndIntegrate_Mid(dataYear, verticalSelection, intStartList, false);
        else
            selectionAndIntegrate(dataYear, verticalSelection, intStartList, false);
        end
    end
    delete(poolobj)
end



%% Outer parfor Version
%% : Sal calculation is significantly long
%% Remark) Requires Significant memory
%{
poolobj = parpool(36, 'IdleTimeout', 1200);
intStartList = 200:50:500;
parfor intStartIdx = 1:numel(intStartList)
    intStart = intStartList(intStartIdx); % Parfor requires increasing 1 index
    fprintf('Target pressure: %d\n', intStart);

    if strcmp(typeTag, 'target')
        selectionAndInterpolateTarget(dataYear, responseTag, intStart, false);
    else
        selectionAndIntegrate(dataYear, verticalSelection, intStart, false);
    end

end
delete(poolobj)
%}
