%% Wrapper function for mergeSelection
%% Run after selection


%% Run Selection
%intStartList = [5, 10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800] % Full upper ocean;
%{
intStartList = [1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1850, 1900] % Mid Ocean
intStartList = [1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1850]%, 1900] % Mid Ocean
%}


isMid = (min(intStartList) > 900)

switch typeTag
    case 'target'
        if isMid
            mergeSelection(typeTag, responseTag, [900, intStartList]);
        else
            mergeSelection(typeTag, responseTag, [intStartList, 900]);
        end
    otherwise
        if isMid
            mergeSelection_Mid(typeTag, responseTag, intStartList, nPart);
        else
            mergeSelection(typeTag, responseTag, intStartList, nPart);
        end
        %{
        for intStart = intStartList
            fprintf('Target pressure: %d\n', intStart);
            verticalSelection = strcat('Relative', num2str(intStart)); %'MidMeso';%'MidMeso';%'UpperOcean';%'Mikael';
            mergeSelection(typeTag, responseTag, verticalSelection, nPart);
        end
        %}
end
