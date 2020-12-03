%% Wrapper function for localMLE and Anomaly Prediction
%% Run after selection / Merge and MeanField
addpath(genpath('../gsw_matlab'));
addpath(genpath('./cbrewer'));
addpath(genpath('./Util'));


%% Run Selection
%intStartList = [10, 15, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800] % Full upper ocean
%intStartList = [intStartList, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1850, 1900] % Mid Ocean
isPlot = false

if is2step && strcmp(typeTag, 'int') % This case is for intlatflux/intlonflux
    fprintf('Target pressure: %d to %d\n', min(intStartList), max(intStartList));
    typeTag = strcat(typeTag, targetVar); %'intlatlon'
    verticalSelection = strcat(num2str(min(intStartList)),'_',num2str(max(intStartList)));

    computeMeanAnomalies(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM);

else
    for intStartIdx = 1:numel(intStartList)
        intStart = intStartList(intStartIdx); % Parfor requires increasing 1 index

        fprintf('Target pressure: %d\n', intStart);
        verticalSelection = strcat('Relative', num2str(intStart)); %'MidMeso';%'MidMeso';%'UpperOcean';%'Mikael';

        if is2step  % This case is for each latflux / lonflux
            computeMeanAnomalies(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM);
        else
%            computeMeanAnomalies(kernelType, month, typeTag, responseTag, verticalSelection, dataYear, meanTag, windowType, windowSize, minNumberOfObs, is2step, isDeriv, targetVar, [], [], isAdjusted, isAbsolute, nAdjust, iterEM);

          if isDeriv
            error("Debiasing works only for non-derivative!")
            %% Run this part after you have both mean Anomaly of lat and lon
%               plotMeanFieldAdjustedAbsolute(month, meanTag, kernelType, typeTag, responseTag, verticalSelection, dataYear, windowSize, minNumberOfObs, is2step, fluxType, eqBorder, nAdjust);
          else
            adjustResidual(month, kernelType, meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, isPlot, isStandardize, [], [], nAdjust, iterEM, isFullMonth);
           
           %% isAdjust = True, nAdjust+1
             divideDataToMonthsSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, [], [], true, false, nAdjust+1, iterEM, isFullMonth);

            extendedDataSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, minNumberOfObs, is2step, true, false, nAdjust+1);

            % Create isAdjust version of datamask as a copy
            switch numel(windowSize)
                case 3
                  windowSizeMean = windowSize(1);
                case 2
                  windowSizeMean = windowSize(1);
                 otherwise
                  windowSizeMean = windowSize;
            end
            srcMaskName = ['./Data/dataMask',verticalSelection,dataYear,[],[],'_',num2str(minNumberOfObs),windowType,'_w',num2str(windowSizeMean),'.mat'];
            descMaskName = ['./Data/dataMask',verticalSelection,dataYear,'Adjusted',[],'_',num2str(minNumberOfObs),windowType,'_w',num2str(windowSizeMean),'.mat'];
            system(['cp ', srcMaskName, ' ', descMaskName]);
            fprintf('Refit LocalMLE afterward!');
          end
        end
        close all
    end
end