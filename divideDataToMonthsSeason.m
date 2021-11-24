function divideDataToMonthsSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, windowType, windowSize, minNumberOfObs, is2step, fluxType, eqBorder, isAdjusted, isAbsolute, nAdjust, iterEM, isFullMonth)
    % Load Data

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

    if windowSizeMean == windowSizeCov
        windowSizeTag = num2str(windowSizeMean)
    else
        windowSizeTag = [num2str(windowSizeMean),'_',num2str(windowSizeCov)]
    end
    windowSizeFullTag = [windowSizeTag, '_', num2str(windowSizeKrig)]

    if isempty(windowType)
        windowType = 'box'; % expected 'spherical'
    end

    switch windowType
        case 'spherical'
          windowTypeTag = 'spherical'
          windowSizeMargined = windowSize * 2;
        otherwise
          windowTypeTag = []
          windowSizeMargined = windowSize;
    end


    if isempty(isAdjusted)
      isAdjusted = false;
    end
    if isempty(isAbsolute)
      isAbsolute = false;
    end

    if isAdjusted | nAdjust > 0
      adjustTag = 'Adjusted';
      adjustNumTag = [adjustTag, num2str(nAdjust)];
      windowSizeTag = windowSizeFullTag;
    else
      adjustTag = [];
      adjustNumTag = adjustTag;
    end

    if isAbsolute
      absoluteTag = 'Absolute';
    else
      absoluteTag = [];
    end    

    if isempty(iterEM) || (iterEM == 0)
        EMTag = []
        prevEMTag = []
    elseif iterEM == 1
        EMTag = ['EM',num2str(iterEM)]
        prevEMTag = []
    else
        EMTag = ['EM',num2str(iterEM)]
        prevEMTag = ['EM',num2str(iterEM-1)]
    end
    
    if isFullMonth && ~isempty(EMTag)
      EMTag = ['Full', EMTag]
    end
    if isFullMonth && ~isempty(prevEMTag)
      prevEMTag = ['Full', prevEMTag]
    end

    if is2step
        srcName = ['./Data/',typeTag, fluxType,responseTag,'Res',verticalSelection,dataYear,adjustTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,'_Eq',num2str(eqBorder),EMTag,'.mat']
%        load(['./Data/',typeTag, fluxType,responseTag,'Res',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'_w',windowSizeTag,'_Eq',num2str(eqBorder),AdjustTag,EMTag,'.mat']);
    else
        %% Absolute is not applicable and for backward-compatibility
        srcName = ['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'Filtered_',num2str(minNumberOfObs),windowTypeTag,'_w',windowSizeTag,EMTag,'.mat']
%        load(['./Data/',typeTag,responseTag,'Res',verticalSelection,dataYear,'Filtered_',num2str(minNumberOfObs),'_w',windowSizeTag,AdjustTag,'.mat']);
    end
    load(srcName);

    nProf = length(profLatAggrSel);
    profYearAggrSel = zeros(1,nProf);
    profMonthAggrSel = zeros(1,nProf);
    for iProf = 1:nProf
        temp = datevec(profJulDayAggrSel(iProf));
        profYearAggrSel(iProf) = temp(1);
        profMonthAggrSel(iProf) = temp(2);
    end

    %% Divide
    yearRange = strsplit(dataYear, '_');
    startYear = str2double(yearRange{2});
    endYear = str2double(yearRange{3});
    for iYear = startYear:endYear
        for iMonth = 1:12
            idx = (profYearAggrSel == iYear & profMonthAggrSel == iMonth);
            profLatAggrMonth = profLatAggrSel(idx);
            profLongAggrMonth = profLongAggrSel(idx);
            profJulDayAggrMonth = profJulDayAggrSel(idx);

            saveName = ['./Data/Monthly/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(iMonth,'%02d'),'_',num2str(iYear),'.mat'];
            if is2step
                switch typeTag
                    case {'intlat', 'intlon'}
                        intFluxResMonth = intFluxRes(idx);
                        save(saveName,...
                      'profLatAggrMonth','profLongAggrMonth','profJulDayAggrMonth','intFluxResMonth');
                    otherwise
                        FluxResMonth = FluxRes(idx);
                        if isempty(targetTempProf)
                            save(saveName,...
                                  'profLatAggrMonth','profLongAggrMonth','profJulDayAggrMonth','FluxResMonth');
                        else
                            targetTempProfMonth = targetTempProf(idx);
                            save(saveName,...
                                  'profLatAggrMonth','profLongAggrMonth','profJulDayAggrMonth','FluxResMonth', 'targetTempProfMonth');
                        end
                end
            else        
%                profFloatIDAggrMonth = profFloatIDAggrSel(idx);
                switch responseTag
                    case 'Temp'
                        targetTempResMonth = targetTempRes(idx);
                        save(saveName,...
                            'profLatAggrMonth','profLongAggrMonth','profJulDayAggrMonth','targetTempResMonth');
                    case 'Dens'
                        intDensResMonth = intDensRes(idx);
                        save(saveName,...
                            'profLatAggrMonth','profLongAggrMonth','profJulDayAggrMonth','intDensResMonth');
                    case 'DUACS'
                        intDensResMonth = targetADTRes(idx); % NAMING CONVENTION
                        save(saveName,...
                            'profLatAggrMonth','profLongAggrMonth','profJulDayAggrMonth','intDensResMonth');
                    case 'ESA'
                        intDensResMonth = targetSSTRes(idx); % NAMING CONVENTION
                        save(saveName,...
                            'profLatAggrMonth','profLongAggrMonth','profJulDayAggrMonth','intDensResMonth');                        
                    case 'Sal'
                        targetSalResMonth = targetSalRes(idx);
                        save(saveName,...
                            'profLatAggrMonth','profLongAggrMonth','profJulDayAggrMonth','targetSalResMonth');
                end
            end

            fprintf("%d observations in %d/%d\n", sum(idx), iYear, iMonth);
        end
    end

end
