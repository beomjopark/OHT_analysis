function extendedDataSeason(meanTag, typeTag, responseTag, verticalSelection, dataYear, minNumberOfObs, is2step, isAdjusted, isAbsolute, nAdjust)
%% Extend monthly data into 3month data by aggregating +-1 month.
%% Always run after divideDataToMonthsSeason
    if isempty(isAdjusted)
      isAdjusted = false;
    end
    if isempty(isAbsolute)
      isAbsolute = false;
    end

    if isAdjusted | nAdjust > 0
      adjustTag = 'Adjusted';
      adjustNumTag = [adjustTag, num2str(nAdjust)];
    else
      adjustTag = [];
      adjustNumTag = adjustTag;
    end

    if isAbsolute
      absoluteTag = 'Absolute';
    else
      absoluteTag = [];
    end   
    
    yearRange = strsplit(dataYear, '_');
    startYear = str2double(yearRange{2});
    endYear = str2double(yearRange{3});
    for iYear = startYear:endYear
        for iMonth = 1:12       
            % Data for iMonth-1
            if iMonth == 1
                if iYear == startYear
                    S1.profLatAggrMonth = [];
                    S1.profLongAggrMonth = [];
                    S1.profJulDayAggrMonth = [];
                    if ~is2step
                    S1.profFloatIDAggrMonth = [];
                    end
                    switch responseTag
                        case 'Flux'
                            switch typeTag
                                case 'lat'
                                    S1.latFluxResMonth = [];
                                case 'lon'
                                    S1.lonFluxResMonth = [];
                                otherwise
                                    S1.intFluxResMonth = [];
                            end
                        case 'Temp'
                            S1.targetTempResMonth = [];
                        case 'Dens'
                            S1.intDensResMonth = [];
                        case 'DUACS'
                            S1.intDensResMonth = [];  %NAME CONVENTION for LOCALMLE
                        case 'Sal'
                            S1.targetSalResMonth = [];
                    end
                else
                    S1 = load(['./Data/Monthly/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(12,'%02d'),'_',num2str(iYear-1),'.mat']);
                end
            else
                S1 = load(['./Data/Monthly/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(iMonth-1,'%02d'),'_',num2str(iYear),'.mat']);
            end
            
            % Data for iMonth
            S2 = load(['./Data/Monthly/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(iMonth,'%02d'),'_',num2str(iYear),'.mat']);
            
            % Data for iMonth+1
            if iMonth == 12
                if iYear == endYear
                    S3.profLatAggrMonth = [];
                    S3.profLongAggrMonth = [];
                    S3.profJulDayAggrMonth = [];
                    if ~is2step
                        S3.profFloatIDAggrMonth = [];
                    end
                    switch responseTag
                        case 'Flux'
                            switch typeTag
                                case 'lat'
                                    S3.latFluxResMonth = [];
                                case 'lon'
                                    S3.lonFluxResMonth = [];
                                otherwise
                                    S3.intFluxResMonth = [];
                            end
                        case 'Temp'
                            S3.targetTempResMonth = [];
                        case 'Dens'
                            S3.intDensResMonth = [];
                        case 'DUACS'
                            S3.intDensResMonth = [];
                        case 'Sal'
                            S3.targetSalResMonth = [];
                    end
                else
                    S3 = load(['./Data/Monthly/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(1,'%02d'),'_',num2str(iYear+1),'.mat']);
                end
            else
                S3 = load(['./Data/Monthly/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(iMonth+1,'%02d'),'_',num2str(iYear),'.mat']);
            end
            
            profLatAggr3Months = [S1.profLatAggrMonth S2.profLatAggrMonth S3.profLatAggrMonth];
            profLongAggr3Months = [S1.profLongAggrMonth S2.profLongAggrMonth S3.profLongAggrMonth];
            profJulDayAggr3Months =  [S1.profJulDayAggrMonth S2.profJulDayAggrMonth S3.profJulDayAggrMonth];
%{
            if ~is2step
                profFloatIDAggr3Months = [S1.profFloatIDAggrMonth S2.profFloatIDAggrMonth S3.profFloatIDAggrMonth];
            end
%}

            
            fprintf("%d observations at %d/%d \n", numel(profLatAggr3Months), iYear, iMonth);

            switch responseTag
                case 'Flux'
                    switch typeTag
                        case 'lat'
                            latFluxRes3Months = [S1.latFluxResMonth S2.latFluxResMonth S3.latFluxResMonth];
                            save(['./Data/Extended/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(iMonth,'%02d'),'_',num2str(iYear),'_extended.mat'],...
                                'latFluxRes3Months','profLatAggr3Months','profLongAggr3Months','profJulDayAggr3Months');
                        case 'lon'
                            lonFluxRes3Months = [S1.lonFluxResMonth S2.lonFluxResMonth S3.lonFluxResMonth];
                            save(['./Data/Extended/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(iMonth,'%02d'),'_',num2str(iYear),'_extended.mat'],...
                                'lonFluxRes3Months','profLatAggr3Months','profLongAggr3Months','profJulDayAggr3Months');
                        otherwise
                            intFluxRes3Months = [S1.intFluxResMonth S2.intFluxResMonth S3.intFluxResMonth];
                            save(['./Data/Extended/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(iMonth,'%02d'),'_',num2str(iYear),'_extended.mat'],...
                                'intFluxRes3Months','profLatAggr3Months','profLongAggr3Months','profJulDayAggr3Months');                            
                    end
                case 'Temp'
                    targetTempRes3Months = [S1.targetTempResMonth S2.targetTempResMonth S3.targetTempResMonth];
                    save(['./Data/Extended/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(iMonth,'%02d'),'_',num2str(iYear),'_extended.mat'],...
                        'targetTempRes3Months','profLatAggr3Months','profLongAggr3Months','profJulDayAggr3Months');
                case 'Dens'
                    intDensRes3Months = [S1.intDensResMonth S2.intDensResMonth S3.intDensResMonth];
                    save(['./Data/Extended/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(iMonth,'%02d'),'_',num2str(iYear),'_extended.mat'],...
                        'intDensRes3Months','profLatAggr3Months','profLongAggr3Months','profJulDayAggr3Months');
                case 'DUACS'
                    intDensRes3Months = [S1.intDensResMonth S2.intDensResMonth S3.intDensResMonth];
                    save(['./Data/Extended/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(iMonth,'%02d'),'_',num2str(iYear),'_extended.mat'],...
                        'intDensRes3Months','profLatAggr3Months','profLongAggr3Months','profJulDayAggr3Months');
                case 'Sal'
                    % For convenience
                    intDensRes3Months = [S1.targetSalResMonth S2.targetSalResMonth S3.targetSalResMonth];
                    save(['./Data/Extended/',typeTag,responseTag,'Res',verticalSelection,dataYear,adjustNumTag,absoluteTag,'SeasonMonth_',num2str(iMonth,'%02d'),'_',num2str(iYear),'_extended.mat'],...
                        'intDensRes3Months','profLatAggr3Months','profLongAggr3Months','profJulDayAggr3Months');
            end
            
        end
    end

end
