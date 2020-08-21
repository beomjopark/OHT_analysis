function mergeSelection_Mid(typeTag, responseTag, verticalSelection, nPart)
    %% Merge the Updated Data(_2017_2018) with legacy
    tag = 'PchipPotTemp';
    dataYear = '_2017_2018';

    if ~isnumeric(verticalSelection)
        switch verticalSelection
            case 'MidMeso'
                targetPres = 600;
            case 'FullMeso'
                targetPres = 200;
        end
        presString = verticalSelection;
    else
        targetPres = verticalSelection;
        nTargetPres = numel(targetPres);
        if nTargetPres > 1
            presString = [num2str(min(targetPres)),'_',num2str(max(targetPres))];
        else
            presString = num2str(targetPres);
        end

        if strcmp(typeTag, 'int')
            presString = ['Relative', presString];
        end
    end

    switch typeTag
        case 'int'
            if nargin > 3 % exists nPart
                S1 = load(['./Data/',typeTag,'TempDensProf',tag,presString,'Part',num2str(1),'_',num2str(nPart(1)),'.mat']);
                for iPart = 2:nPart(1)
                    S1part = load(['./Data/',typeTag,'TempDensProf',tag,presString,'Part',num2str(iPart),'_',num2str(nPart(1)),'.mat']);
                    S1.profLatAggrSel = [S1.profLatAggrSel S1part.profLatAggrSel];
                    S1.profLongAggrSel = [S1.profLongAggrSel S1part.profLongAggrSel];
                    S1.profYearAggrSel = [S1.profYearAggrSel S1part.profYearAggrSel];
                    S1.profJulDayAggrSel = [S1.profJulDayAggrSel S1part.profJulDayAggrSel];
                    S1.profFloatIDAggrSel = [S1.profFloatIDAggrSel S1part.profFloatIDAggrSel];

                    S1.targetDynhProf = [S1.targetDynhProf S1part.targetDynhProf];
                    S1.intTempProf = [S1.intTempProf S1part.intTempProf];
                    S1.intDensProf = [S1.intDensProf S1part.intDensProf];

                    S1.endPresIdx = [S1.endPresIdx S1part.endPresIdx];
                end
            else
                S1 = load(['./Data/',typeTag,'TempDensProf',tag,presString,'.mat']);
            end

            if numel(nPart > 1)
                S2 = load(['./Data/',typeTag,'TempDensProf',tag,presString,dataYear,'Part',num2str(1),'_',num2str(nPart(2)),'.mat']);
                for iPart = 2:nPart(2)
                    S2part = load(['./Data/',typeTag,'TempDensProf',tag,presString,dataYear,'Part',num2str(iPart),'_',num2str(nPart(2)),'.mat']);
                    S2.profLatAggrSel = [S2.profLatAggrSel S2part.profLatAggrSel];
                    S2.profLongAggrSel = [S2.profLongAggrSel S2part.profLongAggrSel];
                    S2.profYearAggrSel = [S2.profYearAggrSel S2part.profYearAggrSel];
                    S2.profJulDayAggrSel = [S2.profJulDayAggrSel S2part.profJulDayAggrSel];
                    S2.profFloatIDAggrSel = [S2.profFloatIDAggrSel S2part.profFloatIDAggrSel];

                    S2.targetDynhProf = [S2.targetDynhProf S2part.targetDynhProf];
                    S2.intTempProf = [S2.intTempProf S2part.intTempProf];
                    S2.intDensProf = [S2.intDensProf S2part.intDensProf];

                    S2.endPresIdx = [S2.endPresIdx S2part.endPresIdx];
                end                
            else
                S2 = load(['./Data/',typeTag,'TempDensProf',tag,presString,dataYear,'.mat']);
            end
            targetDynhProf = [S1.targetDynhProf S2.targetDynhProf];
            intDensProf = [S1.intDensProf S2.intDensProf];
            intTempProf = [S1.intTempProf S2.intTempProf];
            profFloatIDAggrSel = [S1.profFloatIDAggrSel S2.profFloatIDAggrSel];
            profJulDayAggrSel = [S1.profJulDayAggrSel S2.profJulDayAggrSel];
            profLatAggrSel = [S1.profLatAggrSel S2.profLatAggrSel];
            profLongAggrSel = [S1.profLongAggrSel S2.profLongAggrSel];
            profYearAggrSel = [S1.profYearAggrSel S2.profYearAggrSel];

            endPresIdx = [S1.endPresIdx S2.endPresIdx];

            intStart = S1.intStart;
            intEnd = S1.intEnd;

        case 'target'
            S1 = load(['./Data/',typeTag,responseTag,'Prof',tag,presString,'.mat']);
            S2 = load(['./Data/',typeTag,responseTag,'Prof',tag,presString,dataYear,'.mat']);

            switch responseTag
                case 'Temp'
                    targetTempProf = [S1.targetTempProf S2.targetTempProf];
                case 'Dens'
                    targetDensProf = [S1.targetDensProf S2.targetDensProf];
                otherwise % Dynh
                    targetDynhProf = [S1.targetDynhProf S2.targetDynhProf];
            end
            profFloatIDAggrSel = [S1.profFloatIDAggrSel S2.profFloatIDAggrSel];
            profJulDayAggrSel = [S1.profJulDayAggrSel S2.profJulDayAggrSel];
            profLatAggrSel = [S1.profLatAggrSel S2.profLatAggrSel];
            profLongAggrSel = [S1.profLongAggrSel S2.profLongAggrSel];

            endPresIdx = [S1.endPresIdx S2.endPresIdx];
    end


    %% Postprocess bad quality data
    badFloatID = 2901179;
    dropIdx = (profFloatIDAggrSel == badFloatID);
    profLatAggrSel = profLatAggrSel(~dropIdx);
    profLongAggrSel = profLongAggrSel(~dropIdx);
    profJulDayAggrSel = profJulDayAggrSel(~dropIdx);
    profFloatIDAggrSel = profFloatIDAggrSel(~dropIdx);
    endPresIdx = endPresIdx(~dropIdx);
    switch typeTag
        case 'int'
            profYearAggrSel = profYearAggrSel(~dropIdx);
            targetDynhProf = targetDynhProf(:, ~dropIdx);
            intTempProf = intTempProf(:, ~dropIdx);
            intDensProf = intDensProf(:, ~dropIdx);
            save(['./Data/intTempDensProf',tag,presString,...
                '_',num2str(min(profYearAggrSel)),'_',num2str(max(profYearAggrSel)),'.mat'],...
                'targetDynhProf', 'intDensProf', 'intTempProf', 'profFloatIDAggrSel', 'profJulDayAggrSel', 'profLongAggrSel', 'profLatAggrSel', 'profYearAggrSel', 'endPresIdx', 'intEnd', 'intStart');

        case 'target'
            profYearAggrSel = datevec(profJulDayAggrSel);
            profYearAggrSel = profYearAggrSel(:,1);
            saveName = ['./Data/',typeTag,responseTag,'Prof',tag,presString,...
                            '_',num2str(min(profYearAggrSel)),'_',num2str(max(profYearAggrSel)),'.mat'];
            switch responseTag
                case 'Temp'
                    targetTempProf = targetTempProf(:, ~dropIdx);
                    save(saveName,...
                            'targetTempProf', 'profFloatIDAggrSel', 'profJulDayAggrSel', 'profLongAggrSel', 'profLatAggrSel', 'endPresIdx', 'targetPres');
                case 'Dens'
                    targetDensProf = targetDensProf(:, ~dropIdx);
                    save(saveName,...
                            'targetDensProf', 'profFloatIDAggrSel', 'profJulDayAggrSel', 'profLongAggrSel', 'profLatAggrSel', 'endPresIdx', 'targetPres');
                otherwise
                    refPres = S1.refPres;
                    targetDynhProf = targetDynhProf(:, ~dropIdx);
                    save(saveName,...
                            'targetDynhProf', 'profFloatIDAggrSel', 'profJulDayAggrSel', 'profLongAggrSel', 'profLatAggrSel', 'endPresIdx', 'targetPres' ,'refPres');
            end
    end

end
