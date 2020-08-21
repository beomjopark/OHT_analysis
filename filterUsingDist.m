function filterUsingDist(typeTag, responseTag, verticalSelection)
    tag = 'PchipPotTemp';
    dataYear = '_2007_2018';
    
    % Load Data
    if isnumeric(verticalSelection)
        % Specifically for target / intlatlon
        % Multiple pressure inputs
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
        saveName = ['./Data/',typeTag,responseTag,'Prof',tag,presString,dataYear,'.mat']
        load(saveName);
    else
        presString = verticalSelection;
        switch responseTag
            case 'Flux'
                saveName = ['./Data/',typeTag,responseTag,'Prof',verticalSelection,dataYear,'.mat']
                load(saveName);
            otherwise
                saveName = ['./Data/',typeTag,responseTag,'Prof',tag,verticalSelection,dataYear,'.mat']
                load(saveName);
        end
    end

%    intStart = intStart;
%    intEnd = intEnd;

    nProfile = numel(profLatAggrSel);
    profYearAggrSel = zeros(1,nProfile);
    profMonthAggrSel = zeros(1,nProfile);
    for iProf = 1:nProfile
        temp = datevec(profJulDayAggrSel(iProf));
        profYearAggrSel(iProf) = temp(1);
        profMonthAggrSel(iProf) = temp(2);
    end


    min_st_thresh = 1e-06 % QC1 threshold
    disp(['QC1: SP dist < ',num2str(min_st_thresh)])

    profSpatioTempSel = [profLatAggrSel; profLongAggrSel; profJulDayAggrSel];
    drop_idx_total = [];
    for iYear = min(profYearAggrSel):max(profYearAggrSel)
        for iMonth = 1:12
            isSel = (profYearAggrSel == iYear & profMonthAggrSel == iMonth);
            selIdx = find(isSel);
            
            profDist = squareform(pdist(profSpatioTempSel(:,isSel)'));
            profDist(1:sum(isSel)+1:end) = Inf;
%            nClose = sum(profDist < min_st_thresh, 'all');

            [dup_idx, ~] = find(profDist < min_st_thresh);

            if ~isempty(dup_idx)
                dup_sel_conv = selIdx(dup_idx);
                drop_idx_total = [drop_idx_total, dup_sel_conv(2:end)];
            end
        end
    end
                    
    figure;
    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
                
    handle = worldmap('World');
    setm(handle, 'Origin', [0 200 0]);
    tightmap;
    mlabel('off');
    plabel('off');
            
    scatterm(profLatAggrSel(drop_idx_total),profLongAggrSel(drop_idx_total),10,'r', '+'); 
    load coast;
    plotm(lat,long,'k');
    hold 
    title('Spatio-Temporally close');
    drawnow;
    print('-depsc2',['./Figures/','QC1_',typeTag,responseTag,'SpatioTemp',num2str(min_st_thresh),'.eps']);
    
    
    % FILTER
    profLatAggrSel(drop_idx_total) = [];
    profLongAggrSel(drop_idx_total) = [];
    profJulDayAggrSel(drop_idx_total) = [];
    profYearAggrSel(drop_idx_total) = [];
    profMonthAggrSel(drop_idx_total) = [];
    profFloatIDAggrSel(drop_idx_total) = [];
    switch typeTag
        case 'int'
            targetDynhProf(:, drop_idx_total) = [];
            intTempProf(:, drop_idx_total) = [];
            intDensProf(:, drop_idx_total) = [];
            startPresIdx(drop_idx_total) = [];
        case 'target'
            switch responseTag
                case 'Temp'
                    targetTempProf(:, drop_idx_total) = [];
                case 'Dens'
                    targetDensProf(:, drop_idx_total) = [];
                otherwise
                    targetDynhProf(:, drop_idx_total) = [];
            end
    end


    min_temp_thres = 1e-02 % QC2 threshold
    disp(['QC2: Spatially exact but time within ', num2str(min_temp_thres)])
                
    figure;
    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))
    
    drop_idx_total = [];
    cnt = 0;
    for iYear = min(profYearAggrSel):max(profYearAggrSel)
        for iMonth = 1:12
            isSel = (profYearAggrSel == iYear & profMonthAggrSel == iMonth);
            selIdx = find(isSel);

            %dd_spatio = squareform(pdist([profLatAggrSel; profLongAggrSel]'));
            checkMat = (squareform(pdist(profJulDayAggrSel(isSel)')) <= min_temp_thres); % Temporal threshold about 15min
            checkMat(1:size(checkMat, 1)+1:end) = false;

            % Find matching spatial location
            [~,~,ic] = unique([profLatAggrSel(isSel)',profLongAggrSel(isSel)'], 'rows');

            [n_ic, ic_edges] = histcounts(ic, 1:(max(ic)+1));
            ic_edges(end) = [];
            
            % Find close time profile within spatial match
            drop_idx = [];
            for row_unique = ic_edges(n_ic > 1)
                is_dup_rows = (ic == row_unique);
                dup_rows_idx = find(is_dup_rows);
                
                is_temp_not_close = (sum(checkMat(is_dup_rows, is_dup_rows)) == 0);
                if sum(is_temp_not_close) == 0
                    % every profile is close enough: keep the first
                    drop_row_idx = dup_rows_idx(2:end);
                else
                    % Exists Non trivial time close: keep them + first nonunique
                    if sum(is_temp_not_close) ~= numel(is_temp_not_close)
                        sel_dup_idx = find(~ is_temp_not_close);
                        is_temp_not_close(sel_dup_idx(1)) = true;
                    end
                    drop_row_idx = dup_rows_idx(~is_temp_not_close);
                end
                drop_idx = [drop_idx; drop_row_idx];
            end
            drop_idx_total = [drop_idx_total; selIdx(drop_idx)'];

            if ~isempty(drop_idx)
                cnt = cnt + 1;
                disp([num2str(iYear), '/', num2str(iMonth), ': ', num2str(numel(drop_idx))])
                
                handle = worldmap('World');
                setm(handle, 'Origin', [0 200 0]);
                tightmap;
                mlabel('off');
                plabel('off');
            
                scatterm(profLatAggrSel(selIdx(drop_idx)),profLongAggrSel(selIdx(drop_idx)),10,'r', '+'); 
%                diff(profJulDayAggrSel(drop_idx))
                load coast;
                plotm(lat,long,'k');
                hold off
                title([num2str(iYear), '/', num2str(iMonth)])
                drawnow;

                M(cnt) = getframe(gcf);
            end
        end
    end
    myVideo = VideoWriter(['./Figures/','QC2_',typeTag,responseTag,'Temp',num2str(min_temp_thres),'.avi']);    
    myVideo.FrameRate = 5;
    myVideo.Quality = 100;
    open(myVideo);
    writeVideo(myVideo, M);
    close(myVideo);   
    
    % FILTER
    profLatAggrSel(drop_idx_total) = [];
    profLongAggrSel(drop_idx_total) = [];
    profJulDayAggrSel(drop_idx_total) = [];
    profYearAggrSel(drop_idx_total) = [];
    profMonthAggrSel(drop_idx_total) = [];
    profFloatIDAggrSel(drop_idx_total) = [];
    switch typeTag
        case 'int'
            targetDynhProf(:, drop_idx_total) = [];
            intTempProf(:, drop_idx_total) = [];
            intDensProf(:, drop_idx_total) = [];
            save(saveName,...
                'targetDynhProf', 'intDensProf', 'intTempProf', 'profFloatIDAggrSel', 'profJulDayAggrSel', 'profLongAggrSel', 'profLatAggrSel', 'profYearAggrSel', 'startPresIdx', 'intEnd', 'intStart');
        case 'target'
            switch responseTag
                case 'Temp'
                    targetTempProf(:, drop_idx_total) = [];
                    save(saveName,...
                            'targetTempProf', 'profFloatIDAggrSel', 'profJulDayAggrSel', 'profLongAggrSel', 'profLatAggrSel', 'startPresIdx', 'targetPres');
                case 'Dens'
                    targetDensProf(:, drop_idx_total) = [];
                    save(saveName,...
                            'targetDensProf', 'profFloatIDAggrSel', 'profJulDayAggrSel', 'profLongAggrSel', 'profLatAggrSel', 'startPresIdx', 'targetPres');
                otherwise
                    targetDynhProf(:, drop_idx_total) = [];
                    save(saveName,...
                            'targetDynhProf', 'profFloatIDAggrSel', 'profJulDayAggrSel', 'profLongAggrSel', 'profLatAggrSel', 'startPresIdx', 'targetPres' ,'refPres');
            end
    end

    
end

%{
month = 2;
startYear = min(profYearAggrSel);
endYear = max(profYearAggrSel);
adjustNumTag = [];
fitMLE = load(['./Results/','prethresh','/localMLESpaceTime',kernelType,responseTag,'Relative',num2str(targetPres(presIdx)),'Season_',num2str(month,'%02d'),'_',num2str(startYear),'_',num2str(endYear),adjustNumTag,absoluteTag,'_w',num2str(windowSize),'.mat']);

profLatAggrSelRounded = roundHalf(profLatAggrSel);
profLongAggrSelRounded = roundHalf(profLongAggrSel);  
profScaled = zeros(3, nProfile);
parfor iProfile = 1:nProfile
    profIdx = (fitMLE.latGrid == profLatAggrSelRounded(iProfile) & fitMLE.longGrid == profLongAggrSelRounded(iProfile));
    profScaled(1, iProfile) = profLatAggrSel(iProfile) ./ fitMLE.thetaLatOpt(profIdx);
    profScaled(2, iProfile) = profLongAggrSel(iProfile) ./ fitMLE.thetaLongOpt(profIdx);
    profScaled(3, iProfile) = profJulDayAggrSel(iProfile) ./ fitMLE.thetatOpt(profIdx);
end
%}
