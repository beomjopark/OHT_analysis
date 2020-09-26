

  [latGrid,longGrid] = meshgrid(linspace(-89.5,89.5,180), linspace(20.5,379.5,360));
  nGrid = numel(latGrid);
gridList = int64(linspace(1, nGrid, 100));
iGrid = gridList(40)


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
          windowSizeMargined = windowSizeKrig * 2;
        otherwise
          windowTypeTag = []
          windowSizeMargined = windowSizeKrig;
    end

    if isempty(isAdjusted)
      isAdjusted = false;
    end
    if isempty(isAbsolute)
      isAbsolute = false;
    end

    if isAdjusted
      adjustTag = 'Adjusted';
      adjustNumTag = ['Adjusted', num2str(nAdjust)];
      windowSizeTag = windowSizeFullTag;
    else
      adjustTag = [];
      adjustNumTag = [];
    end

    if isAbsolute
      absoluteTag = 'Absolute';
    else
      absoluteTag = [];
    end   

    if isempty(iterEM) || (iterEM == 0)
        EMTag = []
    else
        EMTag = ['EM',num2str(iterEM)]
    end

    isFullMonth = (numel(month) == 12)
    if isFullMonth
      EMOutTag = ['Full', EMTag] %% For anomaly, add Full even for 0 iterEM
    else
      EMOutTag = EMTag;
    end
    if isFullMonth && ~isempty(EMTag)
        EMTag = ['Full', EMTag]
    end

    destFolder = ['anomaly_',typeTag,responseTag,adjustNumTag,absoluteTag,windowTypeTag,'_w',windowSizeFullTag,'_',kernelType,'_',verticalSelection,'Season_',num2str(month,'%02d'),EMOutTag,'Kernel'];

    yearRange = strsplit(dataYear, '_');
    startYear = str2num(yearRange{2});
    endYear = str2num(yearRange{3});

figure;

    set(gcf,'units','centimeters')
    set(gcf,'pos',[0 0 22.5 15])
    set(gcf,'paperunits',get(gcf,'units')) 
    set(gcf,'paperpos',get(gcf,'pos'))

    pos = get(gcf,'position');
    set(gcf,'position',pos + [0 0 0 -4])

    for iYear = startYear:endYear
        for iMonth = 1:12
        clf;

        handle = worldmap('World');
        setm(handle, 'Origin', [0 200 0]);
        tightmap;
        mlabel('off');
        plabel('off');
        
        pos = get(gca,'position');
        set(gca,'position',pos + [-0.13 -0.025 0.2 0])

        predGridAccum = zeros(size(latGrid));
        for iGrid = gridList
            if isDeriv
                srcName = ['./Results/',destFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,targetVar,'Deriv_',...
                    num2str(iMonth,'%02d'),'_',num2str(iYear),'_',num2str(iGrid),'.mat'];
            else
                srcName = ['./Results/',destFolder,'/anomaly',responseTag,verticalSelection,dataYear,'SeasonSpaceTime',kernelType,'_',...
                                num2str(iMonth,'%02d'),'_',num2str(iYear),'_',num2str(iGrid),'.mat'];
            end
            if exist(srcName)
                load(srcName);
                predGridAccum = predGridAccum + predGrid;
            end
        end
        
        predGridAccum(predGridAccum == 0) = NaN;

        % Determine ylim from the first Year/Month
%{
        if (iYear == startYear) && (iMonth == 1)
           ylimit = quantile(mask(:).*predLatGrid(:), [0.005 0.995]);
        end

%}
        surfm(latGrid,longGrid,predGridAccum); 
        
        if isDeriv
        title([targetVar,' Anomaly Smoother: ',...
            num2str(iMonth,'%02d'),'/',num2str(iYear)]);
        caxis([0, 0.8]);
        else
        title(['Anomaly Smoother: ',...
            num2str(iMonth,'%02d'),'/',num2str(iYear)]);
                caxis([0, 1]);
        end

        load coast;
        plotm(lat,long,'k');
        
        cLims = caxis;  
        cb = colorbar;
        cb.FontSize = 14;
        drawnow;
        
        M((iYear-2007)*12+iMonth) = getframe(gcf);
%        print('-depsc2',['./Figures/',destFolder,'/','Anomaly_',num2str(iYear,'%02d'),'_',num2str(iMonth, '%02d'),'.eps']);

%        pause(0.01);
        
        end
    end

mkdir(['./Figures/',destFolder])
if isDeriv
    saveName = [responseTag,verticalSelection,dataYear,'_anom_movie_',targetVar]
else
    saveName = [responseTag,verticalSelection,dataYear,'_anom_movie']
end
myVideo = VideoWriter(['./Figures/',destFolder,'/',saveName,'.avi']);

    myVideo.FrameRate = 5;
    myVideo.Quality = 100;
    open(myVideo);
    writeVideo(myVideo, M);
    close(myVideo);