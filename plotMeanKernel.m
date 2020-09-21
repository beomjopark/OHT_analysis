
nHar = 6;
   iLat = 1 + (2*nHar) + 1;
    iLong = iLat + 1;
    iLatLong = iLat + 2;
    iLat2 = iLat + 3;
    iLong2 = iLat + 4;
    
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
      windowSizeMargined = windowSizeMean * 2;
    otherwise
      windowTypeTag = []
      windowSizeMargined = windowSizeMean;
  end


  % Additional Parameter
  tag = 'PchipPotTemp';
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
    prevEMTag = []
  elseif iterEM == 1
    EMTag = ['EM',num2str(iterEM)]
    prevEMTag = []
  else
    EMTag = ['EM',num2str(iterEM)]
    prevEMTag = ['EM',num2str(iterEM-1)]
  end


  isFullMonth = (numel(month) == 12)
  
  if isFullMonth && ~isempty(EMTag)
    EMTag = ['Full', EMTag]
  end
  if isFullMonth && ~isempty(prevEMTag)
    prevEMTag = ['Full', prevEMTag]
  end
 

 EM = load(['./Results/meanField',responseTag,meanTag,tag,verticalSelection,dataYear,adjustTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_',num2str(minNumberOfObs),EMTag,'.mat']);
 EMKernel = load(['./Results/meanFieldKernel',responseTag,meanTag,tag,verticalSelection,dataYear,adjustTag,absoluteTag,windowTypeTag,'_w',windowSizeTag,'_',num2str(minNumberOfObs),EMTag,'.mat']);

   if is2step
    destFolder = [typeTag,fluxType,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,'_Eq',num2str(eqBorder)];
  else
    destFolder = [typeTag,responseTag,verticalSelection,dataYear,adjustTag,absoluteTag,meanTag];
   end
   destFolder = [destFolder, EMTag, 'Kernel']
  mkdir(['./Figures/',destFolder])

 
 for ii = 1:size(EMKernel.betaGrid,3)
 figure;
  handle = worldmap('World');
  setm(handle, 'Origin', [0 200 0]);
  tightmap;
  mlabel('off');
  plabel('off');

  surfm(EM.latGrid,EM.longGrid,mask.*EMKernel.betaGrid(:,:,ii));

  load coast;
  plotm(lat,long,'k');

  h = colorbar;

  colorbar;

  temp = squeeze(mask.*EMKernel.betaGrid(:,:,ii));
  caxis([quantile(temp(:), 0.03), quantile(temp(:), 0.95)]);
%  caxis([quantile(temp(:), 0.05), quantile(temp(:), 0.95)]);

  cLims = caxis;
  colormap(darkb2r(cLims(1),cLims(2)));
  %cLims = caxis;
  %colormap(darkb2r(cLims(1),cLims(2)));
  if cLims(1) > 0 || cLims(2) <0
          colormap(parula(100));
  end
  
  drawnow;

  set(gcf,'units','centimeters')
  set(gcf,'pos',[0 0 22.5 15])
  set(gcf,'paperunits',get(gcf,'units')) 
  set(gcf,'paperpos',get(gcf,'pos'))
  
  switch ii
      case 1
          betaType = '0'
      case iLat
          betaType = 'Lat'
      case iLong
          betaType = 'Long'
      case iLatLong
          betaType = 'Latlong'
      case iLat2
          betaType = 'Lat2'
      case iLong2
          betaType = 'Long2'
      otherwise
          betaType = num2str(ii-1)
  end
  print('-depsc2',['./Figures/',destFolder,'/','betaKernel',betaType,'_',responseTag,meanTag,tag,verticalSelection,'_',num2str(windowSize),'.eps']);
 end