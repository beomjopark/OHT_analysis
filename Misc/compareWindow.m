%% Compare the impact of the windowsize
addpath('C:\Users\beomjo\Box Sync\Argo\OHT_analysis\Util')
addpath(genpath('C:\Users\beomjo\Box Sync\Argo\gsw_matlab'))

nHar = 6;
iLat = 1 + (2*nHar) + 1;
iLong = iLat + 1;

% Conversion: Deg / m
[latGrid,longGrid] = meshgrid(linspace(-89.5,89.5,180),linspace(20.5,379.5,360));
latDistGrid = 1 ./ distance(latGrid - 0.5, longGrid, latGrid + 0.5, longGrid,...
                        referenceEllipsoid('WGS84', 'm'));
longDistGrid = 1 ./ distance(latGrid, longGrid - 0.5, latGrid, longGrid + 0.5,...
                        referenceEllipsoid('WGS84', 'm'));

% 0. Load estimates
windows = [3, 4, 5]
windowSizes = ["w3", "w4", "w5"]
windowTag = ["w3_3", "w4_4", "w5_5"]
nWindows = numel(windowSizes)
ARGO = cell2struct(cell(nWindows, 1), cellstr(windowSizes));
for iiter = 1:nWindows
    win = windowSizes(iiter);
    ARGO.(win).mean = load(join(['./Results/meanFieldDUACSNoTrendPchipPotTempRelative10_2007_2018spherical_',win,'_20EM5.mat'], ""));
    ARGO.(win).mean.zonAbsVel = (- ARGO.(win).mean.betaGrid(:,:,iLat) .* latDistGrid ./ gsw_f(latGrid)) .* gsw_grav(latGrid);
    ARGO.(win).mean.merAbsVel = ARGO.(win).mean.betaGrid(:,:,iLong) .* longDistGrid ./ gsw_f(latGrid) .* gsw_grav(latGrid);
    ARGO.(win).anomFolder = join(['./Results/anomaly_intDUACSspherical_',windowTag(iiter),'_Matern_Relative10Season_11EM5/'], "");
end

eqBorder = 2%1%3
eqmask = ~(ARGO.(windowSizes(1)).mean.latGrid < eqBorder & ARGO.(windowSizes(1)).mean.latGrid > - eqBorder) .* 1;
eqmask(eqmask == 0) = NaN;



%0.1 Figure
for win = windowSizes
figure
handle = worldmap('World');
setm(handle, 'Origin', [0 200 0]);
tightmap;
mlabel('off');
plabel('off');
surfm(ARGO.(win).mean.latGrid, ARGO.(win).mean.longGrid, ARGO.(win).mean.zonAbsVel .* eqmask); 
ylimit = [-0.3, 0.3]; % raw version
caxis(ylimit);
cLims = caxis;  
colormap(darkb2r(cLims(1),cLims(2)));
colorbar;

load coast;
plotm(lat,long,'k');
title(win)
end


DUACSFolder = 'D:\dataset-duacs-rep-global-merged-allsat-phy-l4';
monList = 1:12
nMon = length(monList)%12;
DUACS = struct('zonVel', struct(), 'merVel', struct(), 'adt', struct());
DUACS.adt.mse = zeros(length(2007:2018) * nMon, nWindows);
DUACS.adt.sdse = zeros(length(2007:2018) * nMon, nWindows);
DUACS.adt.mad = zeros(length(2007:2018) * nMon, nWindows);
DUACS.adt.llk = zeros(length(2007:2018) * nMon, nWindows);

DUACS.zonVel.mse = zeros(length(2007:2018) * nMon, nWindows);
DUACS.zonVel.sdse = zeros(length(2007:2018) * nMon, nWindows);
DUACS.zonVel.mad = zeros(length(2007:2018) * nMon, nWindows);
DUACS.zonVel.llk = zeros(length(2007:2018) * nMon, nWindows);

DUACS.merVel.mse = zeros(length(2007:2018) * nMon, nWindows);
DUACS.merVel.sdse = zeros(length(2007:2018) * nMon, nWindows);
DUACS.merVel.mad = zeros(length(2007:2018) * nMon, nWindows);
DUACS.merVel.llk = zeros(length(2007:2018) * nMon, nWindows);


cnt = 1;
for year = 2007:2018
    disp(num2str(year))
    for mon = monList
        %% Get DUACS absVel & adt grid
        DS = load([DUACSFolder, '\absVel\absVel',num2str(year),num2str(mon, '%02d'),'15_regrid.mat']);
        DSADT = load([DUACSFolder, '\adt\adt',num2str(year),num2str(mon, '%02d'),'15_regrid.mat']);

        %% Get pseudoEstimate
        nanMask = true([size(latGrid), 3]);
        for iiter = 1:nWindows
            iter = windowSizes(iiter);

            % ADT
            S = load(join([ARGO.(iter).anomFolder, 'anomalyDUACSRelative10_2007_2018SeasonSpaceTimeMatern_',num2str(mon, '%02d'),'_',num2str(year),'.mat'], ''));
            JulDayRatio = fromJulDayToYearDay(S.midJulDay) ./ yearLength(S.midJulDay);
            ARGO.(iter).fulladt = ARGO.(iter).mean.betaGrid(:,:,1) ...
                         + ARGO.(iter).mean.betaGrid(:,:,2) .* sin(2*pi*1*JulDayRatio)' + ARGO.(iter).mean.betaGrid(:,:,3) .* cos(2*pi*1*JulDayRatio)' ...
                         + ARGO.(iter).mean.betaGrid(:,:,4) .* sin(2*pi*2*JulDayRatio)' + ARGO.(iter).mean.betaGrid(:,:,5) .* cos(2*pi*2*JulDayRatio)' ...
                         + ARGO.(iter).mean.betaGrid(:,:,6) .* sin(2*pi*3*JulDayRatio)' + ARGO.(iter).mean.betaGrid(:,:,7) .* cos(2*pi*3*JulDayRatio)' ...
                         + ARGO.(iter).mean.betaGrid(:,:,8) .* sin(2*pi*4*JulDayRatio)' + ARGO.(iter).mean.betaGrid(:,:,9) .* cos(2*pi*4*JulDayRatio)' ...
                         + ARGO.(iter).mean.betaGrid(:,:,10) .* sin(2*pi*5*JulDayRatio)' + ARGO.(iter).mean.betaGrid(:,:,11) .* cos(2*pi*5*JulDayRatio)' ...
                         + ARGO.(iter).mean.betaGrid(:,:,12) .* sin(2*pi*6*JulDayRatio)' + ARGO.(iter).mean.betaGrid(:,:,13) .* cos(2*pi*6*JulDayRatio)';

            ARGO.(iter).fulladt = ARGO.(iter).fulladt + S.predGrid;
            ARGO.(iter).fulladtsd = sqrt(S.predVarianceGrid);
            nanMask(:,:,1) = nanMask(:,:,1) & ~isnan(ARGO.(iter).fulladt);

            % ZonVel
            S = load(join([ARGO.(iter).anomFolder, 'anomalyDUACSRelative10_2007_2018SeasonSpaceTimeMaternlatDeriv_',num2str(mon, '%02d'),'_',num2str(year),'.mat'], ''));
            ARGO.(iter).fullzonVel = ARGO.(iter).mean.zonAbsVel + (- S.predGrid .* latDistGrid ./ gsw_f(latGrid) .* gsw_grav(latGrid)); 
            ARGO.(iter).fullzonVelsd = sqrt(S.predVarianceGrid) .* (latDistGrid ./ gsw_f(latGrid) .* gsw_grav(latGrid));
            nanMask(:,:,2) = nanMask(:,:,2) & ~isnan(ARGO.(iter).fullzonVel);

            % Meridional
            S = load(join([ARGO.(iter).anomFolder, 'anomalyDUACSRelative10_2007_2018SeasonSpaceTimeMaternlonDeriv_',num2str(mon, '%02d'),'_',num2str(year),'.mat'], ''));
            ARGO.(iter).fullmerVel = ARGO.(iter).mean.merAbsVel + (S.predGrid .* longDistGrid ./ gsw_f(latGrid) .* gsw_grav(latGrid)); 
            ARGO.(iter).fullmerVelsd = sqrt(S.predVarianceGrid) .* (longDistGrid ./ gsw_f(latGrid) .* gsw_grav(latGrid));
            nanMask(:,:,3) = nanMask(:,:,3) & ~isnan(ARGO.(iter).fullmerVel);            
        end
        nanMask = nanMask .* 1;
        nanMask(nanMask == 0) = nan;

        %% Compute MSE, llk
        % adt
        infMask = true(size(DSADT.adtInterp));
        for iiter = 1:nWindows
            iter = windowSizes(iiter);
            DUACS.adt.mse(cnt, iiter) = (mean((DSADT.adtInterp - ARGO.(iter).fulladt.*nanMask(:,:,1)) .^2, 'all', 'omitnan'));
            DUACS.adt.mad(cnt, iiter) = (median(abs(DSADT.adtInterp - ARGO.(iter).fulladt.*nanMask(:,:,1)), 'all', 'omitnan'));
            ARGO.(iter).llk = log(normpdf(DSADT.adtInterp, ARGO.(iter).fulladt.*nanMask(:,:,1), ARGO.(iter).fulladtsd));
     %       ARGO.(iter).llk(ARGO.(iter).llk < log(1e-20)) = 0;
            infMask = infMask & ~isinf(ARGO.(iter).llk);
        end
        for iiter = 1:nWindows
            iter = windowSizes(iiter);            
            DUACS.adt.llk(cnt, iiter) = sum(ARGO.(iter).llk(infMask), 'all', 'omitnan');
        end
        
        % zonal
        infMask = true(size(DS.ugosInterp));
        for iiter = 1:nWindows
            iter = windowSizes(iiter);
            DUACS.zonVel.mse(cnt, iiter) = (mean((DS.ugosInterp - ARGO.(iter).fullzonVel.*nanMask(:,:,2)) .^2 .* eqmask, 'all', 'omitnan'));
            DUACS.zonVel.mad(cnt, iiter) = (median(abs(DS.ugosInterp - ARGO.(iter).fullzonVel.*nanMask(:,:,2)) .* eqmask, 'all', 'omitnan'));
            ARGO.(iter).llk = log(normpdf(DS.ugosInterp, ARGO.(iter).fullzonVel.*nanMask(:,:,2).* eqmask, ARGO.(iter).fullzonVelsd.* eqmask));
   %         ARGO.(iter).llk(ARGO.(iter).llk < log(1e-20)) = 0;            
            infMask = infMask & ~isinf(ARGO.(iter).llk);
        end
        for iiter = 1:nWindows
            iter = windowSizes(iiter);
            DUACS.zonVel.llk(cnt, iiter) = sum(ARGO.(iter).llk(infMask), 'all', 'omitnan');
        end        
        
%{
        surfm(ARGO5.mean.latGrid, ARGO5.mean.longGrid, (ARGO5.fullzonVel - ARGO0.fullzonVel) .* eqmask); 
        surfm(ARGO5.mean.latGrid, ARGO5.mean.longGrid, S.predVarianceGrid.* eqmask);
%}
        % meridional
        infMask = true(size(DSADT.adtInterp));
        for iiter = 1:nWindows
            iter = windowSizes(iiter);
            DUACS.merVel.mse(cnt, iiter) = (mean((DS.vgosInterp - ARGO.(iter).fullmerVel.*nanMask(:,:,3)) .^2 .* eqmask, 'all', 'omitnan'));
            DUACS.merVel.mad(cnt, iiter) = (median(abs(DS.vgosInterp - ARGO.(iter).fullmerVel.*nanMask(:,:,3)) .* eqmask, 'all', 'omitnan'));
            ARGO.(iter).llk = log(normpdf(DS.vgosInterp, ARGO.(iter).fullmerVel.*nanMask(:,:,3).* eqmask, ARGO.(iter).fullmerVelsd.* eqmask));
  %          ARGO.(iter).llk(ARGO.(iter).llk < log(1e-20)) = 0;
            infMask = infMask & ~isinf(ARGO.(iter).llk);
        end
        for iiter = 1:nWindows
            iter = windowSizes(iiter);
            DUACS.merVel.llk(cnt, iiter) = sum(ARGO.(iter).llk(infMask), 'all', 'omitnan');
        end        
       cnt = cnt + 1;
%        surfm(ARGO5.mean.latGrid, ARGO5.mean.longGrid, temp0 - temp5 .* eqmask);

        % surfm(ARGO5.mean.latGrid, ARGO5.mean.longGrid, DS.ugosInterp .* eqmask);

        % caxis(ylimit);
        % cLims = caxis;  
        % colormap(darkb2r(cLims(1),cLims(2)));
        % colorbar;
        % load coast;
        % plotm(lat,long,'k');
        % title([num2str(year), '-', num2str(mon)])
        % pause(0.01)
    end
end


figure
plot(windows, sqrt(mean(DUACS.adt.mse)), '.-','MarkerSize', 20)
%errorbar(windows, sqrt(mean(DUACS.adt.mse)), sqrt(mean(DUACS.adt.sdse .^ 2)), '-*','MarkerSize',5)
set(gca,'XTick',windows)
xlabel('Window Size [\circ]')
title('RMSE - adt')

figure
plot(windows, -sum(DUACS.adt.llk), '.-','MarkerSize', 20)
set(gca,'XTick',windows)
xlabel('Window Size [\circ]')
title('NLL - adt')

figure
%plot(windows, mean(DUACS.adt.mse), '+-')
%errorbar(windows, sqrt(mean(DUACS.zonVel.mse)), sqrt(mean(DUACS.zonVel.sdse .^ 2)), '-*','MarkerSize',5)
%errorbar(windows, sqrt(mean(DUACS.merVel.mse)), sqrt(mean(DUACS.merVel.sdse .^ 2)), '-*','MarkerSize',5)
plot(windows, sqrt(mean(DUACS.zonVel.mse)), '.-','MarkerSize', 20)
hold on
plot(windows, sqrt(mean(DUACS.merVel.mse)), '.-','MarkerSize', 20)
hold off
legend('zonal', 'meridional')
set(gca,'XTick',windows)
xlabel('Window Size [\circ]')
title('RMSE')

figure
%plot(windows, sum(DUACS.adt.llk), '+-')
hold on
plot(windows, -sum(DUACS.zonVel.llk), '.-','MarkerSize', 20)
plot(windows, -sum(DUACS.merVel.llk), '.-','MarkerSize', 20)
hold off
set(gca,'XTick',windows)
xlabel('Window Size [\circ]')
legend('zonal', 'meridional')
title('NLL')


%%%%%%%%%%%%%%%%
% Second Stage %
%%%%%%%%%%%%%%%%

targetVar = 'lon' % 'lon'

windows = [3, 4, 5];
windowSizes = ["w3", "w4", "w5"]
windowTag = ["w3_3", "w4_4", "w5_5"]
nWindows = numel(windowSizes)

eqBorder = 2%1%3
eqmask = ~(latGrid < eqBorder & latGrid > - eqBorder) .* 1;
eqmask(eqmask == 0) = NaN;

% 0. Load estimates
ARGO = cell2struct(cell(nWindows, 1), cellstr(windowSizes));
for iiter = 1:nWindows
    win = windowSizes(iiter);
    ARGO.(win).mean = load(join(['./Results/meanFieldESA',targetVar,'heatFluxNoTrendPchipPotTempRelative10_2007_2018spherical_',win,'_20_Eq',num2str(eqBorder),'EM5.mat'], ""));
    ARGO.(win).anomFolder = join(['./Results/anomaly_ESA',targetVar,'heatFluxspherical_',windowTag(iiter),'_Eq',num2str(eqBorder),'_Matern_Relative10Season_11EM5/'], "");
    ARGO.(win).timeMean = ARGO.(win).mean.betaGrid(:,:,1);
end


DUACSFolder = ['./Data/DUACS_ESA_',targetVar,'Flux_ArgoRes'];
monList = 1:12%11 %1:12
nMon = length(monList)%12;
DUACS = struct();
DUACS.mse = zeros(length(2007:2018) * nMon, nWindows);
DUACS.mad = zeros(length(2007:2018) * nMon, nWindows);
DUACS.llk = zeros(length(2007:2018) * nMon, nWindows);

tempAnom = zeros([size(latGrid), nWindows]);
nanCount = zeros([size(latGrid), nWindows]);
cnt = 1;
for year = 2007:2018
    disp(num2str(year))
    for mon = monList
        %% Get DUACS OHT grid
        DS = load([DUACSFolder, '/',num2str(year, '%04d'), num2str(mon, '%02d'),'15120000.mat']);

        %% Get pseudoEstimate
        nanMask = true(size(latGrid));
        for iiter = 1:nWindows
            iter = windowSizes(iiter);

            S = load(join([ARGO.(iter).anomFolder, 'anomalyESA',targetVar,'heatFluxRelative10_2007_2018SeasonSpaceTimeMatern_',num2str(mon, '%02d'),'_',num2str(year),'.mat'], ''));
            JulDayRatio = fromJulDayToYearDay(S.midJulDay) ./ yearLength(S.midJulDay);
            ARGO.(iter).fullmean = ARGO.(iter).mean.betaGrid(:,:,1) ...
                         + ARGO.(iter).mean.betaGrid(:,:,2) .* sin(2*pi*1*JulDayRatio)' + ARGO.(iter).mean.betaGrid(:,:,3) .* cos(2*pi*1*JulDayRatio)' ...
                         + ARGO.(iter).mean.betaGrid(:,:,4) .* sin(2*pi*2*JulDayRatio)' + ARGO.(iter).mean.betaGrid(:,:,5) .* cos(2*pi*2*JulDayRatio)' ...
                         + ARGO.(iter).mean.betaGrid(:,:,6) .* sin(2*pi*3*JulDayRatio)' + ARGO.(iter).mean.betaGrid(:,:,7) .* cos(2*pi*3*JulDayRatio)' ...
                         + ARGO.(iter).mean.betaGrid(:,:,8) .* sin(2*pi*4*JulDayRatio)' + ARGO.(iter).mean.betaGrid(:,:,9) .* cos(2*pi*4*JulDayRatio)' ...
                         + ARGO.(iter).mean.betaGrid(:,:,10) .* sin(2*pi*5*JulDayRatio)' + ARGO.(iter).mean.betaGrid(:,:,11) .* cos(2*pi*5*JulDayRatio)' ...
                         + ARGO.(iter).mean.betaGrid(:,:,12) .* sin(2*pi*6*JulDayRatio)' + ARGO.(iter).mean.betaGrid(:,:,13) .* cos(2*pi*6*JulDayRatio)';
            S.predGrid(isnan(S.predGrid)) = 0;
            ARGO.(iter).fullmean = ARGO.(iter).fullmean + S.predGrid;
            nanMask = nanMask & ~isnan(ARGO.(iter).fullmean);
            
            nanCount(:,:,iiter) = nanCount(:,:,iiter) + isnan(S.predGrid);
            tempAnom(:,:,iiter) = tempAnom(:,:,iiter) + S.predGrid;
            
            ARGO.(iter).fullVar = S.predVarianceGrid;
            resDat = load(join(['./Data/','ESA',targetVar, 'heatFluxResRelative10_2007_2018Filtered_20spherical','_',iter,'_Eq',num2str(eqBorder),'EM5.mat'], ''));
            S = load(join([ARGO.(iter).anomFolder, 'anomalyESA',targetVar,'heatFluxRelative10_2007_2018SeasonSpaceTimeMaternLinCov_',num2str(mon, '%02d'),'_',num2str(year),'.mat'], ''));
            ARGO.(iter).fullVar = ARGO.(iter).fullVar + (S.predVarianceGrid ./ resDat.stdRes .^ 2);
        end
        
        nanMask = nanMask .* 1;
        nanMask(nanMask == 0) = nan;
        nanMask = nanMask .* eqmask;

        %% Compute MSE, llk
        infMask = true(size(DS.FluxInterp));
        for iiter = 1:nWindows
            iter = windowSizes(iiter);
            
%{
figure

handle = worldmap('World');
setm(handle, 'Origin', [0 200 0]);
tightmap;
mlabel('off');
plabel('off');

set(gcf, 'Position',  [100, 100, 1600, 900])

pos = get(gca,'position');
set(gca,'position',pos + [-0.13 -0.025 0.2 0])

surfm(latGrid,longGrid, ARGO.(iter).fullmean .* nanMask); 
ylimit = quantile(ARGO.(iter).fullmean(:) .* nanMask(:), [0.1 0.90]); % raw version

cLims = caxis;  
colormap(darkb2r(cLims(1),cLims(2)));
colorbar;
load coast;
plotm(lat,long,'k');
title(iter)            

figure

handle = worldmap('World');
setm(handle, 'Origin', [0 200 0]);
tightmap;
mlabel('off');
plabel('off');

set(gcf, 'Position',  [100, 100, 1600, 900])

pos = get(gca,'position');
set(gca,'position',pos + [-0.13 -0.025 0.2 0])

surfm(latGrid,longGrid, DS.FluxInterp .* nanMask); 
ylimit = quantile(DS.FluxInterp(:) .* nanMask(:), [0.1 0.90]); % raw version

cLims = caxis;  
colormap(darkb2r(cLims(1),cLims(2)));
colorbar;
load coast;
plotm(lat,long,'k');
title('DS.FluxInterp')            
                    %}    
            DUACS.mse(cnt, iiter) = (mean((DS.FluxInterp - ARGO.(iter).fullmean.*nanMask) .^2, 'all', 'omitnan'));
            DUACS.mad(cnt, iiter) = (median(abs(DS.FluxInterp - ARGO.(iter).fullmean.*nanMask), 'all', 'omitnan'));
            ARGO.(iter).llk = log(normpdf(DS.FluxInterp, ARGO.(iter).fullmean.*nanMask, sqrt(ARGO.(iter).fullVar)));
%            ARGO.(iter).llk(ARGO.(iter).llk < log(1e-19)) = 0;
            infMask = infMask & ~isinf(ARGO.(iter).llk);
        end
        for iiter = 1:nWindows
            iter = windowSizes(iiter);            
            DUACS.llk(cnt, iiter) = sum(ARGO.(iter).llk(infMask), 'all', 'omitnan');
        end
        cnt = cnt + 1;        
    end
end

figure
plot(windows, sqrt(mean(DUACS.mse)), '.-','MarkerSize', 20)
%errorbar(windows, sqrt(mean(DUACS.adt.mse)), sqrt(mean(DUACS.adt.sdse .^ 2)), '-*','MarkerSize',5)
set(gca,'XTick',windows)
xlabel('Window Size [\circ]')
title('RMSE - OHT')


figure
plot(windows, -sum(DUACS.llk), '.-','MarkerSize', 20)
%errorbar(windows, sqrt(mean(DUACS.adt.mse)), sqrt(mean(DUACS.adt.sdse .^ 2)), '-*','MarkerSize',5)
set(gca,'XTick',windows)
xlabel('Window Size [\circ]')
title('PNLL - OHT')

%% Combined Figure
%Sfin = struct();
Sfin.(targetVar) = struct();
Sfin.(targetVar).mse = DUACS.mse;
Sfin.(targetVar).llk = DUACS.llk;


figure
plot(windows, sqrt(mean(Sfin.lat.mse)), '.-','MarkerSize', 20)
hold on
plot(windows, sqrt(mean(Sfin.lon.mse)), '.-','MarkerSize', 20)
set(gca,'XTick',windows)
xlabel('Window Size [\circ]')
legend('zonal', 'meridional')
title('RMSE - OHT')


figure
plot(windows, -sum(Sfin.lat.llk), '.-','MarkerSize', 20)
hold on
plot(windows, -sum(Sfin.lon.llk), '.-','MarkerSize', 20)
hold off
%errorbar(windows, sqrt(mean(DUACS.adt.mse)), sqrt(mean(DUACS.adt.sdse .^ 2)), '-*','MarkerSize',5)
set(gca,'XTick',windows)
xlabel('Window Size [\circ]')
legend('zonal', 'meridional')

title('PNLL - OHT')
%{
figure
handle = worldmap('World');
setm(handle, 'Origin', [0 200 0]);
tightmap;
mlabel('off');
plabel('off');

set(gcf, 'Position',  [100, 100, 1600, 900])

pos = get(gca,'position');
set(gca,'position',pos + [-0.13 -0.025 0.2 0])

surfm(latGrid,longGrid, ARGO.(win).timeMean); 
ylimit = quantile(DS.FluxInterp(:) .* nanMask(:), [0.1 0.90]); % raw version

cLims = caxis;  
colormap(darkb2r(cLims(1),cLims(2)));
colorbar;
load coast;
plotm(lat,long,'k');
title('DS.FluxInterp')   
%}