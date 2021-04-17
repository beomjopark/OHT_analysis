%% Plot Nino Index & NPGO Index

% Formatiing
myreadtable = @(filename)readtable(filename,'HeaderLines',0, ...
    'Format','%s%d%f%f','Delimiter','space','MultipleDelimsAsOne',1);
options = weboptions('ContentReader',myreadtable);

urlONI = "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt";
ONI = webread(urlONI, options);

% Convert SEAS to Month
MON = (1:12)';
convTab = addvars(ONI(1:12,'SEAS'), MON);
ONI = join(ONI, convTab, 'Keys','SEAS');
ONI.DATE = datetime(ONI.YR, ONI.MON, 15); %midDay 15
ONI = movevars(ONI, 'MON', 'AFTER', 'YR');
ONI = movevars(ONI, 'DATE', 'BEFORE', 'YR');

ONI = ONI(ONI.YR >= 2007 & ONI.YR<=2018,:);

% Figure
ulev = 0.5;
llev = -0.5;

ONI.ANOM_EL = ONI.ANOM;
ONI.ANOM_EL(round(ONI.ANOM, 1) < ulev) = NaN;
ONI.ANOM_LA = ONI.ANOM;
ONI.ANOM_LA(round(ONI.ANOM, 1) > llev) = NaN;
%{
ONI.ANOM_POS = ONI.ANOM;
ONI.ANOM_POS(ONI.ANOM < 0) = NaN;
ONI.ANOM_NEG = ONI.ANOM;
ONI.ANOM_NEG(ONI.ANOM >= 0) = NaN;
%}

figure;
fig1 = subplot(3,1,3);
plot(ONI.DATE, ONI.ANOM, 'color', 'black', 'LineWidth',2);
%bar(ONI.DATE, ONI.ANOM, 'color', 'black', 'LineWidth',2);
hold on;
plot(ONI.DATE, ONI.ANOM_EL, 'color', 'red', 'LineWidth',2);
plot(ONI.DATE, ONI.ANOM_LA, 'color', 'blue', 'LineWidth',2);
hold off;

ax1 = gca;
ax1.XTick = [ONI.DATE(ONI.MON ==1); ONI.DATE(end)+20];
datetick('x','yyyy','keepticks');
ax1.XGrid = 'on';

ax1.XAxis.MinorTick = 'on';
ax1.XAxis.MinorTickValues = ONI.DATE;

%set(gca, 'XTick', (ONI.DATE(1) : 30 : ONI.DATE(end)) );
yline(0, 'k', 'LineStyle', '-');
yline(ulev, 'r', 'LineStyle',':', 'LineWidth',1.5);
yline(llev, 'b', 'LineStyle',':', 'LineWidth',1.5);
title('Oceanic Nino Index');

%% NPSO
%{
NPGOurl = websave('NPGO', "http://www.o3d.org/npgo/npgo.php");
NPGO = readtable(NPGOurl, 'FileType', 'text',...
    'HeaderLines',28, ...
    'Delimiter','space','MultipleDelimsAsOne',1);

NPGO = NPGO(2:(end-3),1:3);
NPGO.Properties.VariableNames = {'YR', 'MON', 'NPGO'};
NPGO.YR = cellfun(@str2num, NPGO.YR);
NPGO.MON = cellfun(@str2num, NPGO.MON);
NPGO.DATE = datetime(NPGO.YR, NPGO.MON, 15); %midDay 15
NPGO.NPGO = cellfun(@str2num, NPGO.NPGO);

NPGO = NPGO(NPGO.YR >= 2007 & NPGO.YR <=2018,:);

fig2 = subplot(3,1,1);
plot(NPGO.DATE, NPGO.NPGO, 'color', 'black', 'LineWidth',2);

ax2 = gca;
ax2.XTick = [NPGO.DATE(NPGO.MON ==1); NPGO.DATE(end)+20];
ax2.XGrid = 'on';

ax2.XAxis.MinorTick = 'on';
ax2.XAxis.MinorTickValues = NPGO.DATE;
datetick('x','yyyy','keepticks');

yline(0, 'k', 'LineStyle', '-');
%yline(ulev, 'r', 'LineStyle', ':', 'LineWidth',1.5);
%yline(llev, 'b', 'LineStyle', ':', 'LineWidth',1.5);
title('North Pacific Gyro Oscilation(NPGO) Index');


set(gcf,'units','centimeters')
set(gcf,'pos',[0 0 25 30])
set(gcf,'paperunits',get(gcf,'units')) 
set(gcf,'paperpos',get(gcf,'pos'))

%}
%print('-depsc2',['./Figures/ONI_NPGO','.eps']);