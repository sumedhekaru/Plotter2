function radar_plot(sen_set,g)
%% comments

% with data variable 'Reflectivity' 
%   k = 1 - produced 2.4 degree scan
%   k = 2 - 3.4 degrees
%   k = 3 - 4.3 degrees
%   k = 4 - 5.3 degrees
%   k = 5 - 6.2 degrees
%   k = 6 - 

%clc

%% User Inputs
% fn = 'C:\Users\Sumedhe\Desktop\RADAR\sumedhe_test.test';
% %fn = 'C:\Users\Sumedhe\Downloads\sumedhe_test.test';
% fgh = figure;
% eleAngNum = 1;

if nargin < 1    
    fn = 'C:\data\2011-08-05 -- 2011-08-16\data\netCDF\MEL\2011\08\14\21\KMLB20110814_213041_V03.netcdf';
    eleAngNum = 1;
    figure
    
else
    fn = sen_set.radarFn;
    eleAngNum = sen_set.radEleAngInd;    
end

%% Programm

% Determine which number is to load
if eleAngNum == 1
    dataVar = 'Reflectivity_HI';
    distVar = 'distanceR_HI';
    anglVar = 'azimuthR_HI';
    eleAngVar = 'elevationR_HI';
    k = 1;
elseif eleAngNum == 2
        dataVar = 'Reflectivity_HI';
    distVar = 'distanceR_HI';
    anglVar = 'azimuthR_HI';
    eleAngVar = 'elevationR_HI';
    k = 3;
else
    dataVar = 'Reflectivity';
    distVar = 'distanceR';
    anglVar = 'azimuthR';
    eleAngVar = 'elevationR';
    k = eleAngNum-2;
end

% Get the hanlde for netCDF file
ncid1 = netcdf.open(fn,'nowrite');


%ncinfo(fn)
%ncdisp(fn)

% Station latitude, longitude
stLat = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),'StationLatitude','double');
stLon = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),'StationLongitude','double');

% In horizontal coordinates
%[x0, y0] = latlon2xy(28.113003,80.653999);
[x0, y0] = latlon2xy(stLat,-stLon);


%% Get all the data needed
varid = netcdf.inqVarID(ncid1,dataVar);
data = netcdf.getVar(ncid1,varid,'double');
factor = netcdf.getAtt(ncid1,varid,'scale_factor','double');
offset = netcdf.getAtt(ncid1,varid,'add_offset','double');
misVals = netcdf.getAtt(ncid1,varid,'missing_value','double');

% Elevation angles 
varid = netcdf.inqVarID(ncid1,eleAngVar);
eleAngs = netcdf.getVar(ncid1,varid);
eleAng = mean(eleAngs(:,k));

% Distance from the center
varid = netcdf.inqVarID(ncid1,distVar);
R = netcdf.getVar(ncid1,varid);

% Convert angular distance to horizontal distance
R = R*cos(eleAng*pi/180);

% Angle from the north, so x = R*sin(theta)
varid = netcdf.inqVarID(ncid1,anglVar);
azR = netcdf.getVar(ncid1,varid);
ang = azR(:,k)';


% Close the file
netcdf.close(ncid1)

% figure
% hold all
% for i= 1:12
% plot(eleAng(:,i))
% mean(eleAng(:,i))
% end




% Number of radius and number of angles
L1 = length(R);
L2 = length(azR(:,1));

%

%return
% %Test data with a histogram
% figure
% d1 = reshape(data(:,:,k),L1*L2,1);
% d1(d1 == 0) = [];
% d1(d1 == 1) = [];
% 
% subplot(2,1,1)
% hist(d1,-128:128)
% 
% subplot(2,1,2)
% inds = find(d1 < 0);
% d1(inds) = d1(inds)+256;
% d1 = d1*factor+offset;
% d1 = round(d1*2);
% hist(d1,-20:159)
% %return


  
% Color map (To go with Gibson Ridge Analyist color table% "default_BR_full.pal"
% C:\Program Files (x86)\GRLevelX\GR2Analyst_2\ColorTables\default_BR_full.pal
% Color map size is 180x3
cmap = ...
    [(64:(164-64)/39:164)'      (64:(164-64)/39:164)'   (64:(164-64)/39:164)'       % -20 : 19 (40 steps)
     (164:(100-164)/19:100)'    (164:(100-164)/19:100)' (255:(192-255)/19:192)'     % 20 :39 (20 steps)
     (64:(32-64)/19:32)'        (128:(64-128)/19:64)'   (255:(128-255)/19:128)'     % 40:59 (20 steps)
     zeros(20,1)                (255:(128-255)/19:128)' zeros(20,1)                 % 60:79 (20 steps)
     zeros(20,1)+255            (255:(128-255)/19:128)' zeros(20,1)                 % 80:99 (10 steps)
     (255:(160-255)/19:160)'    zeros(20,1)             zeros(20,1)                 % 100:119 (10 steps)
     (255:(128-255)/19:128)'    zeros(20,1)             (255:(128-255)/19:128)'     % 120:149 (10 steps)
     (255:(128-255)/19:128)'    (255:(128-255)/19:128)' (255:(128-255)/19:128)'     % 150:179 (10 steps)
     ]/255;
    

% Get a copy of real data
d1 = reshape(data(:,:,k),L1*L2,1);

% Original Reflectivity data are ranged from -128 to 128
% Note that negative values in the Reflectivity data are not really negative
% but they have to be shifted +256 to apear to be at the correct intensity.
inds = find(d1 < 0);
d1(inds) = d1(inds) + 256;

% Data should be first multiplied by the "factor" and then add the "offset"
d1 = d1*factor+offset;

% Lets find the color (index) of each data point. Note that color map defined has 
% a total of 180 steps. Since the "factor" is usually 0.5, to ovoid
% rounding errors data were multiplied by 2 asssign colors. Values less
% than -20 (really small anyways) were assigned to -20 and values greater 
% than 160 (really big and will never reach) assigned as 160.
cind = round(d1*2);
%cind(cind < -20) = -20;
%cind(cind > 160) = 160;

% Each data point must define as trapizoid with 4 corners

% corner1
% Note that when angle pass from 360 to 1, averaging didn't work and we
% need to specially handle that point
% ang1 = [(ang(1)+ang(L2)) , (ang(2:L2)+ang(1:L2-1))]/2*pi/180;
% [mm, ind] = min(360-ang);
% ang1(ind+1) = (ang(ind)+360+ang(ind+1))/2*pi/180;
% ang2 = [(ang(2:L2)+ang(1:L2-1)) , (ang(1)+ang(L2))]/2*pi/180;
% ang2(ind) = (ang(ind)+360+ang(ind+1))/2*pi/180;


ang2 = [(ang(1)+ang(L2)) , (ang(2:L2)+ang(1:L2-1))]/2*pi/180;
[mm, ind] = min(360-ang);
ang2(ind+1) = (ang(ind)+360+ang(ind+1))/2*pi/180;
ang1 = [(ang(2:L2)+ang(1:L2-1)) , (ang(1)+ang(L2))]/2*pi/180;
ang1(ind) = (ang(ind)+360+ang(ind+1))/2*pi/180;

% figure
% hold all
% plot(ang1,'ro')
% plot(ang2,'ko')
%ang1 = ang1 - 2.5*pi/180;
%ang2 = ang2 - 2.5*pi/180;

% figure
% hold all
% plot(sin(ang1))
% plot(sin(ang2))
% grid on
% return

x1s = reshape([R(1); (R(2:L1)+R(1:L1-1))/2]*sin(ang1),L1*L2,1);
y1s = reshape([R(1); (R(2:L1)+R(1:L1-1))/2]*cos(ang1),L1*L2,1);

% corner2


x2s = reshape([R(1); (R(2:L1)+R(1:L1-1))/2]*sin(ang2),L1*L2,1);
y2s = reshape([R(1); (R(2:L1)+R(1:L1-1))/2]*cos(ang2),L1*L2,1);

% corner3
x3s = reshape([(R(2:L1)+R(1:L1-1))/2; R(L1)]*sin(ang1),L1*L2,1);
y3s = reshape([(R(2:L1)+R(1:L1-1))/2; R(L1)]*cos(ang1),L1*L2,1);

% corner4
x4s = reshape([(R(2:L1)+R(1:L1-1))/2; R(L1)]*sin(ang2),L1*L2,1);
y4s = reshape([(R(2:L1)+R(1:L1-1))/2; R(L1)]*cos(ang2),L1*L2,1);


% Remove missing values
d1 = reshape(data(:,:,k),L1*L2,1);

for i = 1:length(misVals)    
    inds = find(d1 == misVals(i));
    x1s(inds) = []; x2s(inds) = []; x3s(inds) = []; x4s(inds) = [];
    y1s(inds) = []; y2s(inds) = []; y3s(inds) = []; y4s(inds) = [];
    cind(inds) = [];
    d1(inds) = [];
end


% X and Y data for patch command (collect 4 corners in to trapizoids)
xData = ([x1s, x2s, x4s, x3s]+x0)/1000;
yData = ([y1s, y2s, y4s, y3s]+y0)/1000;

%figure(fgh)
% x_limit = xlim;
% y_limit = ylim;

for i = 20:159
    indx = find(cind==i); 
    % fprintf('%5.5i\t\t%5.5i\n',i,length(indx))
    pH = patch(xData(indx,1:4)',yData(indx,1:4)',cmap(i+21,:),'edgecolor','none');
    % Exclude from the legend
    set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
end

% add radara fliename to figure so that we can use it later
figD.fileName = fn;
guidata(gcf,figD)


florida_map
daspect([1 1 1])

% Some time when zooming, PBFA LDAR points disapear. This will solve that
% problem.
set(gcf,'renderer','painters')


%plot(-7.575,5.285,'go')
%plot(-1.1115,-47.3646,'ro') % RADAR origin
%xlim([-9.8 -.8])
%ylim([-1.19 7.9])
% axis([x_limit y_limit])



