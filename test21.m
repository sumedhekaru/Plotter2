function plot_vert_radar_plane2(fn,xxx1,yyy1,xxx2,yyy2,plot_locations)
%% comments

% with data variable 'Reflectivity'
%   k = 1 - produced 2.4 degree scan
%   k = 2 - 3.4 degrees
%   k = 3 - 4.3 degrees
%   k = 4 - 5.3 degrees
%   k = 5 - 6.2 degrees
%   k = 6 -

%clc


% fn  = 'C:/data/2011-07-17 -- 2011-08-04/data//netCDF/2011/07/22/17/KMLB20110722_170114_V03.netcdf';
% xxx1 = -19.243214160979626;
% xxx2 = -15.466771556417356;
% yyy1 = -8.459584386003893;
% yyy2 =-4.683141781441623;

% fn = 'C:/data/2011-07-17 -- 2011-08-04/data//netCDF/2011/07/22/17/KMLB20110722_170114_V03.netcdf';
% xxx1 = -18;
% xxx2 = -15;
% yyy1 = -9;
% yyy2 = -6;

fn = 'C:\Users\sumedhe\Desktop\RADAR\KMLB20110722_170114_V03.netcdf';
xxx1 = -3.05;
xxx2 = -12.3;
yyy1 = -33.92;
yyy2 = -34.88;

if nargin < 6
    plot_locations = 1;
end

%% User Inputs
% fn = 'C:\Users\Sumedhe\Desktop\RADAR\sumedhe_test.test';
% %fn = 'C:\Users\Sumedhe\Downloads\sumedhe_test.test';
% fgh = figure;
% eleAngNum = 1;

% if nargin < 1
%     %fn = 'C:\Users\Sumedhe\Desktop\RADAR\KMLB20110722_162649_V03.netcdf';
%     fn = 'C:\data\2011-07-17 -- 2011-08-04\data\netCDF\2011\07\22\16\KMLB20110722_162649_V03.netcdf';
%     eleAngNum = 1;
% else
%     fn = sen_set.radarFn;
%     eleAngNum = sen_set.radEleAngInd;
% end

%% Programm

% Get the hanlde for netCDF file
ncid1 = netcdf.open(fn);
%ncinfo(fn)
%ncdisp(fn)

% Station latitude, longitude
stLat = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),'StationLatitude','double');
stLon = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),'StationLongitude','double');

% Station elevation
stEle = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),'StationElevationInMeters','double')/1000;


% In horizontal coordinates
%[x0, y0] = latlon2xy(28.113003,80.653999);
[x0, y0] = latlon2xy(stLat,-stLon);

% variable to save vertical data
 vData = repmat(struct('x',[],'y',[],'z',[],'cl',[],'vAng',[]),1,20);
 
%  xxx1 = 36.5;
%  xxx2 = 55.6;
%  yyy1 = -183.1;
%  yyy2 = -182;
%  
%  xxx1 = 45.7;
%  yyy1 = -182.6;
%  xxx2 = 46.0;
%  yyy2 = -177.3;
 
 % Map limits
 xL1 = min([xxx1 xxx2])-10;
 xL2 = max([xxx1 xxx2])+10;
 yL1 = min([yyy1 yyy2])-10;
 yL2 = max([yyy1 yyy2])+10;
 
 % Values to use later for smoothing
 xsF = [];
 ysF = [];
 zsF = [];
 clF = [];
 
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
    
temp = 0;
    
for eleAngNum = 11:11
    
    
    
    % Determine which number (horizontal scan number) is to load
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
        
    %% Get all the data needed 
    %if eleAngNum == 1 || eleAngNum == 3
    varid = netcdf.inqVarID(ncid1,dataVar);
    data = netcdf.getVar(ncid1,varid,'double');
    factor = netcdf.getAtt(ncid1,varid,'scale_factor','double');
    offset = netcdf.getAtt(ncid1,varid,'add_offset','double');
    misVals = netcdf.getAtt(ncid1,varid,'missing_value','double');
    
    % Elevation angles
    varid = netcdf.inqVarID(ncid1,eleAngVar);
    eleAngs = netcdf.getVar(ncid1,varid);
    
    % Distance from the center
    varid = netcdf.inqVarID(ncid1,distVar);
    R = netcdf.getVar(ncid1,varid);
    
    % Angle from the north, so x = R*sin(theta)
    varid = netcdf.inqVarID(ncid1,anglVar);
    azR = netcdf.getVar(ncid1,varid);
            
    %end
    
    % Elevation angles
    try
        eleAng = double(mean(eleAngs(:,k)));
    catch
        break
    end

    
    % Convert angular distance to horizontal distance
    R = double(R)*cos(eleAng*pi/180);
        
    % Azimuthal angles
    ang = double(azR(:,k)');

    
    % Number of radius and number of angles
    L1 = length(R);
    L2 = length(azR(:,1));
    
    
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
    %cind2 = cind; % get a copy of cind
    cind(cind < -20) = -20;
    cind(cind > 160) = 160;
    cind2 = cind;
    
    if eleAngNum == 11
        
        
        fgh = figure;
        %set(fgh,'visible','off')
        hold all
        
        % Each data point must define as trapizoid with 4 corners
        
        % corner1
        % Note that when angle pass from 360 to 1, averaging didn't work and we
        % need to specially handle that point      

        ang1 = [(ang(1)+ang(L2)) , (ang(2:L2)+ang(1:L2-1))]/2*pi/180;
        [mm, ind] = min(360-ang);

        if ind == L2
            nextInd = 1;
        else
            nextInd = ind+1;
        end
        
        ang1(nextInd) = (ang(ind)+360+ang(nextInd ))/2*pi/180;
        ang2 = [(ang(2:L2)+ang(1:L2-1)) , (ang(1)+ang(L2))]/2*pi/180;
        ang2(ind) = (ang(ind)+360+ang(nextInd ))/2*pi/180;
        
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
        
        % Only plot data within xL1 xL2 yL1 and yL2 because large area 
        % plotting will make choppy 3-D figure
        inds = find(xData(:,1) < xL1);
        xData(inds,:) = [];
        yData(inds,:) = [];
        cind(inds) = [];
        
        inds = find(xData(:,1) > xL2);
        xData(inds,:) = [];
        yData(inds,:) = [];
        cind(inds) = [];
        
        inds = find(yData(:,1) < yL1);
        xData(inds,:) = [];
        yData(inds,:) = [];
        cind(inds) = [];
        
        inds = find(yData(:,1) > yL2);
        xData(inds,:) = [];
        yData(inds,:) = [];
        cind(inds) = [];
        
        
        % Plot data
        for i = -10:159
            indx = find(cind==i);
            % fprintf('%5.5i\t\t%5.5i\n',i,length(indx))
            pH = patch(xData(indx,1:4)',yData(indx,1:4)',cmap(i+21,:),'edgecolor','none');
            % Exclude from the legend
            set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        end

        
        % Lets plot florida map
        sen_set = open('sensor_setting.mat');
        
        switch sen_set.map_quality
            case 1; map = load('florida_map2.mat');
            case 2; map = load('florida_map3.mat');
            case 3; map = load('florida_map4.mat');
        end
        
        inds = find(map.x < xL1);
        map.x(inds) = [];
        map.y(inds) = [];
        
        inds = find(map.x > xL2);
        map.x(inds) = [];
        map.y(inds) = [];
        
        inds = find(map.y < yL1);
        map.x(inds) = [];
        map.y(inds) = [];
        
        inds = find(map.y > yL2);
        map.x(inds) = [];
        map.y(inds) = [];
        
        
        hLine = plot(map.x,map.y,'Color',[0.6 0.6 0.6]);
        set(get(get(hLine,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off')
        
        % Save map data for later use
        figD.map = map;
        
    end
    
    
    

   %% Let's work on plotting vertical plane 
   
    % Line inquired defined by two points    
    % Remove x0 and y0 so xy is really w.r.t. to RARAR Site 
    x1 = xxx1 - x0/1000;
    x2 = xxx2 - x0/1000;
    y1 = yyy1 - y0/1000;
    y2 = yyy2 - y0/1000;
    
    %% azimuthal angles of two points (ang measure from horizontal)
    if x1 > 0 && y1 > 0;        ang1 = 90 - atan(y1/x1)*180/pi;
    elseif x1 < 0 && y1 > 0;    ang1 = 270 + atan(-y1/x1)*180/pi;
    elseif x1 < 0 && y1 < 0;    ang1 = 270- atan(y1/x1)*180/pi;
    else                        ang1 = 90 + atan(-y1/x1)*180/pi;
    end
    
    
    if x2 > 0 && y2 > 0;        ang2 = 90 - atan(y2/x2)*180/pi;
    elseif x2 < 0 && y2 > 0;    ang2 = 270 + atan(-y2/x2)*180/pi;
    elseif x2 < 0 && y2 < 0;    ang2 = 270- atan(y2/x2)*180/pi;
    else                        ang2 = 90 + atan(-y2/x2)*180/pi;
    end
   
    % Closesest angle to ang1
    [m , mind1] = min(abs(ang - ang1));
    
    % Closesest angle to ang2
    [m , mind2] = min(abs(ang - ang2));
    
    
    
    %% Radii between two points
    R1 = sqrt(x1^2+y1^2);
    R2 = sqrt(x2^2+y2^2);
    
    % Closesest Radius to R1
    [m , mRind1] = min(abs(R/1000 - R1));
    
    % Closesest Radius to R2
    [m , mRind2] = min(abs(R/1000 - R2));
    
  
    %% Line
    pH = plot([xxx1 xxx2],[yyy1 yyy2],'ro-');
    set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    
    % Closest radius for the first point
    [mm, n1] = min(abs(R/1000 - R1));
    
    % closest angle for the first point
    [mm, n2] = min(abs(ang - ang1));
    
    % Closest radius for the first point
    [mm, n3] = min(abs(R/1000 - R2));
    
    % closest angle for the first point
    [mm, n4] = min(abs(ang - ang2));
    
    if n1 > n3;   rInc = -1;
    elseif n1 < n3; rInc = 1;
    else rInc = 0;    
    end

    
    if n2 > n4;   aInc = -1;
    elseif n2 < n4; aInc = 1;
    else aInc = 0;    
    end
    
  
    if abs(n2 - n4) > L2/2
        aInc = -aInc;
    end
    

    
    vx = [];
    vy = [];
    vc = [];    
    
    done = false;
    
    % for time out
    tOutInd = 0;

    while ~done
        
        r1 = (R(n1 - 1) + R(n1))/2;
        r2 = (R(n1 + 1) + R(n1))/2;
        
        % Index can go N2 to N1 suddenly because it is a circular
        % coordinates
        if n2 == 1
            an1 = ang(L2);      an2 = ang(n2);
            an3 = ang(n2 + 1);  an4 = ang(n2);            
        elseif n2 == L2
            an1 = ang(n2 - 1);  an2 = ang(n2);
            an3 = ang(1);       an4 = ang(n2);
        else
            an1 = ang(n2 - 1);  an2 = ang(n2);
            an3 = ang(n2 + 1);  an4 = ang(n2);
        end
        
        
        % Handle when angle pasing 360 degrees
        if abs(an1-an2) > 180
            a1 = (an1 + an2 - 360)/2;
        else
            a1 = (an1 + an2)/2;
        end
        
        if abs(an3-an4) > 180
            a2 = (an3 + an4 - 360)/2;
        else
            a2 = (an3 + an4)/2;
        end
        
        
        if aInc == -1
            % Consider going theta wall 1
            lx1 = r1/1000*sin(a1*pi/180)+x0/1000;
            ly1 = r1/1000*cos(a1*pi/180)+y0/1000;
            
            lx2 = r2/1000*sin(a1*pi/180)+x0/1000;
            ly2 = r2/1000*cos(a1*pi/180)+y0/1000;
            
            [x, y] = intersection(lx1,ly1,lx2,ly2,xxx1,yyy1,xxx2,yyy2);
            
            if ~isnan(x)
                % negative theta direction
                vx = [vx x];
                vy = [vy y];
                vc = [vc cind2(sub2ind([L1,L2],n1,n2))];
                if n2 == 1
                    n2 = L2;
                else
                    n2 = n2 - 1;
                end
                % plot(x,y,'kp')
            end   
        end
        
        % Consider going theta wall 2        
        if aInc == 1
            lx1 = r1/1000*sin(a2*pi/180)+x0/1000;
            ly1 = r1/1000*cos(a2*pi/180)+y0/1000;
            
            lx2 = r2/1000*sin(a2*pi/180)+x0/1000;
            ly2 = r2/1000*cos(a2*pi/180)+y0/1000;
            
            [x,y] = intersection(lx1,ly1,lx2,ly2,xxx1,yyy1,xxx2,yyy2);
            
            if ~isnan(x)
                % positive theta direction
                vx = [vx x];
                vy = [vy y];
                vc = [vc cind2(sub2ind([L1,L2],n1,n2))];
                if n2 == L2
                    n2 = 1;
                else
                    n2 = n2 + 1;
                end
                % plot(x,y,'kp')
            end
       
        end
        
        % Consider going R wall 1
        if rInc == -1
            lx1 = r1/1000*sin(a1*pi/180)+x0/1000;
            ly1 = r1/1000*cos(a1*pi/180)+y0/1000;
            
            lx2 = r1/1000*sin(a2*pi/180)+x0/1000;
            ly2 = r1/1000*cos(a2*pi/180)+y0/1000;
            
            [x,y] = intersection(lx1,ly1,lx2,ly2,xxx1,yyy1,xxx2,yyy2);
            
             if ~isnan(x)
                 % negative R direction
                 vx = [vx x];
                 vy = [vy y];
                 vc = [vc cind2(sub2ind([L1,L2],n1,n2))];
                 n1 = n1 - 1;
                % plot(x,y,'kp')
            end
  
        end
        
        if rInc == 1
            % Consider going R wall
            lx1 = r2/1000*sin(a1*pi/180)+x0/1000;
            ly1 = r2/1000*cos(a1*pi/180)+y0/1000;
            
            lx2 = r2/1000*sin(a2*pi/180)+x0/1000;
            ly2 = r2/1000*cos(a2*pi/180)+y0/1000;
            
            [x y] = intersection(lx1,ly1,lx2,ly2,xxx1,yyy1,xxx2,yyy2);   
                      
            if ~isnan(x)
                vx = [vx x];
                vy = [vy y];
                vc = [vc cind2(sub2ind([L1,L2],n1,n2))];
                % positive R direction
                n1 = n1 + 1;
                % plot(x,y,'kp')
             end            
        end
        

        tOutInd = tOutInd + 1;
        
        % Exit if you have collected all the data for this layer 
        % or timeout (which is more than 500 points were callected, whichi
        % is highly unlikely)
        if ~isInside(x,y,xxx1,yyy1,xxx2,yyy2) || tOutInd > 500
            done = true;
        end       
        
    end
    

    % Begining and end correction
    vx = [xxx1 vx];
    vy = [yyy1 vy];
    vx(end) = xxx2;
    vy(end) = yyy2;
    
    % Saving data for vertical info
    vD(eleAngNum).vx = vx;
    vD(eleAngNum).vy = vy;
    vD(eleAngNum).vc = vc;
    vD(eleAngNum).a = eleAng;

end

%% Ploting vertical data
    
for k = 1:length(vD)
    
    vx = vD(k).vx;
    vy = vD(k).vy;
    vc = vD(k).vc; 
    
    
%     if k == 1
%         % altitudes of the lowest layer
%         h1 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan(vD(k).a*pi/180);
%         h2 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k+1).a)/2*pi/180);
%     elseif k == length(vD)
%         % altitudes of the top layer
%         h1 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k-1).a)/2*pi/180);
%         h2 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan(vD(k).a*pi/180);
%     else
%         % altitudes of middle layers
%         h1 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k-1).a)/2*pi/180);
%         h2 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k+1).a)/2*pi/180);
%     end
%     
%     h1 = h1 + stEle;
%     h2 = h2 + stEle;

     h1 = zeros(size(vx));
     h2 = h1 + 1;
        
    %None somooth 3D
    for i = 1:length(vx)-1          
        if vc(i)+21 > 20
            pH = fill3([vx(i) vx(i+1) vx(i+1) vx(i)],...
                [vy(i) vy(i+1) vy(i+1) vy(i)],...
                [h1(i) h1(i) h2(i) h2(i)] ,....
                cmap(vc(i)+21,:),'edgecolor','none');
            set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        end
    end
    
    temp = temp + 1;
end
    

% Save vertical patch data in figure handle
figD.vD = vD;

% Other important data to save
figD.x1 = xxx1;   figD.x2 = xxx2;   figD.y1 = yyy1; figD.y2 = yyy2;
figD.x0   = x0/1000;    figD.y0   = y0/1000;
figD.stEle = stEle;
figD.cmap = cmap;

% Add data to figure handle
guidata(gcf,figD)

% % Plot position data
% if plot_locations
%     plot_position_data_3D(xL1,xL2,yL1,yL2)
% end


% Finalize the figure
[~,fname,ext] = fileparts(fn);
titStr = sprintf('%s.%s',fname,ext);

daspect([1 1 1])
xlabel('East (km)')
ylabel('North (km)')
zlabel('Altitude (km)')
title(titStr,'interpreter','none')


% Setup viewing angle
% angle of the line
dx = xxx2 - xxx1;
dy = yyy2 - yyy1;

if dx > 0 && dy > 0;        angV = atan(dy/dx);
elseif dx > 0 && dy < 0;    angV = 2*pi - atan(-dy/dx);
elseif dx < 0 && dy > 0;    angV = pi - atan(dy/-dx);
else                        angV = pi + atan(dy/dx);
end

%view([angV*180/pi,30])
box on


% % add Extra Menu to figure
% f = uimenu('Label','Plotter');
% uimenu(f,'Label','Real 2D plot','Callback','radarTools(2)');
% uimenu(f,'Label','Real 3D plot','Callback','radarTools(3)');
% uimenu(f,'Label','Smoothed 2D plot','Callback','radarTools(4)');
% uimenu(f,'Label','Smoothed 3D plot','Callback','radarTools(5)');
% uimenu(f,'Label','Show Radar Color bar','Callback','radarTools(7)');
% % if sen_set.CC_loc_on_radar
% %     checked = 'on';
% % else
% %     checked = 'off';
% % end
% uimenu(f,'Label','Color coded locations','Callback','radarTools(8)','checked',checked);

% Close the file
netcdf.close(ncid1)




function [x y] = intersection(x1,y1,x2,y2,x3,y3,x4,y4)

x = det([det([x1 ,y1; x2, y2]), x1 - x2; ...
         det([x3 ,y3; x4, y4]), x3 - x4])/...
   det([x1 - x2 , y1 - y2; ...
        x3 - x4 , y3 - y4]);
    
y = det([det([x1 ,y1; x2, y2]), y1 - y2; ...
         det([x3 ,y3; x4, y4]), y3 - y4])/...
   det([x1 - x2 , y1 - y2; ...
        x3 - x4 , y3 - y4]);
  
   
if ~isInside(x,y,x1,y1,x2,y2)
    x = NaN;
    y = NaN;
end
    
 
    
    
function answer = isInside(x,y,x1,y1,x2,y2)

% find point x,y is inside the line conecting (x1,y1), (x2,y2)
answer = true;
% sort data
if x1 > x2
    temp = x1;
    x1 = x2;
    x2 = temp;
end

if y1 > y2
    temp = y1;
    y1 = y2;
    y2 = temp;
end

if x < x1 || x > x2 || y < y1 || y > y2
    answer = false;
end



function [t,x,y,z] = limitData(t,x,y,z,xL1,xL2,yL1,yL2)

    inds = find(x < xL1);
    x(inds) = [];     y(inds) = [];     z(inds) = [];     t(inds) = [];
    
    inds = find(x > xL2);
    x(inds) = [];     y(inds) = [];     z(inds) = [];     t(inds) = [];
    
    inds = find(y < yL1);
    x(inds) = [];     y(inds) = [];     z(inds) = [];     t(inds) = [];
    
    inds = find(y > yL2);
    x(inds) = [];     y(inds) = [];     z(inds) = [];     t(inds) = [];
    

function colorcoded3D(t,x,y,z,g,marker,mz) 

map = colormap;

cL = length(map);

cIndx = round((t-g.t1)/(g.t2 - g.t1)*63)+1;


for i = 1:cL
    inds = find(cIndx == i);
    
    if ~isempty(inds)
        h = plot3(x(inds),y(inds),z(inds),marker,'color',map(i,:),...
            'markerfacecolor',map(i,:),'markersize',mz);
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
end

set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','on')

function plot_position_data_3D(xL1,xL2,yL1,yL2)
%%

% Get plotter2 data
try; h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first! Without plotter2 can''t plot location data.'); return;
end

g = h.g;
sen_set = h.sen_set;

% Time shift for LDAR and Linet
if sen_set.ldar_tshiftOn==1
    sn=sen_set.ldar_tshift_sn;
    x=sen_set.x;
    y=sen_set.y;
    z=sen_set.z;
    x0=x(sn)-sen_set.x0;
    y0=y(sn)-sen_set.y0;
    z0=z(sn)-sen_set.z0;
else
    x0=0;
    y0=0;
    z0=0;
end

if g.mm < 30;    ext=0;
else    ext=30;  end

% Location data folder for local access
loc_dir = sen_set.loc_dir;

% Variable for legend
lg = {};

% Data to save on figure handle for future use figD
figD = guidata(gcf);
figD.plotter2Handle = h;

% Time range of Location data for for later color coded by time
times = [];


% Load  and plot LDAR2 and CGLSS data
if sen_set.ldarOn || sen_set.cglssOn
    
    dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
        -datenum(str2double(g.YYYY),0,0));
    
    ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);
    
    % if file exist load data
    if exist(ldar_fn,'file')
        [CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(sen_set.ldar_r),x0,y0,z0,0,0);
        
        if sen_set.ldarOn 
            [t,x,y,z] = limitData(DLS(:,10),DLS(:,6)/1000,DLS(:,7)/1000,DLS(:,8)/1000,xL1,xL2,yL1,yL2);
            if ~isempty(t) && sum(isnan(t)) ~= length(t)
                % colorcoded3D(t,x,y,z,g,'o',4)
                plot3(x,y,z,'ko','markerfacecolor','k','markersize',5)
                lg = [lg 'LDAR2'];
                figD.LDAR = [t, x, y, z];
                times = [times; t];
            end
        end
        
        if sen_set.cglssOn
            [t,x,y,z] = limitData(CG(:,10),CG(:,6)/1000,CG(:,7)/1000,CG(:,8)/1000,xL1,xL2,yL1,yL2);
            
            if ~isempty(t) && sum(isnan(t)) ~= length(t)
                %colorcoded3D(t,x,y,z,g,'o',4)
                plot3(x,y,z,'ks','markerfacecolor','k','markersize',5)
                lg = [lg 'CGLSS'];
                figD.CGLSS = [t, x, y, z];
                times = [times; t];
            end
        end
        
        
    end
end

% Load PBFA data
if sen_set.pbfaOn    
    pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
    if exist(pbfa_fn,'file')
        PBFA=pbfaExtract(pbfa_fn,g.t1,g.t2,str2double(sen_set.ldar_r),x0,y0,z0,0);
        
        [t,x,y,z] = limitData(PBFA(:,6),PBFA(:,3)/1000,PBFA(:,4)/1000,PBFA(:,5)/1000,xL1,xL2,yL1,yL2);
        
        if ~isempty(t) && sum(isnan(t)) ~= length(t)
            plot3(x,y,z,'kp','markerfacecolor','r','markersize',7)
            lg = [lg 'PBFA'];
            figD.PBFA = [t, x, y, z];
            times = [times; t];
            
        end
        
    end
end
 
% Load PBFA-auto data
if sen_set.nldnOn    
    pbfa2_fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
    if exist(pbfa2_fn,'file')
        PBFAa=pbfaExtract(pbfa2_fn,g.t1,g.t2,str2double(sen_set.ldar_r),x0,y0,z0,0); 
        
        [t,x,y,z] = limitData(PBFAa(:,6),PBFAa(:,3)/1000,PBFAa(:,4)/1000,PBFAa(:,5)/1000,xL1,xL2,yL1,yL2);
        
        if ~isempty(t) && sum(isnan(t)) ~= length(t)
            plot3(x,y,z,'ks','markerfacecolor','r','markersize',5)
            lg = [lg 'PBFA-A'];
            figD.PBFAa = [t, x, y, z];
            times = [times; t];
        end
    end
end

% Load PBFA old data
if sen_set. pbfaOOn   
    pbfa_old_fn=sprintf('%s/PBFA_old/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
    if exist(pbfa_old_fn,'file')
        PBFAO=pbfaExtract(pbfa_old_fn,g.t1,g.t2,str2double(sen_set.ldar_r),x0,y0,z0,0); 
        
        [t,x,y,z] = limitData(PBFAO(:,6),PBFAa(:,3)/1000,PBFAO(:,4)/1000,PBFAO(:,5)/1000,xL1,xL2,yL1,yL2);
        
        if ~isempty(t) && sum(isnan(t)) ~= length(t)
            plot3(x,y,z,'ko','markerfacecolor','m','markersize',5)
            lg = [lg 'PBFA-O'];
            figD.PBFAo = [t, x, y, z];
            times = [times; t];
        end
    end
end

% Load Linet data
if sen_set.linetOn
    
    linet_fn=sprintf('%s/LINET/%s/%s/%s/linet_%s%s%s_%2.2i%2.2i.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
    if exist(linet_fn,'file')
        LINET=linetExtract2(linet_fn,g.t1,g.t2,str2double(sn_set.ldar_r),x0,y0,z0);
        
        [t,x,y,z] = limitData(LINET(:,3),LINET(:,6)/1000,LINET(:,7)/1000,LINET(:,8)/1000,xL1,xL2,yL1,yL2);
        
        if ~isempty(t) && sum(isnan(t)) ~= length(t)
            plot3(x,y,z,'^k','markerfacecolor','k','markersize',5)
            lg = [lg 'LINET'];
            figD.LINET = [t, x, y, z];
            times = [times; t];
        end
    end
end


% Load NLDN data
if sen_set.nldn2On    
    nldn2_fn=sprintf('%s/NLDN2/%s/%s/NLDN2_%s%s%s.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.YYYY{:},...
        g.MM{:},g.DD{:});
    
    if exist(nldn2_fn,'file')
        [NLDNc NLDNg] = nldnExtract(a.nldn2_fn,0,0,g.t1,g.t2,str2double(sen_set.ldar_r),x0,y0,z0,1);
        
        [t,x,y,z] = limitData(NLDNc(:,8),NLDNc(:,2)/1000,NLDNc(:,3)/1000,NLDNc(:,4)/1000,xL1,xL2,yL1,yL2);
        
        if ~isempty(t) && sum(isnan(t)) ~= length(t)
            plot3(x,y,z,'<k','markerfacecolor','r','markersize',5)
            lg = [lg 'NLDN-C'];
            figD.NLDNc = [t, x, y, z];
            times = [times; t];
        end
        
        [t,x,y,z] = limitData(NLDNg(:,8),NLDNg(:,2)/1000,NLDNg(:,3)/1000,NLDNg(:,4)/1000,xL1,xL2,yL1,yL2);
        
        if ~isempty(t) && sum(isnan(t)) ~= length(t)
            plot3(x,y,z,'>k','markerfacecolor','b','markersize',5)
            lg = [lg 'NLDN-G'];
            figD.NLDNg = [t, x, y, z];
            times = [times; t];
        end
        
        
    end
end

if ~isempty(lg)
    legend(lg)
end

figD.maxLocTime = max(times);
figD.minLocTime = min(times);

% Save location data into figure handle for fugure use
guidata(gcf,figD)
