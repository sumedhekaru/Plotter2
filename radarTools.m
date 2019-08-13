function radarTools(toolNum)

switch toolNum
    case 1
        radar_altitude_curser
    case 2
        real_data_2D_plot
    case 3
        real_data_3D_plot
    case 4
        smoothed_data_2D_plot
    case 5
        smoothed_data_3D_plot
    case 6
        rotate_view_180
    case 7
        show_radar_color_table
    case 8
        on_off_color_coded_locations
end


function radar_altitude_curser

sen_set = open('sensor_setting.mat');

ncid1 = netcdf.open(sen_set.radarFn,'nowrite');

% Station latitude, longitude, and elevation
stLat = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),'StationLatitude','double');
stLon = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),'StationLongitude','double');
stEle = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),'StationElevationInMeters','double')/1000;

% In horizontal coordinates
[x0, y0] = latlon2xy(stLat,-stLon);
x0 = x0/1000;
y0 = y0/1000;

eleAng = nan(1,20);

for eleAngNum =  1:20
    % Determine which number (horizontal scan number) is to load
    if eleAngNum == 1
        eleAngVar = 'elevationR_HI';
        k = 1;
    elseif eleAngNum == 2
        eleAngVar = 'elevationR_HI';
        k = 3;
    else
        eleAngVar = 'elevationR';
        k = eleAngNum-2;
    end
        
    
    % Elevation angles
    varid = netcdf.inqVarID(ncid1,eleAngVar);
    eleAngs = netcdf.getVar(ncid1,varid);

    % Elevation angles
    try
        eleAng(eleAngNum) = mean(eleAngs(:,k));
    catch
        break
    end
end

netcdf.close(ncid1);

eleAng = eleAng(1:eleAngNum -1);

h = impoint;

pos = h.getPosition;

%R = sqrt((x0 - pos(1))^2+(y0 - pos(2))^2);
%Z = R*tan(eleAng'*pi/180) + stEle;

 ke_ae = 8494.7;
 r_sq = (pos(2) - x0).^2 + (pos(2) - y0).^2;
 Z = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin(eleAng'*pi/180))-ke_ae + stEle;


fprintf('***************************\n')
fprintf('\t x0 = %6.2f\n',pos(1))
fprintf('\t y0 = %6.2f\n',pos(2))
fprintf('Indx\t Ang \t Altitude\n');
fprintf('***************************\n')

for i = 1:length(eleAng)
   fprintf('%2.2i\t\t%4.1f\t%6.2f\n',i,eleAng(i),Z(i));
end

pause(0.5)
delete(h);


function real_data_2D_plot

d0 = guidata(gcf);
vD = d0.vD;
x0 = d0.x0;
y0 = d0.y0;
cmap = d0.cmap;

figure
hold all

for k = 1:length(vD)
    
    vx = vD(k).vx;
    vy = vD(k).vy;
    vc = vD(k).vc; 
    
   
    % Parameter need to calculate the beam height (the equation for height
    % calculations is came from 
    % http://commons.wikimedia.org/wiki/File:Radar-hauteur-en.svg
    ke_ae = 8494.7;
    r_sq = (vx - x0).^2 + (vy - y0).^2;
    
    if k == 1
        % altitudes of the lowest layer
        %h1 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan(vD(k).a*pi/180);
        %h2 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k+1).a)/2*pi/180);
        
        h1 = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin(vD(k).a*pi/180))-ke_ae;
        h2 = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin((vD(k).a+vD(k+1).a)/2*pi/180))-ke_ae;        
        
    elseif k == length(vD)
        % altitudes of the top layer
        %h1 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k-1).a)/2*pi/180);
        %h2 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan(vD(k).a*pi/180);
        
        h1 = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin((vD(k).a+vD(k-1).a)/2*pi/180))-ke_ae;
        h2 = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin(vD(k).a*pi/180))-ke_ae;
    else
        % altitudes of middle layers
        % h1 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k-1).a)/2*pi/180);
        % h2 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k+1).a)/2*pi/180);        
        
        h1 = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin((vD(k).a+vD(k-1).a)/2*pi/180))-ke_ae;
        h2 = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin((vD(k).a+vD(k+1).a)/2*pi/180))-ke_ae;
    end
    
    h1 = h1 + d0.stEle;
    h2 = h2 + d0.stEle;
    
    % None smooth 2D
    r = sqrt((vx - d0.x1).^2 + (vy - d0.y1).^2);
    for i = 1:length(vx)-1
        if vc(i)+21 > 20
            pH = fill([r(i) r(i+1) r(i+1) r(i)],...
                [h1(i) h1(i+1) h2(i+1) h2(i)] ,....
                cmap(vc(i)+21,:),'edgecolor','none');
            set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        end
    end
end

% angle of the line
dx = d0.x2 - d0.x1;
dy = d0.y2 - d0.y1;

if dx > 0 && dy > 0;        ang1 = atan(dy/dx);
elseif dx > 0 && dy < 0;    ang1 = 2*pi - atan(-dy/dx);
elseif dx < 0 && dy > 0;    ang1 = pi - atan(dy/-dx);
else                        ang1 = pi + atan(dy/dx);
end


plot2Dloations(d0)
daspect([1 1 1])
%view([0, 0])
box on

xlabel([sprintf('Distance (km) from (%0.1f, %0.1f) km along ',d0.x1,d0.y1) '\theta = ' ...
    sprintf('%0.1f',ang1*180/pi) char(176)])
ylabel('Altitude (km)')

% Title
try
h = d0.plotter2Handle;
[rdn, rfn ,extn ] = fileparts(h.sen_set.radarFn);
title(sprintf('%s.%s\n%s -- %s UT', rfn,extn,sec2hhmmss(h.g.t1),sec2hhmmss(h.g.t2)),...
    'Interpreter','none')
end
set(gcf,'renderer','painters')



function plot2Dloations(d0)

lg = {};

% Do we need to plot color coded locations?
h =guidata(findall(0,'Tag','plotter2'));
cc = h.sen_set.CC_loc_on_radar;

try
    t = d0.LDAR(:,1);
    x = d0.LDAR(:,2);
    y = d0.LDAR(:,3);
    z = d0.LDAR(:,4);    
    r = convert2D(d0,x,y);
    if cc && length(t) > 1
        plot_timeCoded2D(t,r,z,d0.minLocTime,d0.maxLocTime,'o',5)
    else
        plot(r,z,'ko','markerfacecolor','k','markersize',5) 
    end
    lg = [lg 'LDAR2'];
end

try
    t = d0.CGLSS(:,1);
    x = d0.CGLSS(:,2);
    y = d0.CGLSS(:,3);
    z = d0.CGLSS(:,4);    
    r = convert2D(d0,x,y); 
    if cc && length(t) > 1
        plot_timeCoded2D(t,r,z,d0.minLocTime,d0.maxLocTime,'s',5)
    else
        plot(r,z,'ks','markerfacecolor','k','markersize',5)
    end
    lg = [lg 'CGLSS'];
end

try
    t = d0.PBFA(:,1);
    x = d0.PBFA(:,2);
    y = d0.PBFA(:,3);
    z = d0.PBFA(:,4);
    r = convert2D(d0,x,y);

    if cc  && length(t) > 1
        plot_timeCoded2D(t,r,z,d0.minLocTime,d0.maxLocTime,'p',5)
    else
        plot(r,z,'kp','markerfacecolor','r','markersize',5)
    end
    lg = [lg 'PBFA'];
end

try
    t = d0.PBFAa(:,1);
    x = d0.PBFAa(:,2);
    y = d0.PBFAa(:,3);
    z = d0.PBFAa(:,4);    
    r = convert2D(d0,x,y);
    if cc && length(t) > 1
        plot_timeCoded2D(t,r,z,d0.minLocTime,d0.maxLocTime,'*',5)
    else
        plot(r,z,'ks','markerfacecolor','r','markersize',5) 
    end
    lg = [lg 'PBFA-A'];
end

try
    t = d0.PBFAo(:,1);
    x = d0.PBFAo(:,2);
    y = d0.PBFAo(:,3);
    z = d0.PBFAo(:,4);    
    r = convert2D(d0,x,y); 
    if cc && length(t) > 1
        plot_timeCoded2D(t,r,z,d0.minLocTime,d0.maxLocTime,'d',5)
    else
        plot(r,z,'kd','markerfacecolor','r','markersize',5) 
    end
    
    lg = [lg 'PBFA-O'];
end

try
    t = d0.LINET(:,1);
    x = d0.LINET(:,2);
    y = d0.LINET(:,3);
    z = d0.LINET(:,4);    
    r = convert2D(d0,x,y); 
    if cc && length(t) > 1
        plot_timeCoded2D(t,r,z,d0.minLocTime,d0.maxLocTime,'^',5)
    else
        plot3(r,z,'k^','markerfacecolor','k','markersize',5)
    end
    lg = [lg 'LINET'];
end

try
    t = d0.NLDNc(:,1);
    x = d0.NLDNc(:,2);
    y = d0.NLDNc(:,3);
    z = d0.NLDNc(:,4);    
    r = convert2D(d0,x,y);    
    if cc && length(t) > 1
        plot_timeCoded2D(t,r,z,d0.minLocTime,d0.maxLocTime,'<',5)
    else
        plot(r,z,'k<','markerfacecolor','r','markersize',5)
    end
    lg = [lg 'NLDN-C'];
end

try
    t = d0.NLDNg(:,1);
    x = d0.NLDNg(:,2);
    y = d0.NLDNg(:,3);
    z = d0.NLDNg(:,4);    
    [x,y] = convert2D(d0,x,y);    
    if cc && length(t) > 1
        plot_timeCoded2D(t,r,z,d0.minLocTime,d0.maxLocTime,'>',5)
    else
        plot(r,z,'k>','markerfacecolor','b','markersize',5)
    end
    lg = [lg 'NLDN-G'];
end

if ~isempty(lg)
    legend(lg)
end

try
    % Put the colorbar
    t1 = d0.minLocTime;
    t2 = d0.maxLocTime;
    dt = t2 - t1;
    
    if cc
        cbh = colorbar('YTick',0:0.1:1);
        
        if dt > 0.10
            str = round(((t1:dt/10:t2)-t1)'*1000)/1000;
            set(cbh,'YTickLabel',num2str(str))
            set(get(cbh,'YLabel'),'String',sprintf('Seconds after %0.6f UT',t1));
        elseif dt > 0.0001
            str = round(((t1:dt/10:t2)-t1)'*1e6)/1000;
            set(cbh,'YTickLabel',num2str(str))
            set(get(cbh,'YLabel'),'String',sprintf('mili seconds after %0.6f UT',t1));
        else
            str = round(((t1:dt/10:t2)-t1)'*1e9)/1000;
            set(cbh,'YTickLabel',num2str(str))
            set(get(cbh,'YLabel'),'String',sprintf('micro seconds after %0.6f UT',t1));
        end
    end
end



function r = convert2D(d0,x,y)

x1 = d0.x1;
x2 = d0.x2;
y1 = d0.y1;
y2 = d0.y2;

% slope and intercept of the line (need later)
m = (y2 - y1)/(x2 - x1);
b = (x2*y1 - x1*y2)/(x2 - x1);

% Snap the location coordinates to line
x = (m*y+x-m*b)/(m*m + 1);
y = (m*m*y+m*x+b)/(m*m + 1);

% Make them 2D
r = sqrt((x - x1).^2 + (y - y1).^2);


function real_data_3D_plot

d0 = guidata(gcf);
vD = d0.vD;
x0 = d0.x0;
y0 = d0.y0;
cmap = d0.cmap;

figure
hold all

for k = 1:length(vD)
    
    vx = vD(k).vx;
    vy = vD(k).vy;
    vc = vD(k).vc; 
    
    
% Parameter need to calculate the beam height (the equation for height
    % calculations is came from 
    % http://commons.wikimedia.org/wiki/File:Radar-hauteur-en.svg
    ke_ae = 8494.7;
    r_sq = (vx - x0).^2 + (vy - y0).^2;
    
    if k == 1
        % altitudes of the lowest layer
        %h1 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan(vD(k).a*pi/180);
        %h2 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k+1).a)/2*pi/180);
        
        h1 = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin(vD(k).a*pi/180))-ke_ae;
        h2 = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin((vD(k).a+vD(k+1).a)/2*pi/180))-ke_ae;        
        
    elseif k == length(vD)
        % altitudes of the top layer
        %h1 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k-1).a)/2*pi/180);
        %h2 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan(vD(k).a*pi/180);
        
        h1 = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin((vD(k).a+vD(k-1).a)/2*pi/180))-ke_ae;
        h2 = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin(vD(k).a*pi/180))-ke_ae;
    else
        % altitudes of middle layers
        % h1 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k-1).a)/2*pi/180);
        % h2 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k+1).a)/2*pi/180);        
        
        h1 = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin((vD(k).a+vD(k-1).a)/2*pi/180))-ke_ae;
        h2 = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin((vD(k).a+vD(k+1).a)/2*pi/180))-ke_ae;
    end
    
    h1 = h1 + d0.stEle;
    h2 = h2 + d0.stEle;
        
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
    
end

% Plot map
hLine = plot(d0.map.x,d0.map.y,'Color',[0.6 0.6 0.6]);
set(get(get(hLine,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off')



% angle of the line
dx = d0.x2 - d0.x1;
dy = d0.y2 - d0.y1;

if dx > 0 && dy > 0;        ang1 = atan(dy/dx);
elseif dx > 0 && dy < 0;    ang1 = 2*pi - atan(-dy/dx);
elseif dx < 0 && dy > 0;    ang1 = pi - atan(dy/-dx);
else                        ang1 = pi + atan(dy/dx);
end

% Plot the line
pH = plot([d0.x1,d0.x2],[d0.y1,d0.y2],'k','lineWidth',2);
set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')


% Plot LDAR, PBFA ,etc.
plot3Dloations(d0)
daspect([1 1 1])
view([ang1*180/pi,30])
box on

xlabel('East (km)')
ylabel('North (km)')
zlabel('Altitude (km)')


% Title
h = d0.plotter2Handle;
[rdn, rfn ,extn ] = fileparts(h.sen_set.radarFn);
title(sprintf('%s.%s\n%s -- %s UT', rfn,extn,sec2hhmmss(h.g.t1),sec2hhmmss(h.g.t2)),...
    'Interpreter','none')

% Add viewing angle to figure data 
figD.viewAng = ang1*180/pi;
guidata(gcf,figD)

% Extra menu
f = uimenu('Label','Plotter');
uimenu(f,'Label','Rotate view 180 degress','Callback','radarTools(6)');
set(gcf,'renderer','painters')



function plot3Dloations(d0)

lg = {};

% Do we need to plot color coded locations?
h =guidata(findall(0,'Tag','plotter2'));
cc = h.sen_set.CC_loc_on_radar;

try
    t = d0.LDAR(:,1);
    x = d0.LDAR(:,2);
    y = d0.LDAR(:,3);
    z = d0.LDAR(:,4);  
    if cc
        plot_timeCoded3D(t,x,y,z,d0.minLocTime,d0.maxLocTime,'o',5)
    else
        plot3(x,y,z,'ko','markerfacecolor','k','markersize',5) 
    end
    lg = [lg 'LDAR2'];
end

try
    t = d0.CGLSS(:,1);
    x = d0.CGLSS(:,2);
    y = d0.CGLSS(:,3);
    z = d0.CGLSS(:,4);
    if cc
        plot_timeCoded3D(t,x,y,z,d0.minLocTime,d0.maxLocTime,'s',5)
    else
        plot3(x,y,z,'ks','markerfacecolor','k','markersize',5)
    end
    lg = [lg 'CGLSS'];
end

try
    t = d0.PBFA(:,1);
    x = d0.PBFA(:,2);
    y = d0.PBFA(:,3);
    z = d0.PBFA(:,4);
    if cc
        plot_timeCoded3D(t,x,y,z,d0.minLocTime,d0.maxLocTime,'p',5)
    else
        plot3(x,y,z,'kp','markerfacecolor','r','markersize',5)
    end
    lg = [lg 'PBFA'];
end

try
    t = d0.PBFAa(:,1);
    x = d0.PBFAa(:,2);
    y = d0.PBFAa(:,3);
    z = d0.PBFAa(:,4);    
    if cc
        plot_timeCoded3D(t,x,y,z,d0.minLocTime,d0.maxLocTime,'*',5)
    else
        plot3(x,y,z,'k*','markerfacecolor','r','markersize',5) 
    end
    lg = [lg 'PBFA-A'];
end

try
    t = d0.PBFAo(:,1);
    x = d0.PBFAo(:,2);
    y = d0.PBFAo(:,3);
    z = d0.PBFAo(:,4); 
    if cc
        plot_timeCoded3D(t,x,y,z,d0.minLocTime,d0.maxLocTime,'d',5)
    else
        plot3(x,y,z,'kd','markerfacecolor','m','markersize',5) 
    end
    lg = [lg 'PBFA-O'];
end

try
    t = d0.LINET(:,1);
    x = d0.LINET(:,2);
    y = d0.LINET(:,3);
    z = d0.LINET(:,4);
    if cc
        plot_timeCoded3D(t,x,y,z,d0.minLocTime,d0.maxLocTime,'^',5)
    else
        plot3(x,y,z,'k^','markerfacecolor','k','markersize',5)
    end
    lg = [lg 'LINET'];
end

try
    t = d0.NLDNc(:,1);
    x = d0.NLDNc(:,2);
    y = d0.NLDNc(:,3);
    z = d0.NLDNc(:,4); 
    if cc
        plot_timeCoded3D(t,x,y,z,d0.minLocTime,d0.maxLocTime,'<',5)
    else
        plot3(x,y,z,'k<','markerfacecolor','r','markersize',5)
    end
    lg = [lg 'NLDN-C'];
end

try
    t = d0.NLDNg(:,1);
    x = d0.NLDNg(:,2);
    y = d0.NLDNg(:,3);
    z = d0.NLDNg(:,4);
    if cc
        plot_timeCoded3D(t,x,y,z,d0.minLocTime,d0.maxLocTime,'>',5)
    else
        plot3(x,y,z,'k>','markerfacecolor','b','markersize',5)
    end
    lg = [lg 'NLDN-G'];
end

if ~isempty(lg)
    legend(lg)
end

t1 = d0.minLocTime;
t2 = d0.maxLocTime;
dt = t2 - t1;

if cc
    cbh = colorbar('YTick',0:0.1:1);
    
    if dt > 0.10
        str = round(((t1:dt/10:t2)-t1)'*1000)/1000;
        set(cbh,'YTickLabel',num2str(str))
        set(get(cbh,'YLabel'),'String',sprintf('Seconds after %0.6f UT',t1));
    elseif dt > 0.0001
        str = round(((t1:dt/10:t2)-t1)'*1e6)/1000;
        set(cbh,'YTickLabel',num2str(str))
        set(get(cbh,'YLabel'),'String',sprintf('mili seconds after %0.6f UT',t1));
    else
        str = round(((t1:dt/10:t2)-t1)'*1e9)/1000;
        set(cbh,'YTickLabel',num2str(str))
        set(get(cbh,'YLabel'),'String',sprintf('micro seconds after %0.6f UT',t1));
    end
end

function smoothed_data_2D_plot

d0 = guidata(gcf);
vD = d0.vD;
x0 = d0.x0;
y0 = d0.y0;
cmap = d0.cmap;


vx = [];
vy = [];
vc = [];
vh  = [];

for k = 1:length(vD)
    vxt = (vD(k).vx(1:end-1)+vD(k).vx(2:end))/2;
    vyt = (vD(k).vy(1:end-1)+vD(k).vy(2:end))/2;
    %vht = sqrt((vxt - x0).^2 + (vyt - y0).^2)*tan(vD(k).a*pi/180)+d0.stEle;
    
    r_sq = (vxt - x0).^2 + (vyt - y0).^2;
    ke_ae  = 8494.7;
    vht = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin(vD(k).a*pi/180))-ke_ae++d0.stEle;
    
    
    vx = [vx vxt];
    vy = [vy vyt];
    vc = [vc vD(k).vc];
    vh =  [vh vht];
end

r = sqrt((vx - d0.x1).^2 + (vy - d0.y1).^2);
    
% Find out vertical resolution
dh = sqrt((vD(1).vx(1) - x0).^2 + (vD(1).vy(1) - y0).^2)*tan(vD(2).a*pi/180) - ...
    sqrt((vD(1).vx(1) - x0).^2 + (vD(1).vy(1) - y0).^2)*tan(vD(1).a*pi/180);
iRes = dh/2;
    

% Interpolate data
[rs, zs] = meshgrid(min(r(:)):iRes:max(r(:)),min(vh(:)):iRes:max(vh(:)));

% I found some nans (why?, I don't know) in r and so removed those values
indnans = find(isnan(r));
if ~isempty(indnans)
    r(indnans) = [];
    vh(indnans) = [];
    vc(indnans) = [];
end

F = TriScatteredInterp(r',vh',vc','natural');
vq = F(rs,zs);

tfg = figure(7925);
set(tfg,'visible','off')
[cd, ch2] = contourf(rs,zs,round(vq),1:180,'edgecolor','none');
    
figure
hold all

L = length(cd);

done = 0;
i = 1;


while ~done
    
    % color value
    cv = cd(1,i);
    
    % number of data points
    N = cd(2,i);
    
    if cv > 0
        % fill polygon
        %figure(100)
        pH = patch(cd(1,i+1:i+N),cd(2,i+1:i+N),cmap(cv+21,:),'edgecolor','none');
        
        %Exclude from the legend
        set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
    
    if i +N+1 < L
        i = i+N+1;
    else
        done = 1;
    end
    
end

% angle of the line
dx = d0.x2 - d0.x1;
dy = d0.y2 - d0.y1;

if dx > 0 && dy > 0;        ang1 = atan(dy/dx);
elseif dx > 0 && dy < 0;    ang1 = 2*pi - atan(-dy/dx);
elseif dx < 0 && dy > 0;    ang1 = pi - atan(dy/-dx);
else                        ang1 = pi + atan(dy/dx);
end

% Plot LDAR2, PBFA, etc
plot2Dloations(d0)
daspect([1 1 1])
box on

xlabel([sprintf('Distance (km) from (%0.1f, %0.1f) km along ',d0.x1,d0.y1) '\theta = ' ...
    sprintf('%0.1f',ang1*180/pi) char(176)])
ylabel('Altitude (km)')

% Title
try
    h = d0.plotter2Handle;
    [rdn, rfn ,extn ] = fileparts(h.sen_set.radarFn);
    title(sprintf('%s.%s\n%s -- %s UT', rfn,extn,sec2hhmmss(h.g.t1),sec2hhmmss(h.g.t2)),...
        'Interpreter','none')
end
set(gcf,'renderer','painters')

function smoothed_data_3D_plot
clc

d0 = guidata(gcf);
vD = d0.vD;
x0 = d0.x0;
y0 = d0.y0;
cmap = d0.cmap;


vx = [];
vy = [];
vc = [];
vh  = [];

for k = 1:length(vD)
    vxt = (vD(k).vx(1:end-1)+vD(k).vx(2:end))/2;
    vyt = (vD(k).vy(1:end-1)+vD(k).vy(2:end))/2;
    %vht = sqrt((vxt - x0).^2 + (vyt - y0).^2)*tan(vD(k).a*pi/180)+d0.stEle;
    
    r_sq = (vxt - x0).^2 + (vyt - y0).^2;
    ke_ae  = 8494.7;
    vht = sqrt(r_sq + ke_ae^2 + 2*sqrt(r_sq)*ke_ae*sin(vD(k).a*pi/180))-ke_ae++d0.stEle;
    
    
    vx = [vx vxt];
    vy = [vy vyt];
    vc = [vc vD(k).vc];
    vh =  [vh vht];
end

r = sqrt((vx - d0.x1).^2 + (vy - d0.y1).^2);
    
% Find out vertical resolution
dh = sqrt((vD(1).vx(1) - x0).^2 + (vD(1).vy(1) - y0).^2)*tan(vD(2).a*pi/180) - ...
    sqrt((vD(1).vx(1) - x0).^2 + (vD(1).vy(1) - y0).^2)*tan(vD(1).a*pi/180);
iRes = dh/2;
    

% Interpolate data
[rs, zs] = meshgrid(min(r(:)):iRes:max(r(:)),min(vh(:)):iRes:max(vh(:)));


F = TriScatteredInterp(r',vh',vc','natural');
vq = F(rs,zs);

tfg = figure(7925);
set(tfg,'visible','off')
[cd, ch2] = contourf(rs,zs,round(vq),1:180,'edgecolor','none');


% angle of the line
dx = d0.x2 - d0.x1;
dy = d0.y2 - d0.y1;

if dx > 0 && dy > 0;        ang1 = atan(dy/dx);
elseif dx > 0 && dy < 0;    ang1 = 2*pi - atan(-dy/dx);
elseif dx < 0 && dy > 0;    ang1 = pi - atan(dy/-dx);
else                        ang1 = pi + atan(dy/dx);
end
    

figure
hold all

L = length(cd);

done = 0;
i = 1;

while ~done
    
    % color value
    cv = cd(1,i);
    
    % number of data points
    N = cd(2,i);
    
    if cv > 0
        % fill polygon
        xs = d0.x1 + cd(1,i+1:i+N)*cos(ang1);
        ys = d0.y1 + cd(1,i+1:i+N)*sin(ang1);
        zs = cd(2,i+1:i+N);
        
        pH = patch(xs,ys,zs,cmap(cv,:),'edgecolor','none');
        
        %Exclude from the legend
        set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
    
    if i +N+1 < L
        i = i+N+1;
    else
        done = 1;
    end
    
end

% Plot map
hLine = plot(d0.map.x,d0.map.y,'Color',[0.6 0.6 0.6]);
set(get(get(hLine,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off')

pH = plot([d0.x1,d0.x2],[d0.y1,d0.y2],'k','lineWidth',2);
set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

% Plot LDAR, PBFA ,etc.
plot3Dloations(d0)
daspect([1 1 1])
view([ang1*180/pi,30])
box on

xlabel('East (km)')
ylabel('North (km)')
zlabel('Altitude (km)')

% Title
h = d0.plotter2Handle;
[rdn, rfn ,extn ] = fileparts(h.sen_set.radarFn);
title(sprintf('%s.%s\n%s -- %s UT', rfn,extn,sec2hhmmss(h.g.t1),sec2hhmmss(h.g.t2)),...
    'Interpreter','none')

figD.viewAng = ang1*180/pi;
guidata(gcf,figD)

% Extra menu
f = uimenu('Label','Plotter');
uimenu(f,'Label','Rotate view 180 degress','Callback','radarTools(6)');

set(gcf,'renderer','painters')

function rotate_view_180

d0 = guidata(gcf);

ang = d0.viewAng;

if ang < 180
    newAng = ang + 180;
    view([newAng,30])
else
    newAng = ang - 180;
    view(newAng,30)
end

figD.viewAng = newAng;

guidata(gcf,figD)

function plot_timeCoded2D(t,r,z,minLocTime,maxLocTime,marker,mz)

cmap = colormap('jet');
L = length(cmap);

% color index of 

cind = round((t - minLocTime)/(maxLocTime - minLocTime)*(L-1)+1);

for i = 1:L
    
    inds = find(cind == i);
    if ~isempty(inds)
        pH = plot(r(inds),z(inds),marker,'color',cmap(i,:),...
            'markersize',mz,'markerfacecolor',cmap(i,:));
        set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
end

set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','on')


function show_radar_color_table
clc

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
    
    figure
    colormap(cmap)
    cbh = colorbar;
    ylabel(cbh,'dBZ')
    set(cbh,'YTick',1:20:181,'YTickLabel',-10:10:80)
        
    cbh = colorbar('Location','SouthOutside');
    xlabel(cbh,'dBZ')
    set(cbh,'XTick',1:20:181,'XTickLabel',-10:10:80)
    
    tools2fig
    

    
function on_off_color_coded_locations

gcbo = findall(findall(gcf,'type','uimenu'),'label','Color coded locations');

h=guidata(findall(0,'Tag','plotter2'));
sen_set = h.sen_set;

if strcmp(get(gcbo, 'Checked'),'on')
    set(gcbo, 'Checked', 'off');
    sen_set.CC_loc_on_radar = 0;
else
    set(gcbo, 'Checked', 'on');
    sen_set.CC_loc_on_radar = 1;
end

h.sen_set = sen_set;
save('sensor_setting.mat','-Struct','sen_set');
guidata(findall(0,'Tag','plotter2'), h)


function plot_timeCoded3D(t,x,y,z,minLocTime,maxLocTime,marker,mz)

cmap = colormap('jet');
L = length(cmap);

% color index of 

cind = round((t - minLocTime)/(maxLocTime - minLocTime)*(L-1)+1);

for i = 1:L
    
    inds = find(cind == i);
    
    if ~isempty(inds)
        
        pH = plot3(x(inds),y(inds),z(inds),marker,'color',cmap(i,:),...
            'markersize',mz,'markerfacecolor',cmap(i,:));
        set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        
    end
end

set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','on')
