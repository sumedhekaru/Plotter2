function flash_initiation_xy_plot

%% User inputs
%fn = '\\SUMEDHE-HP\Users\sumedhe\Desktop\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-details.xlsx';
% fn = '\\SUMEDHE-HP\Users\sumedhe\Desktop\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011.xlsx';
%fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-08012010.xlsx';
%fn = 'C:\data\Share\Nadeeka\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011.xlsx';
%fn = '\\SUMEDHE-HP\Users\sumedhe\Desktop\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-08012010.xlsx';
fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011.xlsx';

% sheet number
shn = 6;

% Marker size
mz = 4;

% Turn on plots
both_sources = 1;     % Plot both 1st and 2nd source together
first_source = 0;     % Only plot the 1st source
secnd_source = 0;     % Only plot the 2nd source
individ_flash = 0;    % Plot individual flash identified by flash number fln
% flash num
fln = 5;

% Plot flashes recrsively
individ_individ = 0;


% Readar file name
%radarFn = 'C:\data\2011-07-17 -- 2011-08-04\data\netCDF\2011\07\22\17\KMLB20110722_171020_V03.netcdf';
radarFn = '\\SADAQ7\data\2010-07-01 -- 2010-08-19\data\netCDF\2010\08\01\18\KMLB20100801_180453_V03.netcdf';

angNum = 8;

%% Start program
[num,txt,~] = xlsread(fn,shn);

ID = num(:,3)
t = num(:,6);
x = num(:,7)/1000;
y = num(:,8)/1000;
z = num(:,9)/1000;
type = txt(:,6);

fg = figure(342541);
set(fg,'Visible','off')
map=colormap;
miv=min(t);
mav=max(t);

clrstep = (mav-miv)/size(map,1);
rang= miv:clrstep:mav-clrstep;

ind1 = find(num(:,4) == 1.0);
sen_set = open('sensor_setting.mat');
sen_set.radarFn = radarFn;
sen_set.radEleAngInd = angNum;

% Date
YYYY = fn(end-8:end-5);
DD = fn(end-10:end-9);
MM = fn(end-12:end-11);



%% Plotting 1st and 2nd source
if both_sources
    figure; hold all
    radar_plot(sen_set)
       
    
    for k = 1:length(ind1);
        
        i = ind1(k);
        nc = sum(rang(1:end-1) <= t(i))+1;
        
        plot(x(i),y(i),'o','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
        plot(x(i+1),y(i+1),'s','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
        text(x(i),y(i),num2str(ID(i)))
        text(x(i+1),y(i+1),num2str(ID(i)))
        
        
    end
    
    
    
    
    plot_sensors(sen_set)
    
    xdd = [x; -60.1; 3.4];
    ydd = [y; -53.1; 49.3];
    xlim([min(xdd) max(xdd)]*1.3)
    ylim([min(ydd) max(ydd)]*1.3)
    
    florida_map
    box on
    daspect([ 1 1 1])
    colormap;
    
end

%% Plot only the first source
if first_source
    figure; hold all
    
    for k = 1:length(ind1);
        
        i = ind1(k);
        nc = sum(rang(1:end-1) <= t(i))+1;
        
        plot(x(i),y(i),'o','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
        text(x(i),y(i),num2str(ID(i)))
        
    end
    
    radar_plot(sen_set)
    
    plot_sensors(sen_set)
    
    xdd = [x; -60.1; 3.4];
    ydd = [y; -53.1; 49.3];
    xlim([min(xdd) max(xdd)]*1.3)
    ylim([min(ydd) max(ydd)]*1.3)
    
    florida_map
    box on
    daspect([ 1 1 1])
    colormap;
    
end

%% Plot only the second source
if secnd_source
    figure; hold all
    
    for k = 1:length(ind1);
        
        i = ind1(k);
        nc = sum(rang(1:end-1) <= t(i))+1;
        
        plot(x(i+1),y(i+1),'s','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
        text(x(i+1),y(i+1),num2str(ID(i)))
        
    end
    
    radar_plot(sen_set)
    
    plot_sensors(sen_set)
    
    xdd = [x; -60.1; 3.4];
    ydd = [y; -53.1; 49.3];
    xlim([min(xdd) max(xdd)]*1.3)
    ylim([min(ydd) max(ydd)]*1.3)
    
    florida_map
    box on
    daspect([ 1 1 1])
    colormap;
end

%% Plotting specific flash
if individ_flash
    figure; hold all;
    radar_plot(sen_set)
    
    i = find(ID == fln);
    
    plot(x(i),y(i),'ro','markerfacecolor','r','MarkerSize',mz);
    plot(x(i+1),y(i+1),'ks','markerfacecolor','k','MarkerSize',mz);
    
    plot_sensors(sen_set)
    
    xdd = [x(i); x(i+1); -60.1; 3.4];
    ydd = [y(i); y(i+1); -53.1; 49.3];
    xlim([min(xdd) max(xdd)]*1.3)
    ylim([min(ydd) max(ydd)]*1.3)
    
    florida_map
    box on
    daspect([ 1 1 1])
    colormap;
        
    title(['Flash Number ' num2str(ID(i))])
    % Type includes the header. Therefore i should go to i+1
    legend(sprintf('1-%s',type{i+1}),sprintf('2-%s',type{i+2}))
    
    
    % Vertical Radar 2D plot
    verticleRadarPlot(fn,t(i),x(i),y(i),z(i),type{i+1},...
        x(i+1),y(i+1),z(i+1),type{i+2},ID(i))
    
end

%% Plotting individual flash
if individ_individ
    
    for fln = 1:max(ID)
        figure; hold all;
        radar_plot(sen_set)
        
        i = find(ID == fln);
        
        plot(x(i),y(i),'ro','markerfacecolor','r','MarkerSize',mz);
        plot(x(i+1),y(i+1),'ks','markerfacecolor','k','MarkerSize',mz);
        
        plot_sensors(sen_set)
        
        xdd = [x(i); x(i+1); -60.1; 3.4];
        ydd = [y(i); y(i+1); -53.1; 49.3];
        xlim([min(xdd) max(xdd)]*1.3)
        ylim([min(ydd) max(ydd)]*1.3)
        
        florida_map
        box on
        daspect([ 1 1 1])
        colormap;
        
        title(['Flash Number ' num2str(ID(i))])
        % Type includes the header. Therefore i should go to i+1
        legend(sprintf('1-%s',type{i+1}),sprintf('2-%s',type{i+2}))
        
        
        % Vertical Radar 2D plot
        verticleRadarPlot(fn,t(i),x(i),y(i),z(i),type{i+1},...
            x(i+1),y(i+1),z(i+1),type{i+2},ID(i))
                
        % save png figure
        save_bf = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Radar\2010-08-01\VS\';
        save_fn = sprintf('%s/%s-%s-%s-%i-begining-VS2D.png',save_bf,YYYY,MM,DD,fln);
        print(gcf,save_fn,'-dpng')
        delete(gcf)
        delete(gcf)
        delete(gcf)
       
    end
    
end

delete(fg)

function plot_sensors(sen_set)
    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot(sen_set.x(ns)/1000,sen_set.y(ns)/1000,'kp','markerfacecolor','r','MarkerSize',4)
            text(sen_set.x(ns)/1000,sen_set.y(ns)/1000,sen_set.sen_IDs{ns})
        end
    end
    

function verticleRadarPlot(Exl_fn,t,x1,y1,z1,type1,x2,y2,z2,type2,flashNum)

% Date
YYYY = Exl_fn(end-8:end-5);
DD = Exl_fn(end-10:end-9);
MM = Exl_fn(end-12:end-11);

% Base folder according to date
date = str2num([YYYY, MM, DD]);

if date <= 20100819
    %bf = 'C:\data\2010-07-01 -- 2010-08-19\data\';
    bf = '\\SADAQ7\data\2010-07-01 -- 2010-08-19\data\';
elseif date <= 20110716
    %bf = 'C:\data\2011-07-07 -- 2011-07-16\data\';
    bf = '\\SADAQ7\data\2011-07-07 -- 2011-07-16\data\';
elseif date <=20110804
    %bf = 'C:\data\2011-07-17 -- 2011-08-04\data\';
    bf = '\\SADAQ7\data\2011-07-17 -- 2011-08-04\data\';
elseif date <= 20110816 
    %bf = 'C:\data\2011-08-05 -- 2011-08-16\data\';
    bf = '\\SADAQ7\data\2011-08-05 -- 2011-08-16\data\';
else
    disp('Date outside the file range, please check the file date')
    return
end

% Get radar files (the first file will be the closest)
[rfn rff] = getRadarFiles(bf,YYYY,MM,DD,t);

% defulat +/- distance
dr = 10;

if x2 >= x1 && y2 > y1
    ang = atan((y2-y1)/(x2-x1));  
elseif x2 > x1 && y2 <= y1
    ang = 2*pi - atan(-(y2-y1)/(x2-x1));
elseif x2 < x1 && y2 >= y1
    ang = pi - atan((y2 - y1)/-(x2 - x1));
else 
    ang = pi + atan((y2 -y1)/(x2-x1));
end

xA = x1 - dr*cos(ang);
yA = y1 - dr*sin(ang);
xB = x2 + dr*cos(ang);
yB = y2 + dr*sin(ang);
  
% figure
% hold all
% plot(x1,y1,'ro')
% plot(x2,y2,'ks')
% 
% plot([xA xB],[yA yB])
% daspect([ 1 1 1])
% return

plot_vert_radar_plane(char(rff(1)),xA,yA,xB,yB,0)

% get figure handle
fgH = gcf;
d0 = guidata(fgH);

% plot two points
plot(x1,y1,'ro','markerfacecolor','r')
plot(x2,y2,'ks','markerfacecolor','k')
legend(['1-' type1], ['2-' type2])


view(2)

% Plot soothed vertical 2D radar plot
radarTools(4)

% consider 
xx1 = d0.x1+d0.x0;
xx2 = d0.x2+d0.x0;
yy1 = d0.y1+d0.y0;
yy2 = d0.y2+d0.y0;

% r1 and r2 for two points
r1 = sqrt((x1 - xx1).^2 + (y1 - yy1).^2);
r2 = sqrt((x2 - xx1).^2 + (y2 - yy1).^2);

% plot two points
plot(r1,z1,'ro','markerfacecolor','r')
plot(r2,z2,'ks','markerfacecolor','k')
legend(['1-' type1], ['2-' type2])

%delete(fgH)

title(sprintf('Flash #: %i  -  %s',flashNum,char(rfn(1))),'Interpreter','none')





    
% find the radar file name
function [radarFiles , radarFullFiles] = getRadarFiles(bf,YYYY,MM,DD,t)

[hms, hh , mm] = sec2hhmmss(t);

% Radar file folder
rdn = sprintf('%s/netCDF/%s/%s/%s/%2.2i/',...
    bf,YYYY,MM,DD,hh);

files = dir(rdn);

L = length(files);
times = nan(1,L-2);

radarFiles = {};
radarFullFiles = {};

if L > 2  
    for i = 3:L
        fn =  files(i).name;
        times(i-2) = str2double(fn(14:15))*3600+str2double(fn(16:17))*60+str2double(fn(18:19));
        radarFiles = [radarFiles; fn];
        radarFullFiles = [radarFullFiles; [rdn fn]];
        %return
    end
end

% if within the first 10 mins of the hour, let's load files from the
% previous hour
if mm < 10
    dNum1 = datenum(str2double(YYYY),str2double(MM),str2double(DD),hh,0,0);
    
    % substract an hour
    dNum2 = dNum1 - 1/24;
    
    dvec = datevec(dNum2);
    
    % Radar file folder
    rdn = sprintf('%s/netCDF/%4.4i/%2.2i/%2.2i/%2.2i/',...
        bf,dvec(1),dvec(2),dvec(3),dvec(4));
    
    files = dir(rdn);
    
    L = length(files);
    
    % We just need the last file of this folder
    if L > 2
        fn =  files(L).name;
        times =[times, str2double(fn(14:15))*3600+str2double(fn(16:17))*60+str2double(fn(18:19))];
        radarFiles = [radarFiles; fn];
        radarFullFiles = [radarFullFiles; [rdn fn]];
    end
end

% if within the first 10 mins of the hour, let's load files from the
% previous hour
if mm > 50
    dNum1 = datenum(str2double(YYYY),str2double(MM),str2double(DD),hh,0,0);
    
    % substract an hour
    dNum2 = dNum1 + 1/24;
    
    dvec = datevec(dNum2);
    
    % Radar file folder
    rdn = sprintf('%s/netCDF/%4.4i/%2.2i/%2.2i/%2.2i/',...
        bf,dvec(1),dvec(2),dvec(3),dvec(4));
    
    files = dir(rdn);
    
    % We just need the first file of this folder
    if L > 2
        fn =  files(3).name;
        times =[times, str2double(fn(14:15))*3600+str2double(fn(16:17))*60+str2double(fn(18:19))];
        radarFiles = [radarFiles; fn];
        radarFullFiles = [radarFullFiles; [rdn fn]];
    end
end

dts = abs(times - t);
[dts I] = sort(dts);
radarFiles = radarFiles(I);
radarFullFiles = radarFullFiles(I);
    
    



