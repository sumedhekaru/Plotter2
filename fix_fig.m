function fix_fig
% Fixing Figure for posters 
% Click on the graph just before run this
try
    a=open('fix_fig.mat');
catch
    a.fixtitle = false;
    a.nx=5;
    a.ny=5;
    a.blw=1;
    a.plw=1;
    a.tfs=12;
    a.xyfs=10;
    a.tickfs=10;
    a.y2max=12500;
    a.y2min=0;
    a.RYaxisOn=0;
    a.x_units = 'ms';
    a.Ldar_fs = 10;   
    a.zero_offset = 0;
end
    
 
[choise, button] = settingsdlg(...
    'title' , 'Set This Plot for Publishing...',...
    'separator' , 'Title Only',...
    {'Fix only the title?'; 'fixtitle'}, [logical(a.fixtitle) logical(a.fixtitle)],  ...
    {'Number of X grids','nx'}, a.nx, ...    
    {'X axis units (s,ms,us)','x_units'}, a.x_units, ... 
    {'Number of Y grids','ny'}, a.ny, ...  
    {'Box Line Width','blw'}, a.blw, ...
    {'Plot Line Width','plw'}, a.plw, ...
    {'Title Font Size','tfs'}, a.tfs, ...
    {'x,y Label Font Size','xyfs'},a.xyfs,...
    {'x,y Tick Label Font Size','tickfs'}, a.tickfs, ...
    {'LDAR font Size', 'Ldar_fs'}, a.Ldar_fs ,...
    {'Zero offset time','zero_offset'},a.zero_offset , ...
    {'Is there Rigt Y axis?'; 'RYaxisOn'}, logical(a.RYaxisOn), ...
    {'Right Y min','y2min'}, a.y2min, ...
    {'Right Y max','y2max'},a.y2max,...   
    'WindowWidth' , 400,...
    'ControlWidth', 230);

if strcmp(button,'ok')
    a=choise;
    save('fix_fig.mat','-struct','a')
else
    return
end

clc

% If user asking to fix only the title, let's fix it
if a.fixtitle
    
    str=get(get(gca,'title'),'string');
    str1 = str(1,:);
    str2 = str(2,:);
    
    n = strfind(str1,'UT:');
    
    xrange = xlim;
    time = sec2hhmmss(xrange(1));
    
    str1(1,n+4:n+11) = time(1:8);
    
    str2 = deblank(str2);
    str = sprintf('%s\n%s',str1,str2);
    title(str);
    
    return
end


% Fix title
tstr=get(get(gca,'title'),'string');
try
    date=tstr(1:1,1:10);
catch
    date ='';
end

date(strfind(date,'-')) = '/';
lim = xlim
Time=sec2hhmmss(lim(1));
tstring=sprintf('%s     UT: %s',date,Time);


%title(tstring,'FontWeight','bold','fontsize',a.tfs)
title('')


switch a.x_units
    case 's'
        multiplier = 1;
    case 'ms'
        multiplier = 1000;
    case 'us'
        multiplier = 1000000;
    otherwise
        multiplier = 1000;
        
end


% Is user need to offset the x-axis
if a.zero_offset ~= 0
    
    fgx = figure('visible','off');
    plot((lim - a.zero_offset)*multiplier,[1,1])
    xlim((lim - a.zero_offset)*multiplier)
    xtic = (get(gca,'xtick')/multiplier + a.zero_offset);
    xticklabel = get(gca,'xticklabel');
    delete(fgx)
    
    set(gca,'xtick',xtic,'xticklabel',xticklabel)
    xlabStr = sprintf(' after %s UT on %s',a.zero_offset,date);
else
    xlabStr = sprintf(' after %s UT on %s',Time,date);
end

% Find out x axis multiplier and fix the x label
switch a.x_units
    case 's'
        multiplier = 1;
        % Fix x label
        set(get(gca,'xlabel'),'string',['seconds' xlabStr],'fontsize',a.xyfs)
    case 'ms'
        multiplier = 1000;
        % Fix x label
        set(get(gca,'xlabel'),'string',['milli seconds' xlabStr],'fontsize',a.xyfs)
    case 'us'
        multiplier = 1000000;
        % Fix x label
        set(get(gca,'xlabel'),'string',['micro seconds' xlabStr],'fontsize',a.xyfs)
    otherwise
        multiplier = 1000;
        % Fix x label
        set(get(gca,'xlabel'),'string',['milli seconds' xlabStr] ,'fontsize',a.xyfs)
        fprintf('Can''t understand x units - Use defaults..... - ms\n')
end





% Set font size
set(gca,'FontSize',a.tickfs)

% Changing box linewidth 
box on
set(gca,'LineWidth',a.blw)


% Fix Y label
ylstr= get(get(gca,'ylabel'),'string');

set(get(gca,'ylabel'),'string',ylstr,'fontsize',a.xyfs)


% Fixing left y ticks
lim1=ylim;
set(gca,'YTick',lim1(1):(lim1(2)-lim1(1))/a.ny:lim1(2))




% Changing the right y axiz (if available)
if a.RYaxisOn==1
    h=findobj(gcf,'Type','axes');
    set(get(h(2),'ylabel'),'string','Altitude (m)','fontsize',a.xyfs)    
    set(h(2),'FontSize',a.tickfs)
    set(h(2),'YLim',[a.y2min,a.y2max],'YTick',0:a.y2max/a.ny:a.y2max,...
        'YTickLabel',0:a.y2max/a.ny:a.y2max)
end







% Find LDAR shifted sensor
try
    LDAR_shifted = tstr(2:2,end-2:end);
    
    h=annotation('textbox',[ 0.663, 0.938, 0.10,0.05]);
    str=sprintf('LDAR times are \narrival times at %s site',LDAR_shifted);
    set(h,'String',str)
    set(h,'FontSize',a.Ldar_fs)
    set(h,'FitBoxToText','on')
    set(h,'VerticalAlignment','middle')
    set(h,'HorizontalAlignment','Center')
    set(h,'LineWidth',a.blw)
   
end



% Fixing xlabels

if a.zero_offset == 0
    set(gca,'XTick',lim(1):(lim(2)-lim(1))/a.nx:lim(2))
    xticlabels=(0:(lim(2)-lim(1))/a.nx:(lim(2)-lim(1)))*multiplier;
    set(gca,'XTickLabel',xticlabels)
end

% Changing line width
set(get(gca,'Children'), 'LineWidth',a.plw)

%set(gca,'XMinorTick','on')
%set(gca,'YMinorTick','on')
