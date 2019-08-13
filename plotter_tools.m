function plotter_tools(opt)

% Get all childerens of the current figure
% If there is plotyy AX(3) and AX(4) will be left and right axes
% If there is only plot command AX(3) will be the handle for the left axis

AX=findall(gcf,'Type','axes');

switch opt
   case 1
      % Updating Y ticks
      update_ygrids(AX)
   case 2
      % Updating X grids
      update_xgrids(AX)
%       errordlg('This is not available Yet. Please use manual update instead!',...
%         'Axis update error','modal')
   case 3
      update_t_range(AX)
      
   case 4
        find_del_t
   case 5
        hhmmss_curser
   case 6
       ssssss_curser
    case 7
        % zoom off
        zoom xon
    case 8
        % zoom off
        zoom yon
    case 9
        zoom off
        zoom on
        %h = zoom;
        %set(h,'Motion','both','Enable','on');    
   case 10
        find_del_v()
     
    case 11
        vedio_framing()
        
    case 12
        % make equal x y axises
        lm=axis;
        ml=max(lm(2)-lm(1),lm(4)-lm(3));
        try
            zratio=(lm(6)-lm(5))/ml;
        catch
            zratio = 1;
        end
        
        daspect([1 1 zratio])
        
    case 13
        % Fill image
        daspect('auto')
        
    case 14
        % This will find the distance of two points
        find_distance
        
    case 15
        % This will just fix the title (add the begining time to it)
        fix_title
        
    case 16
        % This will save the figure
        save_figure
        
    case 17
        CGLSS_data_curser
    case 18
        PBFA_data_curser
    case 19
        LDAR2_data_curser
    case 20
        PBFAA_data_curser
    case 21
        PBFAO_data_curser
    case 22
        NLDN2_data_curser
    case 23
        LINET_data_curser   
    case 24
        PBFA_auto_calculations
    case 25
        flash_ba(AX); % Find flash before and after
    case 26
        location_data_curser
    case 27
        vertical_radar_plane
    case 28
        name_subplots_abc
    case 29
        update_xgrids(AX)
        update_ygrids(AX)
    case 30
        delete_pbfa_point
   otherwise
      errordlg('Unknown Method!',...
        'Plotter tool error','modal')
      return
end

%% #######################################################################
function update_ygrids(AX)
try
%     if length(AX)==2
%         ylimits = get(AX(2),'YLim');
%         %xinc = (xlimits(2)-xlimits(1))/5;
%         yinc = (ylimits(2)-ylimits(1))/5;
%         
%         set(AX(3),'YTick',[ylimits(1):yinc:ylimits(2)])
%     else
        ylimits = get(AX(end),'YLim');
        %xinc = (xlimits(2)-xlimits(1))/5;
        yinc = (ylimits(2)-ylimits(1))/5;
        
        set(AX(end),'YTick',[ylimits(1):yinc:ylimits(2)])
        
        ylimits = get(AX(end-1),'YLim');
        %xinc = (xlimits(2)-xlimits(1))/5;
        yinc = (ylimits(2)-ylimits(1))/5;
        
        set(AX(end-1),'YTick',ylimits(1):yinc:ylimits(2), ...
               'YTickLabel',num2str(round(ylimits(1):yinc:ylimits(2))'))
        
%     end
catch
    errordlg('Unknown error occured. May be wrong reference for a plot!',...
        'Axis update error','modal')
end

%% ######################################################################
function update_xgrids(AX)
try
    extendticklabel(AX(end),'x',12)
catch
        errordlg('Unknown error occured. May be wrong reference for a plot!',...
            'Axis update error','modal')
end

%% ######################################################################

function update_t_range(AX)

try
    h=guidata(findall(0,'Tag','plotter2'));

    %handle to t1 and t2 in plotter 
    ht1=findall(0,'Tag','t1');
    ht2=findall(0,'Tag','t2');

    xlimits = get(AX(end),'XLim');

    h.g.t1=xlimits(1);
    h.g.t2=xlimits(2);
    
    t1=sprintf('%.6f',h.g.t1);
    t2=sprintf('%.6f',h.g.t2);
    
    set(ht1,'String',t1)
    set(ht2,'String',t2)
    
    % Updating hh, mm
    hh = floor(h.g.t1/3600);
    mm = floor((h.g.t1 - hh*3600)/300)*5;
    
    hh = sprintf('%02i',hh);    
    l =get(h.hour,'String');
    value = find(strcmp(hh,l));
    
    h.g.hhn = value;
    h.g.hh  = str2double(l(value));
    set(h.hour,'Value',value)
    
    
    mm = sprintf('%02i',mm);    
    l =get(h.minute,'String');
    value = find(strcmp(mm,l));
    
    h.g.mmn = value(1);
    h.g.mm  = str2double(l(value(1)));
    set(h.minute,'Value',value(1));
 
    % Set time maximum and minimum values
    t1minimum=h.g.hh*3600+60*h.g.mm;
    t2maximum=t1minimum+300;
    
    h.g.t1min=t1minimum;
    h.g.t2max=t2maximum;

    t1min_s=sprintf('t Min = %is',t1minimum);
    t2max_s=sprintf('t Max = %is',t2maximum);
    set(h.t1min,'String',t1min_s)
    set(h.t2max,'String',t2max_s)    
    
    guidata(findall(0,'Tag','plotter2'), h)

catch
    errordlg('Unknown error occured. May be wrong reference for a plot!',...
        'Time Range Error','modal')
end


%% ####################################################################
function find_del_t()
%try
    %msg=sprintf('Select Two Points with the left Mouse clicks\n Value automatically copied to the clip board');
    %h=msgbox(msg ,'Delta t');
    %pause(2)
    %close(h)
    [x1,y1]=ginput(1);
    pause(0.25)
    
    %     h=msgbox('Select the 2nd point with the mouse','Delta t')
    %     pause(1)
    %     close(h)
    [x2,y2]=ginput(1);
    del_t=x2-x1;
    
    if abs(del_t) < 0.001       
        msg=sprintf('Delta t = %0.1f us',del_t*1e6);
    elseif abs(del_t) < 0.1
        msg=sprintf('Delta t = %0.4f ms',del_t*1e3);
    else
        msg=sprintf('Delta t = %0.7f s',del_t);
    end
      
    clipboard('copy', del_t)
    msgbox(msg,'Delta t','modal');
%catch
%     errordlg('An unknown error occured!',...
%         'Plotter tool error','modal')
% end
    

function hhmmss_curser()
    %fg=findall(gcf,'Type','figure')
    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@curser)

function ssssss_curser()
    %fg=findall(gcf,'Type','figure')
    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@curser2)
    
    
function find_del_v()
try
    %msg=sprintf('Select Two Points with the left Mouse clicks\n Value automatically copied to the clip board');
    %h=msgbox(msg ,'Delta v');
    %pause(2)
    %close(h)
    [x1,y1]=ginput(1);
    pause(0.25)
    
    %     h=msgbox('Select the 2nd point with the mouse','Delta t')
    %     pause(1)
    %     close(h)
    [x2,y2]=ginput(1);
    del_y=sprintf('%.4f',(y2-y1));
    
    msg=sprintf(['Delta y = ', del_y , ' V']);
    clipboard('copy', del_y)
    h=msgbox(msg,'Delta y','modal');
catch
    errordlg('An unknown error occured!',...
        'Plotter tool error','modal')
end

    
function vedio_framing()
    disp('Working on')
        
function find_distance()
try
%     msg=sprintf('Select Two Points with the left Mouse clicks\n Value automatically copied to the clip board');
%     h=msgbox(msg ,'Delta t');
%     pause(2)
%     close(h)
    [x1,y1]=ginput(1);
    pause(0.25)
    
    %     h=msgbox('Select the 2nd point with the mouse','Delta t')
    %     pause(1)
    %     close(h)
    [x2,y2]=ginput(1);
    del_D=sprintf('%.2f',sqrt((x2-x1)^2+(y2-y1)^2));
    
    msg=sprintf(['Distance = ', del_D , ' km']);
    clipboard('copy', del_D)
    h=msgbox(msg,'Distance','modal');
catch
    errordlg('An unknown error occured!',...
        'Plotter tool error','modal')
end
    

function fix_title

    str=get(get(gca,'title'),'string');
    str1 = str(1,:);
    try
        str2 = str(2,:);
    catch
        str2 ='';
    end
    
    n = strfind(str1,'UT:');
    
    xrange = xlim;
    time = sec2hhmmss(xrange(1));
    
    str1(1,n+4:n+11) = time(1:8);
    
    str2 = deblank(str2);
    if isempty(str2)
        str = str1;
    else
        str = sprintf('%s\n%s',str1,str2);
    end
    title(str);
    

function save_figure

h=guidata(findall(0,'Tag','plotter2'));
sen_set = h.sen_set;


% Default folder 
df = sen_set.working_dir;
% saving graph types
try
    gt = sen_set.save_plot_types;
catch
    gt = [0 0 0 0 0 0 0 0 0 0];
    %gt(1) - fig
    %gt(2) - eps
    %gt(3) - emf
    %gt(4) - jpg
    %gt(5) - png
    %gt(6) - pdf
    %gt(7) - png no background
    %gt(8-10) - for future use...
end

% Default file name
str=get(get(gca,'title'),'string');
xrange = xlim;
t1 = xrange(1);
try
    
file=sprintf('%s//%s%s%s_%5.5us_%3.3ums_',...
        df,str(1,1:4),str(1,6:7),str(1,9:10),...
        floor(t1),floor((t1-floor(t1))*1000));
catch
    file = '';
end
 
[fn,pn] = uiputfile('*','Save As..',file); 


if fn(1)~=0
    
    % Ask with file types to save
    [choise button] = settingsdlg(...
        'description','Choose plot types', ...
        'title' , 'Figure types',...
        {'Matlab figure (*.fig)'; 'fig'}, logical(gt(1)),...
        {'EPS file (*.eps)'; 'eps'}, logical(gt(2)),...
        {'Enhansed meta file (*.emf)'; 'emf'}, logical(gt(3)) ,...
        {'JPEG (*.jpg)';'jpg'},logical(gt(4)) , ...
        {'PNG (*.png)';'png'},logical(gt(5)), ...
        {'PNG without background(*.png)';'npng'},logical(gt(7)), ...
        {'PDF (*.pdf)';'pdf'},logical(gt(6)),...
        {'TIFF (*.tiff)';'tiff'},logical(gt(8)),...
        {'Resolution (dpi)';'dpi'}, sen_set.figure_dpi);
    
    gt = [choise.fig choise.eps choise.emf ...
          choise.jpg choise.png choise.pdf ...
          0 choise.tiff 0 0];
    
   
      set(gcf,'PaperPositionMode','auto')
    
      
    if strcmp(button,'ok')
        % add figure properties to the handle structure
        add_fig_data(gcf)        
        
        if choise.fig
            fname = [pn fn '.fig'];
            save_fig_helper(fname,'None')
            gt(1) = 1;
        end
        
        if choise.eps
            fname = [pn fn '.eps'];            
            save_fig_helper(fname,'-depsc2')
            gt(2) = 1;
        end
        
        if choise.emf
            fname = [pn fn '.emf'];
            %save_fig_helper(fname,'-dmeta')
            %print(gcf,fname,'-dmeta','-painters')
            print -depsc2 -painters fname
            gt(3) = 1;
        end
        
        if choise.jpg
            fname = [pn fn '.jpg'];
            save_fig_helper(fname,'-djpeg',choise.dpi)
            gt(4) = 1;
        end
        
        if choise.png
            fname = [pn fn '.png'];
            save_fig_helper(fname,'-dpng',choise.dpi)
            gt(5) = 1;
        end
        
        if choise.npng
            fname = [pn fn '.png'];
            save_fig_helper(fname,'-dpng',choise.dpi,0)
            gt(7) = 1;
        end
        
        if choise.pdf
            fname = [pn fn '.pdf'];
            save_fig_helper(fname,'-dpdf')
            gt(6) = 1;
        end
        
        if choise.tiff
            fname = [pn fn '.tiff'];
            save_fig_helper(fname,'-dtiffn',choise.dpi)
            gt(8) = 1;
        end
    end
       
    sen_set.working_dir = pn;
    sen_set.save_plot_types = gt;
    sen_set.figure_dpi = choise.dpi;
    
    h.sen_set = sen_set;
    guidata(findall(0,'Tag','plotter2'), h)
   
end

function add_fig_data(figh)

if nargin < 2
    figh = gcf;
end

% get all the axis handles in the figure
hAllAxes = findobj(figh,'type','axes');
hLeg = findobj(hAllAxes,'tag','legend');
hAxes = setdiff(hAllAxes,hLeg);

% get position of the figure
figProp.fig_position = get(figh,'position');

figProp.hAxes = hAxes;
% get locations of all the axes
for i = 1:length(hAxes);
    figProp.position(i).position = get(hAxes(i),'position');
end

% is there current data in this figure
data = guidata(figh);

% add figure properties to the figure
data.figProp = figProp;
guidata(figh,data)


function save_fig_helper(fname,opt,dpi,bg)
   clc
    if exist(fname,'file')
        button = questdlg(['The file ' fname ' exist!'],...
            'Graph save','Replace','Exit','Exit');
    else
        button = 'Replace';
    end
    
    switch button
        case 'Exit'
            return
        case 'Replace'
            
            if nargin > 3
                % turn off annoying warning
                warning('off','MATLAB:hg:ColorSpec_None')
                % png witout background
                background = get(gcf, 'color');
                % specify transparent background
                set(gcf,'color','none');
                % create output file
                set(gcf,'InvertHardCopy','off');
                print('-dpng', 'notTransparent.png',sprintf('-r%i',dpi));
                % read image data back in
                cdata = imread('notTransparent.png');
                % write it back out - setting transparency info
                imwrite(cdata, fname, 'png', 'BitDepth', 16, 'transparency', background)
                % Remove the temorory file
                delete('notTransparent.png')
                % set the background color back
                set(gcf,'color',background);
                % turn on annoying warning back
                warning('on','MATLAB:hg:ColorSpec_None')
                
            elseif nargin > 2
                disp(fname)
                print(gcf,fname,opt,sprintf('-r%i',dpi))
            else
                if strcmp(opt,'None')
                    saveas(gcf,fname)
                else
                    print(gcf,fname,opt)
                end
            end
    end
     
    
function location_data_curser 
    %fg=findall(gcf,'Type','figure');
    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@loc_data_curser)
    
function CGLSS_data_curser
    %fg=findall(gcf,'Type','figure');
    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@cglss_curser)

function  PBFA_data_curser
    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@pbfa_curser)
    
function LDAR2_data_curser
    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@ldar2_curser)
    
function  PBFAA_data_curser
    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@pbfaA_curser)

function  PBFAO_data_curser
    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@pbfaO_curser)
    
function  NLDN2_data_curser
    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@NLDN2_curser)
    
function  LINET_data_curser
    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@LINET_curser)
    
function PBFA_auto_calculations
try
    h=guidata(findall(0,'Tag','plotter2'));
    PBFA_auto5(h)
catch
    disp('Please run plotter2 first!')
end

%% Flash Before and After
function flash_ba(AX)

datacursormode off;
datacursormode on;
 
%Get the handle to the data cursor.
menu = findall(get(gcf,'Children'),'Type','uicontextmenu');
menuCallback = get(menu,'Callback');
dataCursor = menuCallback{2};


info = getCursorInfo(dataCursor);

while isempty(info)
    info = getCursorInfo(dataCursor);
    pause(0.1)
end

datacursormode off;

t = info.Position(1);
z = info.Position(2);

delete(findall(gcf,'Type','hggroup'));

% Let's see it is the PBFA location
try; h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

settings = h.sen_set;
g = h.g;


use_xyz = settings.nearest.use_xyz;
R = settings.nearest.R;
dt = settings.nearest.dt;
idt = [settings.nearest.idt1 settings.nearest.idt2]*1e-6;
    

if g.mm < 30;    ext=0;
else    ext=30;     end

if settings.ldar_tshiftOn
    x0 = settings.x(settings.ldar_tshift_sn);
    y0 = settings.y(settings.ldar_tshift_sn);
    z0 = settings.z(settings.ldar_tshift_sn);
else
    x0 = 0 ; y0 = 0; z0 = 0;
end


% Assume clicked point is PBFA
pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);


% Load PBFA data
if exist(pbfa_fn,'file')

    pbfa = pbfaExtract(pbfa_fn,t-dt,t+dt,str2double(settings.ldar_r),x0,y0,z0,0);  
    
    ind = find(pbfa(:,6)== t);
    ind2 = find(pbfa(:,5)==z);
    
    indx = intersect(ind,ind2);
    
    if ~isempty(indx)
        ind = indx(1);
        
        fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
            pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5),pbfa(ind,8),pbfa(ind,7))
        txyz = [pbfa(ind,6),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5)];
        find_flash_before_after(txyz,R,dt,idt,use_xyz)
        
        return
    end
end


% Assume clicked point is PBFA auto data
pbfa_fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);


% Load PBFA data
if exist(pbfa_fn,'file')

    pbfa = pbfaExtract(pbfa_fn,t-dt,t+dt,str2double(settings.ldar_r),x0,y0,z0,0);
    
    
    ind = find(pbfa(:,6)== t);
    ind2 = find(pbfa(:,5)==z);
    
    indx = intersect(ind,ind2);
    
    if ~isempty(indx)
        ind = indx(1);
        
        fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
            pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5),pbfa(ind,8),pbfa(ind,7))
        txyz = [pbfa(ind,6),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5)];
        find_flash_before_after(txyz,R,dt,idt,use_xyz)
        
        return
    end
end

% Assume clicked point is LDAR2 data
dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
        -datenum(str2double(g.YYYY),0,0));
    
ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);


if exist(ldar_fn,'file')
    
    % Load PBFA data  
    [CG,CAL,ldar]= ldarExtract2(ldar_fn,t-dt,t+dt,str2double(settings.ldar_r),x0,y0,z0,0,0);
    
   
    %x = pbfa(:,6)-t;
    ind = find(ldar(:,10)==t);
    ind2 = find(ldar(:,8)==z);

    indx = intersect(ind,ind2);


    if ~isempty(indx)
        ind = indx(1);
        
        fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\n',...
            ldar(ind,1),ldar(ind,6),ldar(ind,7),ldar(ind,8))
        txyz = [ldar(ind,10),ldar(ind,6),ldar(ind,7),ldar(ind,8)];
        find_flash_before_after(txyz,R,dt,idt,use_xyz)
        
        return
    end    
end


% Assume clicked point is CGLSS data
if exist(ldar_fn,'file')
    
    %x = pbfa(:,6)-t;
    ind = find(CG(:,10)==t);
    ind2 = find(CG(:,8)==z);

    indx = intersect(ind,ind2);


    if ~isempty(indx)
        ind = indx(1);
        
        fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\n',...
            CG(ind,1),CG(ind,6),CG(ind,7),CG(ind,8))
        txyz = [CG(ind,10),CG(ind,6),CG(ind,7),CG(ind,8)];
        find_flash_before_after(txyz,R,dt,idt,use_xyz)
        
        return
    end    
end

function vertical_radar_plane
  fgh = gcf;
  h = imline;
  pos = h.getPosition;
  delete(h);
  hold all
  plot(pos(1:2),pos(3:4),'-sk','LineWidth',1)
  sen_set = open('sensor_setting.mat');
  
  % get the file name from this figure and if not let's get it from the
  % settngs
  try 
      figD = guidata(fgh);
      
      if exist(figD.fileName,'file')
          fn = figD.fileName;
      else
          fn = sen_set.radarFn;
      end
  catch
      fn = sen_set.radarFn;
  end
  
  % Plot vertical radar plane
  plot_vert_radar_plane2(fn,pos(1),pos(3),pos(2),pos(4))
  
function name_subplots_abc

figH = gcf;
ax = findall(figH,'type','axes');
% Only get axis, not legend, etc...
ax = ax(~ismember(get(ax,'Tag'),{'legend','Colobar','scribeOverlay'}));
pt1 = get(ax,'Position');

xs = zeros(size(ax));
ys = xs;


for i = 1:length(ax);
   vals = pt1{i}(:);
   xs(i) = vals(1);
   ys(i) = vals(2);
end

L = sqrt(xs.^2+(1-ys).^2);

%[L,inds] = sort(L,'descend');
[L,inds] = sort(L);


for i = inds'
    xl  = get(ax(i),'xlim');
    yl  = get(ax(i),'ylim');
    axes(ax(i))
    text(xl(1)+range(xl)*0.01,yl(2)-range(yl)*0.05, sprintf('(%s)',char(97+i-1)),'fontsize',14)

end


function delete_pbfa_point

% Get location
[t, v] = ginput(1);
fh = gcf;

% find all the lines
h = findobj(fh,'Type','line');

Rs = NaN(1,length(h));
inds = Rs;

xl = range(xlim);
yl = range(ylim);

for i = 1:length(h)
    name = get(h(i),'DisplayName');
    
    if strcmp(name,'PBFA') || strcmp(name,'PBFA-A') || strcmp(name,'PBFA-O')
        % let's get x and y data
        ts = get(h(i),'xdata');
        vs = get(h(i),'ydata');
        
        data(i).ts = ts;
        data(i).vs = vs;
        
        dist = ((ts -t)/xl).^2 + ((vs - v)/yl).^2;
        [Rs(i), inds(i)] = min(dist);
        type{i} = name;
        
    end
end

% Snap the data point
[mm, i] = nanmin(Rs);

if mm > 0.001
    button = questdlg('Coudn''t pickup the point','Deleting PBFA','Try again','Cancel','Try again');
    
    if strcmp(button,'Try again')
        plotter_tools(30)
    else
        return
    end
else
    X = data(i).ts(inds(i));
    Y = data(i).vs(inds(i));
    
    hold all
    hL = plot(X,Y,'ro','markersize',10);
    
    % Massage
    button1 = questdlg(sprintf('Do you really want to delete the circled point?\nPlease close this to cancel.'), ...
        'Deelting PBFA','Local','Server','Local+Server' ,'Local+Server');
    
    disp('Deleting PBFA point ... ')
    
    delete(hL)
    [fn1, backfile1, done1 ,fn2, backfile2, done2]= delete_PBFA_from_file(X,type{i},button1);
    
    
    % Let's user give one more chane to undo
    if done1 == 1 || done2 == 1
        button = questdlg('You still have a chance to revert the delete. Do you want to revert?',...
            'Deleting PBFA','Yes','No','No');
        
        if strcmp(button,'Yes')
            % Restore files
            
            if done1 == 1
                status = movefile(backfile1,fn1);
                
                if status
                    disp('Server PBFA file restored.')
                else
                    msg('Server PBFA file cound''t restore.')
                end
            end
            
            if done2 == 1
                status = movefile(backfile2,fn2);
                
                if status
                    disp('Local PBFA file restored.')
                else
                    msg('Local PBFA file cound''t restore.')
                end
            end
            
            
        else 
            % delete backup files
            if done1 
                delete(backfile1)
            end
            
            if done2
                delete(backfile2)
            end
        end
    end   
            
    
    
end









  
 
