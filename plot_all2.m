%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot all2(g)
%
%   This function writtent to specially work with plotter GUI
%   The necessary input argument g contains all the controll to the
%   plotter.
%
%   Modification History
%       2010-05-13 Created by Sumedhe Karunarathne
%       2011-07-04 Renamed as plotter 2 and add ability to plot 20 sensors
%       data
%       2011-07-04 Time shift info will dispaly if it is used. Voltage info
%       will not dilplayed
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_all(g)
%clc

% creating struct to store things
a=struct;

global AX



%% Creating the directory path for all files
% Opening sensor settings
settings=open('sensor_setting.mat');
% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

%Is manual Time shift for each sensor activated
if settings.man_tshiftOn==1
    tshift=settings.t_shift;
else
    tshift=zeros(1,60);
end

%Is offset voltages should be included?
if settings.vshiftOn==1
    vshift=settings.vshift;
else
    vshift=zeros(1,60);
end

%Is manual gain should be included?
if settings.gainOn==1
    gain=settings.gain;
else
    gain=zeros(1,60)+1;
end

% Is calibration should be included
if settings.plot_calibrated == 1
    factor = g.factor;
else
    factor = zeros(1,60)+1;
end

% Generating LP and CH Legend names
lp_legend=cell(1,60);
ch_legend=cell(1,60);

% Variable to store tshift info
tshift_str =[];

% Witch sensors has time shift turn on?
tmp = g.lpgraphs + g.chgraphs;

for i=1:20;
    lp_legend{i*3-2}=[settings.sen_IDs{i} ':lp1'];
    lp_legend{i*3-1}=[settings.sen_IDs{i} ':lp2'];
    lp_legend{(i*3)}=[settings.sen_IDs{i} ':lp3'];
    
    ch_legend{i*3-2}=[settings.sen_IDs{i} ':ch1'];
    ch_legend{i*3-1}=[settings.sen_IDs{i} ':ch2'];
    ch_legend{i*3}=[settings.sen_IDs{i} ':ch3'];
    
    if sum(tmp(i*3-2:i*3))> 0 && tshift(i*3)~= 0
        tmps = sprintf(' %s:%0.6fs', settings.sen_IDs{i},tshift(i*3));
        tshift_str = [ tshift_str tmps];
    end
end

% the variable for real legend
lg={};



% is plot empty
empty_plot=true;

%% Generating file names for lp graphs and check those files are availble

% checking witch graphs are on
lp_on=g.lpgraphs;

% Save absant file extentions
a.absant_fn={};

for i=1:60
    if lp_on(i)==1
        % finding the file extention number
        ext=mod(i,3);
        if ext==0
            ext=3;
        end
        
        % Finding the stattion ID
        sid=settings.sen_IDs{ceil(i/3)};
        
        
        % finding the file name
        filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.lp%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        % Check whether the file is exists
        if exist(filename,'file')~=0
            % Check whether the file is exists
            a.lp_fn{i}=filename;
        else
            % If file is not exist don't store the file name
            a.lp_fn{i}='';
            % Absant file
            a.absant_fn=[a.absant_fn filename];
        end
    else
        a.lp_fn{i}='';
    end
end


%% Generating file names for ch (with header files) graphs and check those files are availble

% checking witch graphs are on
ch_on=g.chgraphs;

for i=1:60
    if ch_on(i)==1
        % finding the file extention number
        ext=mod(i,3);
        if ext==0
            ext=3;
        end
        
        % Finding the stattion ID
        sid=settings.sen_IDs{ceil(i/3)};
        
        % finding the file name
        filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        
        hfilename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        
        % Check whether the ch file and h file are exists
        if exist(filename,'file')==0
            a.absant_fn=[a.absant_fn filename];
        end
        
        if exist(hfilename,'file')==0
            a.absant_fn=[a.absant_fn hfilename];
        end
        
        % If file is not exist don't store the file name
        if exist(filename,'file')==0 || exist(hfilename,'file')==0
            a.ch_fn{i}='';
            a.h_fn{i}='';
        else
            a.ch_fn{i}=filename;
            a.h_fn{i}=hfilename;
        end
    else
        a.ch_fn{i}='';
        a.h_fn{i}='';
    end
end


%% Generating file name for LDAR data
if g.ldar==1
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
    
    dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
        -datenum(str2double(g.YYYY),0,0));
    
    ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);
    
    if exist(ldar_fn,'file')~=0
        % Check whether the file is exists
        a.ldar_fn=ldar_fn;
    else
        % If file is not exist don't store the file name
        a.ldar_fn='';
        % Absant file
        a.absant_fn=[a.absant_fn ldar_fn];
    end
else
    a.ldar_fn='';
end

%% Generating file name for LINET
% Did not completed yet!!!

if g.linet==1
    
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
    
    linet_fn=sprintf('%s/LINET/%s/%s/%s/linet_%s%s%s_%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
    if exist(linet_fn,'file')~=0
        % Check whether the file is exists
        a.linet_fn=linet_fn;
    else
        % If file is not exist don't store the file name
        a.linet_fn='';
        % Absant file
        a.absant_fn=[a.absant_fn linet_fn];
    end
else
    a.linet_fn='';
end


%% Generating file name for PBFA
% Did not completed yet!!!

if g.pbfa==1
    
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
    
    pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
    if exist(pbfa_fn,'file')~=0
        % Check whether the file is exists
        a.pbfa_fn=pbfa_fn;
    else
        % If file is not exist don't store the file name
        a.pbfa_fn='';
        % Absant file
        a.absant_fn=[a.absant_fn pbfa_fn];
    end
else
    a.pbfa_fn='';
end


%% Letting user know about missing file

a.absant_fn=sort(a.absant_fn);

if isempty (a.absant_fn)==0
    answer=questdlg(a.absant_fn,'Files not found!','OK','Stop!','OK');
    if strcmp(answer,'Stop!')
        return
    end
end


%% If there is nothing to plot Exit
if sum(strcmp(a.lp_fn,''))==60 && ...
        sum(strcmp(a.ch_fn,''))== 60 && ...
        strcmp(a.ldar_fn,'')== 1 && ...
        strcmp(a.linet_fn,'')==1 && ...
        strcmp(a.pbfa_fn,'')==1
    errordlg('No data files were found. All Plot commads failed!', ...
        'Plotter Error','model')
    return
end



%% Loading and plotting LP data

wbh= waitbar(0,'Please wait...','Name','Plotter Busy');

fg=figure;
set(fg,'visible','off')
hold all

% Store one data set for plot yy
ldar_t=[];
ldar_v=[];

for i=1:60
    % for loop used becuase there may be 15 lp plots
    if strcmp(a.lp_fn{i},'')==0
        [tn,vn]=SA_Extract1(a.lp_fn{i},g.t1,g.t2,tshift(i));
        ldar_t=tn;
        ldar_v=vn;
        plot(tn,vn*gain(i)*factor(i)+vshift(i))
        lg=[lg lp_legend{i}];
        empty_plot=false;
        try
            waitbar(i*0.005,wbh,'Loading slow antenna data','Name','Plotter Busy')
        catch
            return
        end
    end
end




%% Loading and plotting ch data

t=[];   % Variable for time
y=[];   % Variable for voltage
empty_trigs=[]; % Variable for empty triggers

for i=1:60
    
    try
        waitbar(0.3+i*0.005,wbh,'Loading Fast Antenna data','Name','Plotter Busy')
    catch
        return
    end
    
    t=[];   % Variable for time
    y=[];   % Variable for voltage
    % Plotting up to 60 fast antenna files
    if strcmp(a.ch_fn{i},'')==0
        % Load corresponding header file times
        fId = fopen(a.h_fn{i}, 'r');
        trigs  = fread(fId, inf, 'double') ;
        fclose( fId );
        % Introduce time shift before filtering time range
        trigs = trigs + tshift(i);
        % Finding the triggers between the given time range
        lol=length(trigs)- nnz(trigs>g.t1)+1;       % Index of the lower matrix element
        ul=nnz(trigs<g.t2);                      % Index of the upper matrix element
        % Let's find two more triggers from both ends
        if lol>1
            lol=lol-1;
        end
        
        if ul < length(trigs)
            ul=ul+1;
        end
        % Triggers between given time range and remove time shift
        trigs=trigs(lol:ul)-tshift(i);
        
        % Loading fast antenna data
        for j=1:length(trigs)
            sFa = epp_load_trigfile_time(a.ch_fn{i},trigs(j));
            %%%%%%%%%%%%%%%%% Reducing Sampling Rate %%%%%%%%%%%%%%%%%%%%%
            switch settings.chSumMethod
                case 1 % User needs downsampling
                    % find the current frequency
                    f=round(1/(sFa.t_i(2)-sFa.t_i(1)));
                    % Number of intervals
                    ni = round(f/settings.chfreq);
                    
                    if ni > 1
                        sFa.t_i = downsample(sFa.t_i,ni);
                        sFa.y_i = downsample(sFa.y_i,ni);
                    end
                    
                case 2 % User needs summing
                    % find the current frequency
                    f=round(1/(sFa.t_i(2)-sFa.t_i(1)));
                    % Number of intervals
                    ni = round(f/settings.chfreq);
                    
%                     if ni > 1
%                         len = floor(length(sFa.t_i)/ni);
%                         
%                         tempt=nan(1,len);
%                         tempy=nan(1,len);
%                         
%                         
%                         for ind = 1:len
%                             tempt(ind)=mean(sFa.t_i((ind-1)*ni+1:ind*ni));
%                             tempy(ind)=mean(sFa.y_i((ind-1)*ni+1:ind*ni));
%                         end
%                         
%                         sFa.t_i = tempt;
%                         sFa.y_i = tempy;
%                         
%                     end
                    
                    if ni > 1
                        len = floor(length(sFa.t_i)/ni);
                        
                        temp1 = nan(ni,len);
                        temp2 = temp1;
                        
                        for indn = 1:ni
                            tempt= downsample(sFa.t_i,ni,indn-1);
                            tempy= downsample(sFa.y_i,ni,indn-1);
                            temp1(indn,:)=tempt(1:len);
                            temp2(indn,:)=tempy(1:len);
                            clear tempt tempy
    
                        end
                        
                        sFa.t_i=mean(temp1);
                        sFa.y_i=mean(temp2);
                    end
                    
                otherwise
                    % User needs nothing
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Keep a track of empty trigs
            if isempty(sFa.y_i)
                empty_trigs=[empty_trigs,trigs(j)]
            else
                y=[y,NaN,sFa.y_i];
                t=[t,NaN,sFa.t_i];
            end
        end
        tshift(i);
        t=t+tshift(i);
        % Finding data in given time range
        lol=length(t)- nnz(t>g.t1)+1 ;      % Index of the lower matrix element
        ul=nnz(t<g.t2);                     % Index of the upper matrix element
        y=y(lol:ul);
        t=t(lol:ul);
        if isempty(t)==0
            ldar_t=t;
            ldar_v=y;
            plot(t,y*gain(i)*factor(i)+vshift(i))
            lg=[lg ch_legend{i}];
            empty_plot=false;
        end
    end
end

if settings.plot_calibrated == 1
    ylabel('E (V/m)')
else
    ylabel('Voltage (V)')
end

%% Time shift for LDAR and Linet
if settings.ldar_tshiftOn==1
    sn=settings.ldar_tshift_sn;
    x=settings.x;
    y=settings.y;
    z=settings.z;
    x0=x(sn)-settings.x0;
    y0=y(sn)-settings.y0;
    z0=z(sn)-settings.z0;
else
    x0=0;
    y0=0;
    z0=0;
end

%% %%%%%%%%%%%%%%   PLOTTING LINET / LDAR / PBFA %%%%%%%%%%%%%%%%%%%%%%%%%
ld_lin='';
x1=[];
y1=[];

if strcmp(a.ldar_fn,'')==0
    ld_lin=[ld_lin '-LDAR-'];
    try
        waitbar(0.7,wbh,'Loading LDAR data','Name','Plotter Busy')
    catch
        return
    end
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0,0);
    
    x1=[DLS(:,10);CG(:,10)];
    y1=[DLS(:,8);CG(:,8)];
    %     %For future use to extract ldar data in to the workspace
    %     assignin('base', 'lt', x1)
    %     assignin('base', 'lx', [DLS(:,6);CG(:,6)])
    %     assignin('base', 'ly', [DLS(:,7);CG(:,7)])
    %     assignin('base', 'lz', y1)
    
    
end

if strcmp(a.linet_fn,'')==0
    ld_lin=[ld_lin '-LINET-'];
    try
        waitbar(0.73,wbh,'Loading LINET data','Name','Plotter Busy')
    catch
        return
    end
    LINET=linetExtract(linet_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0);
    % LINET Column DISCRIPTION
    % 1 - time
    % 2 - Distance from the time correcting sensor
    % 3 - time shifted time
    % 6 - x distance
    % 7 - y distance
    % 8 - z distance
    x1=[x1; LINET(:,3)];
    y1=[y1; LINET(:,8)];
    
end

if strcmp(a.pbfa_fn,'')==0
    ld_lin=[ld_lin '-PBFA-'];
    try
        waitbar(0.76,wbh,'Loading PBFA data','Name','Plotter Busy')
    catch
        return
    end
    PBFA=pbfaExtract(pbfa_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0);
    
    % PBFA column Discription
    %   1 - distance to the pulse from the given sensor
    %   2 - Occuring time
    %   3 - x
    %   4 - y
    %   5 - z
    %   6 - Detection time at the sensor
    
    
    x1=[x1; PBFA(:,6)];
    y1=[y1; PBFA(:,5)];
    
end





if isempty(x1)==0
    try
        waitbar(0.78,wbh,'Plotting graphs','Name','Plotter Busy')
    catch
        return
    end
%     if isempty (ldar_t)==0
%         %[AX,H1,H2]=plotyy(ldar_t,ldar_v,x1,y1);
%         % commented the above line and entered the following line
%         %entered the following line to fix a auto scale problem
%         [AX,H1,H2]=plotyy(nan,nan,x1,y1);
%         empty_plot=false;
%     else
%         ldar_t=x1;
%         ldar_v=x1-x1;
%         %[AX,H1,H2]=plotyy(ldar_t,ldar_v,x1,y1);
%         [AX,H1,H2]=plotyy(nan,nan,x1,y1);
%         empty_plot=false;
%     end
    set(H1,'LineStyle','none')
    set(H2,'LineStyle','none')
    set(H2,'Marker','p','MarkerEdgeColor','g','MarkerFaceColor','g')
    
    if settings.plot_calibrated == 1
        set(get(AX(1),'Ylabel'),'String','E (V/m)')
    else
        set(get(AX(1),'Ylabel'),'String','Voltage (V)')
    end
    
    set(get(AX(2),'Ylabel'),'String','Altitude (m)')
    set(AX,'xlim',[g.t1 g.t2])
    linkaxes(AX,'x')
    set(AX(1),'ylimmode','auto')
    %mmm=ceil(max([y1;y2])/1000)*1000;
    mmm=12500;
    set(AX(2),'yLim',[0 mmm],'YTick',0:2000:mmm);
    set(AX(2),'xtick',[])
    set(AX, 'YColor', [0 0 0])
    
    
    %xlimits = get(ax1,'XLim');
    ylimits = get(AX(1),'YLim');
    %xinc = (xlimits(2)-xlimits(1))/5;
    yinc = (ylimits(2)-ylimits(1))/5;
    
    set(AX(1),'YTick',ylimits(1):yinc:ylimits(2))
    
    ylimits = get(AX(2),'YLim');
    %xinc = (xlimits(2)-xlimits(1))/5;
    yinc = (ylimits(2)-ylimits(1))/5;
    
    set(AX(2),'YTick',ylimits(1):yinc:ylimits(2))
    
    %lg=[lg 'LDAR'];
else
    AX=get(gcf,'CurrentAxes');
end

try
    waitbar(0.8,wbh,'All plotting finished!','Name','Plotter Busy')
catch
    return
end

% Create additional tool menu in the plot window

try
    waitbar(0.85,wbh,'Creating Extra Plotting Tools','Name','Plotter Busy')
catch
    return
end

f = uimenu('Label','Plotter');
uimenu(f,'Label','Update Y Grids','Callback','plotter_tools(1)','Accelerator','Y');
uimenu(f,'Label','Update X Grids','Callback','plotter_tools(2)','Accelerator','X');
uimenu(f,'Label','Update Time Range','Callback','plotter_tools(3)','Accelerator','T');
uimenu(f,'Label','Find Delta t','Callback','plotter_tools(4)','Accelerator','D');
uimenu(f,'Label','Find Delta y','Callback','plotter_tools(10)','Accelerator','F');
uimenu(f,'Label','hh:mm:ss Data Curser','Callback','plotter_tools(5)','Accelerator','K');
uimenu(f,'Label','ss.ssssss Data Curser','Callback','plotter_tools(6)','Accelerator','M');
uimenu(f,'Label','X Zoom!','Callback','plotter_tools(7)','Accelerator','G');
uimenu(f,'Label','Y Zoom!','Callback','plotter_tools(8)','Accelerator','H');
uimenu(f,'Label','XY Zoom!','Callback','plotter_tools(9)','Accelerator','J');
uimenu(f,'Label','Find Position V5','Callback','location_v5');
uimenu(f,'Label','Find Position V7','Callback','location_v7','Accelerator','L');
uimenu(f,'Label','Pulse Modeling','Callback','pulse4','Accelerator','B');
uimenu(f,'Label','Vedio Framing','callback','plotter_tools(11)');
uimenu(f,'Label','Fix for Publishing','callback','fix_fig');
uimenu(f,'Label','Find Rise & Fall Times','callback','rise_fall_times','Accelerator','R');




try
    
    waitbar(0.9,wbh,'Final Preparation','Name','Plotter Busy')
catch
    return
end

if empty_plot==false
    box on
    grid on
    title_str=sprintf('%s-%s-%s    %s   UT: %2.2i:%2.2i:%2.2i  \n ', ...
        g.YYYY{:},g.MM{:},g.DD{:},ld_lin,g.hh,g.mm,g.ss);
    
    if ~isempty(tshift_str)
        title_str = [title_str  tshift_str];
    end
    
    if settings.ldar_tshiftOn==1
        title_str=[title_str '  LDAR:' settings.sen_IDs{sn}];
    end
    
    try
        waitbar(0.95,wbh,'Final Preparation','Name','Plotter Busy')
    catch
        return
    end
    
    
    if isempty(lg)==0
        legend(lg)
    end
    title(title_str)
    
    %     xlab=sprintf('Time (s) \n V shifts : %.2f:%.2f:%.2f -:- %.2f:%.2f:%.2f -:- %.2f:%.2f:%.2f -:- %.2f:%.2f:%.2f -:- %.2f:%.2f:%.2f',...
    %         vshift(1),vshift(2),vshift(3),vshift(4),vshift(5),vshift(6),vshift(7),vshift(8),vshift(9),vshift(10),vshift(11),vshift(12),vshift(13),vshift(14),vshift(15));
    %     xlab=sprintf('%s \n Gains : %.1f:%.1f:%.1f -:- %.1f:%.1f:%.1f -:- %.1f:%.1f:%.1f -:- %.1f:%.1f:%.1f -:- %.1f:%.1f:%.1f',...
    %         xlab,gain(1),gain(2),gain(3),gain(4),gain(5),gain(6),gain(7),gain(8),gain(9),gain(10),gain(11),gain(12),gain(13),gain(14),gain(15));
    
    xlabel('Time (s)')
    
    try
        waitbar(1,wbh,'Final Preparation','Name','Plotter Busy')
    catch
        return
    end
    % set(AX,'PlotBoxAspectRatio',[1 1 1])
    set(AX,'xlim',[g.t1 g.t2])
    
    % Delete wait bar just before plot is visible
    try
        delete(wbh)
    catch
        % Do nothing
    end
    set(fg,'visible','on')
    
else
    try
        delete(wbh)
    catch
        % Do nothing
    end
    errordlg('Plot was empty. May be no data in the given time range!', ...
        'Plotter Error','model')
    
end




function fix_yticks()

ylimits = get(AX1,'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(AX1,'YTick',ylimits(1):yinc:ylimits(2))

ylimits = get(AX2,'YLim');
yinc = (ylimits(2)-ylimits(1))/5;

set(AX2,'YTick',ylimits(1):yinc:ylimits(2))
