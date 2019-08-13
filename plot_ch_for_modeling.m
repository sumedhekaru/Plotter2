function plot_ch_for_modeling

% Let's load the nevest plotter values
try
    h=guidata(findall(0,'Tag','plotter2'));
    % Lets setup calibration values before plot
    h=calibration(h);
catch
    errordlg(sprintf('E-field plotting parameters are \ncoming from Plotter2. \n\nPlease run plotter2 first.\n'),'Plotter2 not found')
    return
end

g = h.g;

a.absant_fn=[];


clc

%LDAR LINET PBFA marker size
mz = 4;

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

wbh= waitbar(0,'Please wait...','Name','Plotter Busy');

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

%% Letting user know about missing file

a.absant_fn=sort(a.absant_fn);

if isempty (a.absant_fn)==0
    answer=questdlg(a.absant_fn,'Files not found!','OK','Stop!','OK');
    if strcmp(answer,'Stop!')
        return
    end
end


%% If there is nothing to plot Exit
if sum(strcmp(a.ch_fn,''))== 60
    errordlg('No data files were found. All Plot commads failed!', ...
        'Plotter Error','model')
    return
end

%% Loading and plotting ch data

global fg
fg=figure;
set(fg,'visible','off')
hold all


% Subplot count 
sbcnt = 0;

% Subplot axises
global ax



ng = 60 - sum(strcmp(a.ch_fn,'')); % Number of graphs (subplots)
ax = nan(1,ng);

% Find number of raws and columns
if      ng == 1;    raw = 1;    col = 1;
elseif  ng == 2;    raw = 1;    col = 2;
elseif  ng <= 4;    raw = 2;    col = 2;
elseif  ng <= 8;    raw = 2;    col = 4;
elseif  ng == 9;    raw = 3;    col = 3;
elseif  ng == 10;   raw = 2;    col =5;
else
    disp('this can not handle more than 10 plots. Exiting...')
    return
end

for i=1:60
    
    try
        waitbar(0.3+i*0.005,wbh,'Loading Fast Antenna data','Name','Plotter Busy')
    catch
        return
    end

    % Plotting up to 60 fast antenna files
    if strcmp(a.ch_fn{i},'')==0
        sbcnt = sbcnt +1;
        subplot(raw,col,sbcnt)
        ax(sbcnt) = gca;

        [t y ch_freq_str] = FA_Extract1(a.ch_fn{i},a.h_fn{i},g.t1,g.t2,tshift(i),settings);
        
        if isempty(t)==0           
            plot(t,y*gain(i)*factor(i)+vshift(i))
            legend(ch_legend{i})
            lg=[lg ch_legend{i}];
            xlim([g.t1 g.t2])
            empty_plot=false;
        end
    end
end

%%
if settings.plot_calibrated == 1
    ylabel('E (V/m)')
else
    ylabel('Voltage (V)')
end

% Create additional tool menu in the plot window

try
    waitbar(0.85,wbh,'Creating Extra Plotting Tools','Name','Plotter Busy')
catch
    return
end

linkaxes(ax,'x')

subplot(raw,col,1)

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
uimenu(f,'Label','Calculate PBFA points','Callback','peak_modifier','Accelerator','L');
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
    %grid on
    title_str=sprintf('%s-%s-%s       UT: %2.2i:%2.2i:%2.2i  \n ', ...
        g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss);
    
    if ~isempty(tshift_str)
        title_str = [title_str  tshift_str];
    end
    
%     if settings.ldar_tshiftOn==1
%         title_str=[title_str '  LDAR:' settings.sen_IDs{sn}];
%     end
    
    try
        waitbar(0.95,wbh,'Final Preparation','Name','Plotter Busy')
    catch
        return
    end
    
    
%     if isempty(lg)==0
%         legend(lg)
%     end
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
    set(gca,'xlim',[g.t1 g.t2])
    
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


