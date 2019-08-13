function varargout = peak_modifier(varargin)
% PEAK_MODIFIER MATLAB code for peak_modifier.fig
%      PEAK_MODIFIER, by itself, creates a new PEAK_MODIFIER or raises the existing
%      singleton*.
%
%      H = PEAK_MODIFIER returns the handle to a new PEAK_MODIFIER or the handle to
%      the existing singleton*.
%
%      PEAK_MODIFIER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PEAK_MODIFIER.M with the given input arguments.
%
%      PEAK_MODIFIER('Property','Value',...) creates a new PEAK_MODIFIER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before peak_modifier_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to peak_modifier_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help peak_modifier

% Last Modified by GUIDE v2.5 24-Sep-2013 13:57:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @peak_modifier_OpeningFcn, ...
                   'gui_OutputFcn',  @peak_modifier_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before peak_modifier is made visible.
function peak_modifier_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to peak_modifier (see VARARGIN)
%clc

try
    %check whether user asking recalculate peaks
    handles.g = varargin{1};
    handles.pn = handles.g.pn; % peak number
    set(handles.peak_nu,'String',num2str(handles.pn))
    set(handles.n_of_itterations,'String',num2str(handles.g.N))
    set(handles.estimate_error,'Value',handles.g.error_cal)
   
catch
    h=guidata(findall(0,'Tag','plotter2'));
    handles.g = h.g;
    handles.g.N = 100; % Number of iterations for error calculations
    handles.g.t_accu = 0.0000001;
    handles.g.error_cal = 0; % Calculate error is off by default
    handles.pn = 1; % Default peak number
    
end

try
    set(handles.t_inc,'String',handles.g.t_inc)
end

handles.sen_set = open('sensor_setting.mat');

% set user
set(handles.user,'string',handles.sen_set.users,'value',handles.sen_set.userN)

handles.add_btn = 0;
handles.del_btn = 0;

cla(handles.axes1)
%cla(handles.axes2)
%cla(handles.axes3)
%cla(handles.axes4)

handles=load_data(handles);
% Choose default command line output for peak_modifier
handles.output = hObject;

set(handles.info,'String','Important Info Will Display Here!!')

% Open settings file
handles.sen_set = open('sensor_setting.mat');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes peak_modifier wait for user response (see UIRESUME)
% uiwait(handles.peak_modifier);


% --- Outputs from this function are returned to the command line.
function varargout = peak_modifier_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in del_peak.
function del_peak_Callback(hObject, eventdata, handles)

    if handles.add_btn || handles.del_btn
        datacursormode off;
        set(handles.info,'String','Please stop current operation by pressing Reset')
        blink_reset(handles)
        return
    end
    
    handles.del_btn = 1;
    guidata(hObject, handles)
    
    set(handles.info,'String','Select Peak (Red star) to delete')

    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@cur, 'DisplayStyle','window')
    
    x=getCursorInfo(dcm_obj);
    
    while isempty(x) || isempty(find(handles.t == x.Position(1)))
         x=getCursorInfo(dcm_obj);
         pause(0.01)
         if strcmp(get(dcm_obj,'Enable'),'off')
             %disp('Data Curser off - sumedhe')
             return
         end         
    end
    
    index = find(handles.t == x.Position(1));
    
    if ~isempty(index)
        %handles.sen_num(index) = [];
        handles.t(index)       = NaN;
        handles.y(index)       = NaN;
        handles=plot_peaks(handles);
        set(handles.info,'String','Last delete was successful')        
   
    else
        set(handles.info,'String','Couldn''t find the peak!! Please Retry.')        
    end
   
    % Recalculate and show results
    temp_cal_pbfa(handles)
    
    datacursormode off
    handles.del_btn = 0;
    guidata(hObject, handles)
    del_peak_Callback(hObject, eventdata, handles)
    
    

% --- Executes on button press in add_point.
function add_point_Callback(hObject, eventdata, handles)

    if handles.add_btn || handles.del_btn
        datacursormode off;
        set(handles.info,'String','Please stop current operation by pressing Reset')
        blink_reset(handles)
        return
    end
    
    handles.add_btn = 1;
    guidata(hObject, handles)

   
    set(handles.info,'String','Click on plotted data to select a peak')

    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@cur,'DisplayStyle','window')
    
    x=getCursorInfo(dcm_obj);
        
    while isempty(x) ...
            || ( isempty(find(handles.h1 == x.Target,1))...
            && isempty(find(handles.h2 == x.Target,1))...
            && isempty(find(handles.h3 == x.Target,1))...
            && isempty(find(handles.h4 == x.Target,1)))
         x=getCursorInfo(dcm_obj);
         if strcmp(get(dcm_obj,'Enable'),'off')
             %disp('Data Curser off - sumedhe')
             return
         end
         pause(0.01)
    end
    
    %Lets check whether user add the peak for a graph that we don't already
    % have a peak.
    index1 = find(handles.h4 == x.Target);
    
    if isempty(index1);     index1 = find(handles.h1 == x.Target);  end    
    if isempty(index1);     index1 = find(handles.h3 == x.Target);  end
    if isempty(index1);     index1 = find(handles.h2 == x.Target);  end
    
    if ~isnan(handles.t(index1))
        set(handles.info,'String',...
            'A peak already exsist for the selected plot. Please delete that peak and try again.')
        datacursormode off
        return
    end
    
    
    handles.t(index1) = x.Position(1);
    
    y = get(handles.h1(index1),'Ydata');
    handles.y(index1)       = y(x.DataIndex);
    
    guidata(hObject, handles)
    handles=plot_peaks(handles);
    
    datacursormode off;
    set(handles.info,'String','A peak added successfuly')
    
    handles.add_btn = 0;    
    guidata(hObject, handles)
    
    % Calculate and show PBFA for current peaks
    temp_cal_pbfa(handles)
    
    % if we have more points to add, let's continue
    if sum(isnan(handles.t))>=1
        add_point_Callback(hObject, eventdata, handles)
    end
    
    
    
% --- Executes on button press in cal_pbfa.
function cal_pbfa_Callback(hObject, eventdata, handles)

%%%  Input Arguments for pbfa_finder %%%

    arg.inner = handles.sen_num;
    arg.outer = [];
    arg.t_out = [];
    arg.t_in  = handles.t + handles.t_shifts1(handles.sen_num);
    arg.sen_set = handles.sen_set;
    
    % Removing NaNs
    indx = isnan(arg.t_in);
    arg.inner(indx) = [];
    arg.t_in(indx) = [];
    
      
    if length([arg.t_in arg.t_out]) < 5
        errordlg('Atleast 5 sensors need to calculate a PBFA point','PBFA error')
    else
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Method 3 (calculate z using sensors within 6 - 30 km)
       arg.method = 3;

       [x1,y1,z1,t1]=pbfa_finder(arg);

       %N1 = length([arg.t_in arg.t_out])+230;
       % +10 is to understand manual calculation (not automatic)
       
       Ip1 = findPeakCurr(handles,t1,x1,y1,z1,1);
   
       % Ki -sqrd
       ki_sqrd1 = cal_ki_sqrd(t1,x1,y1,z1,arg,handles);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Method 4 (calculate z using all the sensors at > 6 km)
       arg.method = 4;
       
       [x2,y2,z2,t2]=pbfa_finder(arg);
       
       %N2 = length([arg.t_in arg.t_out])+240;
       % +10 is to understand manual calculation (not automatic)
       
       Ip2 = findPeakCurr(handles,t2,x2,y2,z2,0);
   
       % Ki -sqrd
       ki_sqrd2 = cal_ki_sqrd(t2,x2,y2,z2,arg,handles);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Method 5 (Optimize using levenberg-marquardt algorythm)
       arg.method = 5;
       
       [x3,y3,z3,t3]=pbfa_finder(arg);
       
       %N3 = length([arg.t_in arg.t_out])+250;
       % +10 is to understand manual calculation (not automatic)
       
       Ip3 = findPeakCurr(handles,t3,x3,y3,z3,0);
       
       % Ki -sqrd
       ki_sqrd3 = cal_ki_sqrd(t3,x3,y3,z3,arg,handles);
       
       str = sprintf('         Method3  Mothod4  Method5\n=====================================\n');
       str = sprintf('%st (s):   %0.7f    %0.7f    %0.7f\n\n',str,t1,t2,t3);
       
       str = sprintf('%sx (m):   %0.0f    %0.0f    %0.0f\n',str,x1,x2,x3);
       str = sprintf('%sy (m):   %0.0f    %0.0f    %0.0f\n',str,y1,y2,y3);
       str = sprintf('%sz (m):   %0.0f    %0.0f    %0.0f\n\n',str,z1,z2,z3);
       
       str = sprintf('%ski-sqrd:   %0.2f    %0.2f    %0.2f\n',str,ki_sqrd1,ki_sqrd2,ki_sqrd3);
       str = sprintf('%sIp (kA):   %0.2f    %0.2f    %0.2f\n\n',str,Ip1,Ip2,Ip3);
    
       handles.arg = arg;
       guidata(hObject, handles)
       
       if handles.g.error_cal
                     
           errors =pbfa_error_new(handles,handles.g.N);
           
           str = sprintf('%sdt(us):   %0.2f    %0.2f    %0.2f\n',str,errors.dt*1e6);           
           str = sprintf('%sdx(m) :   %0.0f    %0.0f    %0.0f\n',str,errors.dx(1),errors.dx(2),errors.dx(3));
           str = sprintf('%sdy(m) :   %0.0f    %0.0f    %0.0f\n',str,errors.dy(1),errors.dy(2),errors.dy(3));
           str = sprintf('%sdz(m) :   %0.0f    %0.0f    %0.0f\n\n',str,errors.dz(1),errors.dz(2),errors.dz(3));
          
       end
       
       % Check if user want to remove any sensors to improve the results
       removes = remove_check(arg,handles);
       
       str = sprintf('%sImprove ki-sqrd by removing (M5):\n',str);
       for i = 1:length(arg.t_in);
           str = sprintf('%s(%s) %0.2f    ',...
               str,handles.sen_set.sen_IDs{arg.inner(i)},removes.s_ki(i));
       end
       
        %handles = optimize_results(handles);
        %return
       
        button = questdlg(str,'PBFA Finder','Add M3','Add M4','Add M5', 'Add M3');
        
        switch button
            case ''
                % do nothing
            case 'Add M3'
                  add_to_file2(handles,3)
           
            case 'Add M4'
                  add_to_file2(handles,4)
                
            case 'Add M5'
                  add_to_file2(handles,5)
                
        end      
           
    end
    
    
                
            
% --- Executes on button press in add_to_file.
function add_to_file_Callback(hObject, eventdata, handles)

    % calculate pbfa point with 
    handles = cal_to_add2(handles);

    %exit if pbfa has not calculated
    if ~isfield(handles,'pbfa');
       set(handles.info,'String','Have to calculate before add to the file.')
       errordlg('Have to calculate before add')
       return
    end
        
    g=handles.g;
    
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
    %Derectory Name
    dir=sprintf('%s/PBFA/%s/%s/%s/',...
        handles.sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:});
    
    %file name
    fn=sprintf('pbfa_%s%s%s_%2.2i%2.2i.txt',...
        g.YYYY{:},g.MM{:},g.DD{:},g.hh,ext);
    
    %file name to keep peak informations
    fn2=sprintf('pbfa_%s%s%s_%2.2i%2.2i_Peaks.txt',...
        g.YYYY{:},g.MM{:},g.DD{:},g.hh,ext);
    
    % backup file name
    bfn=sprintf('pbfa_%s%s%s_%2.2i%2.2i_backup_%s.txt',...
        g.YYYY{:},g.MM{:},g.DD{:},g.hh,ext,datestr(now,'yyyymmdd_HHMMSS'));
    
 
    % create the folder if it is not excist
    if ~exist(dir,'dir')
        mkdir(dir)
    end
    
    % If not exsist let's write the header
    if ~exist([dir fn],'file')
        fid = fopen([dir fn],'a+');
        fprintf(fid,'ID\tt\t\tx\t\ty\t\tz');
        fprintf(fid,'\n===============================================================');
        fclose(fid);
    else
        % Let's make a backup so that we don't want to rgret
        copyfile([dir fn],[dir bfn])
    end
    
    % Lets read all the data in the file
    fid = fopen([dir fn]);
    data=textscan(fid,'%f %f %f %f %f %f %f %s %f %f %f %f %f','HeaderLines',2);
    fclose(fid);
    
    % make the last column to same as number of sensors
    data{8}=data{6};
    
    data = cell2mat(data);
    
    % find the index
    index = max(data(:,1));
    if isempty(index); index = 0; end
    index = index + 1;
     
    % Check whether user trying to add duplicate point
    index1 = find(data(:,2) == round(handles.pbfa.t*1e8)/1e8);
    
    if length(index1 > 1)
        index1 = index1(end); % Reduce confusions if there is more than 1 duplicated data points 
    end
    
    if ~isempty(index1)
        if (data(index1,3) == round(handles.pbfa.x*10)/10 && ...
                data(index1,4) == round(handles.pbfa.y*10)/10 && ...
                data(index1,5) == round(handles.pbfa.z*10)/10)
            btn = questdlg('Do you really want to add this (DUPLICATE) point?',...
                'Duplicate PBFA Entry','Yes','No','No') ;            
            if ~strcmp(btn,'Yes')
                return
            end
        end
    end
    
    sen_str = sprintf('%i',handles.arg.inner(1));
    for i=2:length(handles.arg.inner)
        sen_str = sprintf('%s,%i',sen_str,handles.arg.inner(i));
    end

    
    
    % Let's add the current data point to the end
    fid = fopen([dir fn],'a+');
    fprintf(fid,'\n%i\t%5.8f\t%.1f\t%.1f\t%.1f\t%i\t%0.1f\t%s\t%0.1f\t%i\t%i\t%i\t%0.1f',...
        index,handles.pbfa.t,handles.pbfa.x,handles.pbfa.y,handles.pbfa.z,...
        handles.pbfa.N,handles.pbfa.Ip,sen_str,handles.pbfa.dt*1e6,round(handles.pbfa.dx), ...
        round(handles.pbfa.dy),round(handles.pbfa.dz),handles.pbfa.ki_sqrd);
    fclose(fid);
    
    % Store peak info
    fid = fopen([dir fn2],'a+');
    
    fprintf(fid,'\n%i',index);    
    t_in = handles.t + handles.t_shifts1(handles.sen_num);
    j = 1;
    
    for i = 1:20
        try
            if i == handles.sen_num(j)
                fprintf(fid,'\t%0.8f',t_in(j));
                j = j + 1;
            else
                fprintf(fid,'\tNaN');
            end
        catch
            fprintf(fid,'\tNaN');
        end
    end
    
    fclose(fid);
    
        
   
    
    
    


% --- Executes on button press in recal_peaks.
function recal_peaks_Callback(hObject, eventdata, handles)
g=handles.g;
% Peak number
g.pn = handles.pn;
delete(gcf)
peak_modifier(g)

function handles = load_data(handles)
g = handles.g;

%clc

% creating struct to store things
a=struct;

a.absant_fn = [];

%% Creating the directory path for all files
% Opening sensor settings
settings=open('sensor_setting.mat');
handles.settings = settings;


% Load ch1 data
handles = load_ch1_data(handles);

% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

%Is manual Time shift for each sensor activated
if settings.man_tshiftOn==1
    tshift=settings.t_shift;
else
    tshift=zeros(1,60);
end


handles=tshifts(handles);
tshift=tshift-handles.t_shifts;


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

%% Generating file names for ch (with header files) graphs and check those files are availble

% checking witch graphs are on
ch_on=g.chgraphs;

if sum(ch_on) < 4
    ch_on = [0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 ...
             0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 ...
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0] ;
end

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


%% Letting user know about missing file

a.absant_fn=sort(a.absant_fn);

if isempty (a.absant_fn)==0
    answer=questdlg(a.absant_fn,'Files not found!','OK','Stop!','OK');
    if strcmp(answer,'Stop!')
        return
    end
end


%% If there is nothing to plot Exit
if sum(strcmp(a.ch_fn,''))== 60 && ...
        strcmp(a.ldar_fn,'')== 1 
    errordlg('No data files were found. All Plot commads failed!', ...
        'Plotter Error','model')
    return
end

wbh= waitbar(0,'Please wait...','Name','Plotter Busy');

% fg=figure;
% set(fg,'visible','off')
hold all


%% Loading and plotting ch data
handles.sen_num = [];
handles.t       = [];
handles.y       = [];
handles.real    = [];
handles.imag    = [];
handles.mag     = [];
handles.maxdV   = [];
handles.h1      = []; %store ch graph handles
handles.h2      = []; %store real part graph handles
handles.h3      = []; %store imaginary part graph handles
handles.h4      = []; %store magnitude graph handles

for i=1:60
    
    try
        waitbar(0.6+i*0.005,wbh,'Loading Fast Antenna data')
    catch
        return
    end
    
    if strcmp(a.ch_fn{i},'')==0

        [t y ch_freq_str] = FA_Extract1(a.ch_fn{i},a.h_fn{i},g.t1,g.t2,tshift(i),settings,i);
        
        if isempty(t)==0
            % plot(t,y*gain(i)*factor(i)+vshift(i))
            ch_legend_str = [ch_legend{i} ch_freq_str];
            lg=[lg ch_legend_str];
            empty_plot=false;
        end
    
       
        if isempty(t)==0
            
            % Lets find the time accuracy (may need for error cal)
            
            t_accu = t(2)-t(1);
            
            if handles.g.t_accu < t_accu
                handles.g.t_accu = t_accu;
            end
            
            
            y_actual = y*gain(i)*factor(i)+vshift(i);

            
            % Replace NaNs with zero so that filtfilt will work fine
            y_actual(isnan(y_actual))=0;
            
            % filter out low frequencies
            Fs = 1/(t(2)-t(1));
            [z,p] = butter(5,10000/(Fs/2),'high'); % Create a low-pass butterworth filter;
            % [z,p,k] = butter(n,Wn) designs an designs an order n lowpass digital Butterworth filter with normalized
            % cutoff frequency Wn. It returns the zeros and poles in length n column
            % vectors z and p, and the gain in the scalar k
            
            smoothy = filtfilt(z,p,y_actual);    % filter the data.
            %smoothy(1,10)     
            % Hilbert transform of voltages
            %y_hil = hilbert(smoothy);
            %y_real = real(y_hil);
            %y_imag = imag(y_hil);
            %y_mag = sqrt(imag(y_hil).^2+real(y_hil).^2);
            
                       
             clear peaks_i
             %[ peaks_i(:,1),  peaks_i(:,2)] = peakfinder(y_mag);

             try
                 % For positive peaks
                 %[peaks_i(:,1),  peaks_i(:,2)] = peakfinder(smoothy,0.01,0.01,1);
                 %peaks_i = sortrows(peaks_i,-2);
%                For negative peaks
                  [peaks_i(:,1),  peaks_i(:,2)] = peakfinder(smoothy,0.01,0.01,-1);
                  peaks_i = sortrows(peaks_i,2);
             catch
                 peaks_i = [];
             end
             
             try
                 ysum = y(peaks_i-100:peaks_i+100);
             catch
                 ysum = y;
             end
             
             
             % Sum data to 1MHz to get peakV values to calculate peak currents
             % find the current frequency
            f=round(1/(t(2)-t(1)));
            % Number of intervals            
            ni = round(f/1e6);
                                   
            if ni > 1
                len = floor(length(ysum)/ni);
                
                temp1 = nan(ni,len);
                temp2 = temp1;
                
                for indn = 1:ni
                    
                    tempt= downsample(t,ni,indn-1);                    
                    temp1(indn,:)=tempt(1:len);             
                    
                    tempy= downsample(ysum,ni,indn-1);
                    temp2(indn,:)=tempy(1:len);
                    
                    clear tempt tempy
                    
                end

                temp2=mean(temp2);
                temp1=mean(temp1);
            
            end
            
           

            % Get the Hilbert transform
            [yF yH] = hill_tra(temp1,temp2,2000);
            
            % Derivative of voltage (or the difference of the voltage)
            temp2 = temp2(2:end)-temp2(1:end-1);
            

            
            maxdV1 = max(abs(yH));
            maxdV2 = maxdV1*gain(i)*factor(i);
            
             if length(peaks_i) >= handles.pn
                peakLoc = peaks_i(handles.pn,1);
               % peakMag = peaks_i(handles.pn,2);
                
                % Store Peak info for future use
                handles.sen_num = [handles.sen_num ceil(i/3)];
                handles.t       = [handles.t t(peakLoc)];
                %handles.y       = [handles.y y_actual(peakLoc)];
                handles.y       = [handles.y smoothy(peakLoc)];
                %handles.real    = [handles.real y_real(peakLoc)];
                %handles.imag    = [handles.imag y_imag(peakLoc)];
                %handles.mag     = [handles.mag y_mag(peakLoc)];
                
                if maxdV2 > 1 
                    handles.maxdV = [handles.maxdV maxdV1];
                else
                    handles.maxdV = [handles.maxdV NaN];
                end
                
                
             else
                % Store Peak info for future use
                handles.sen_num = [handles.sen_num ceil(i/3)];
                handles.t       = [handles.t NaN];
                handles.y       = [handles.y NaN];
                handles.maxdV  = [handles.maxdV NaN];
                %handles.real    = [handles.real NaN];
                %handles.imag    = [handles.imag NaN];
                %handles.mag     = [handles.mag NaN];
                %handles.maxdV   = [handles.maxdV NaN];
             end
             
            % Pliotting data
             %h=plot(handles.axes1,t,y_actual);
             h=plot(handles.axes1,t,smoothy);
             handles.h1 = [handles.h1 h];

            
            % Plotting Real Part         
            %h=plot(handles.axes2,t,y_real);
            %handles.h2 = [handles.h2 h];

            % Plotting Imaginary Part
            %h=plot(handles.axes3,t,y_imag);
            %handles.h3 = [handles.h3 h];            
                        
            % Plotting Magnitude        
            %h=plot(handles.axes4,t,y_mag);
            %handles.h4 = [handles.h4 h];

                      
           empty_plot=false;
        end
    end
end

axes(handles.axes1);
if settings.plot_calibrated == 1
    ylabel('E (V/m)')
else
    ylabel('Voltage (V)')
end


if ~isempty(lg)
  legend(handles.axes1,lg)
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

% Create additional tool menu
tools2fig



try
    
    waitbar(0.9,wbh,'Final Preparation','Name','Plotter Busy')
catch
    return
end

if empty_plot==false
    
    % find all axes handle of type 'axes' and empty tag
    %all_ha = findobj( gcf, 'type', 'axes' );
    all_ha = handles.axes1; %handles.axes2 handles.axes3 handles.axes4];
    linkaxes( all_ha, 'x' );
    
    
    axes(handles.axes1)
    box on
    grid on
    title_str=sprintf('%s-%s-%s    UT: %2.2i:%2.2i:%2.2i  \n ', ...
        g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss);
    
    if ~isempty(tshift_str)
        title_str = [title_str  tshift_str];
    end
    
   
    try
        waitbar(0.95,wbh,'Final Preparation','Name','Plotter Busy')
    catch
        return
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
    set(all_ha,'xlim',[g.t1 g.t2])
    
    % Delete wait bar just before plot is visible
    try
        delete(wbh)
    catch
        % Do nothing
    end

    
    %set(fg,'visible','on')
    
else
    try
        delete(wbh)
    catch
        % Do nothing
    end
    errordlg('Plot was empty. May be no data in the given time range!', ...
        'Plotter Error','model')
    
end

handles=plot_peaks(handles);

% Clculate PBFA and show results for currents peaks
temp_cal_pbfa(handles)





function handles=plot_peaks(handles)

try 
    delete(handles.h5)
end
   
handles.h5 = [];

% Plot Peaks
axes(handles.axes1);
h=plot(handles.t, handles.y,'pr','markerFaceColor','r');
handles.h5 = [handles.h5 h];

% axes(handles.axes2);
% h=plot(handles.t, handles.real,'pr','markerFaceColor','r');
% handles.h5 = [handles.h5 h];

% axes(handles.axes3);
% h=plot(handles.t, handles.imag,'pr','markerFaceColor','r');
% handles.h5 = [handles.h5 h];

% axes(handles.axes4);
% h=plot(handles.t, handles.mag,'pr','markerFaceColor','r');
% handles.h5 = [handles.h5 h];

% Show those peak times in the text box

str='';

for i = 1:length(handles.t)
    if ~isnan(handles.t(i))
        str{i} = sprintf('%s : %5.8fs',...
                handles.sen_set.sen_IDs{handles.sen_num(i)},handles.t(i));
    else
         str{i} = sprintf('%s : Choose a peak',...
                handles.sen_set.sen_IDs{handles.sen_num(i)});
    end
        
end
 
% Display Times
set(handles.time_disp,'String',str')


function output_txt = cur(obj,event_obj,handles)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');

t = pos(1);
%y = pos(2);

%assignin('base', 'peak_modifier_t0', t)
output_txt = sprintf('%5.8fs',t);
% datacursormode off;


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)

pan on
pause(0.05)
set(handles.info,'String','Important info will display here.')
handles.del_btn = 0;
handles.add_btn = 0;
guidata(hObject, handles)

pan off

function blink_reset(handles)

    while handles.del_btn || handles.add_btn        
        set(handles.reset,'BackgroundColor',[1 0 0])
        pause(0.5)
        set(handles.reset,'BackgroundColor',[.941 .941 .941])
        pause(0.5)
        handles=guidata(findall(0,'Tag','peak_modifier'));
    end
    

function t_inc_Callback(hObject, eventdata, handles)
% hObject    handle to t_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_inc as text
%        str2double(get(hObject,'String')) returns contents of t_inc as a double


% --- Executes during object creation, after setting all properties.
function t_inc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)

handles.g.t1 = handles.g.t2;
handles.g.t2 = handles.g.t1+str2double(get(handles.t_inc,'String'))/1e6;
guidata(hObject, handles)
g=handles.g;
g.t_inc = get(handles.t_inc,'String');
g.pn =  handles.pn;
delete(gcf)
peak_modifier2(g)


% --- Executes on button press in prev.
function prev_Callback(hObject, eventdata, handles)

handles.g.t2 = handles.g.t1;
handles.g.t1 = handles.g.t2-str2double(get(handles.t_inc,'String'))/1e6;
guidata(hObject, handles)
g=handles.g;
g.t_inc = get(handles.t_inc,'String');
g.pn = handles.pn;
delete(gcf)
peak_modifier2(g)



function peak_nu_Callback(hObject, eventdata, handles)
handles.pn=str2double(get(hObject,'String')); 
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function peak_nu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peak_nu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = tshifts(handles)
% This function will find closest LDAR point to the begining time of the
% plot. Then using the position of that LDAR point, this will give a time
% shift to each graph. This way peaks might be alligned with each other
% almost.

%% Generating file name for LDAR data
    g = handles.g;
    settings = handles.sen_set;
    
    handles.t_shifts = zeros(1,60);
    handles.t_shifts1 = zeros(1,20);
    
    x0=0;
    y0=0;
    z0=0;
    
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
    
    dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
        -datenum(str2double(g.YYYY),0,0));
    
    ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);
    
    % Let's load extra 5s worth data
    try
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2+5,str2double(settings.ldar_r),...
        x0,y0,z0,0);
    catch
        disp('Warning:: LDAR file not found')
        return
    end
        
    ldar(:,1)=[CG(:,10);DLS(:,10)];         % time
    ldar(:,2)=([CG(:,6);DLS(:,6)]);         % East
    ldar(:,3)=([CG(:,7);DLS(:,7)]);         % North
    ldar(:,4)=([CG(:,8);DLS(:,8)]);         % Altitude
    
    ldar = sortrows(ldar,1);
    
    % Let's return if there is no ldar within 5 seconds selected
    if isnan(ldar(1,1))
        disp('Warning:: There are no LDAR2 data within 5seconds')
        return
    end
    
    fprintf('\nNext LDAR point is at: \n\t\tt=%.7f\t\tx=%0.1fm\t\ty=%0.fm\t\tz=%0.1fm\n',...
                    ldar(1,1),ldar(1,2),ldar(1,3),ldar(1,4))
                
    % Manually set time shift by sumedhe
    %cln = MFileLineNr();
    %fprintf('\n\nSumedhe has sat manual time shift for the flash he is working on.\n')
    %fprintf('To remove it please comment line numbers %i to %i in peak_modifier.m\n\n',cln,cln+6) 
    %ldar(1,2) = -24814;
    %ldar(1,3) = 19899;
    %ldar(1,4) = 4539;

                
    t_shifts=round((sqrt((settings.x-ldar(1,2)).^2 +(settings.y-ldar(1,3)).^2 +...
                    (settings.z-ldar(1,4)).^2)/3e8)*1e7)/1e7;
    
    handles.t_shifts1 = t_shifts;
    
    for i=1:60
        handles.t_shifts(i) = t_shifts(ceil(i/3));
    end
    
    
    
    
    
    
    
    
    


% --- Executes on button press in estimate_error.
function estimate_error_Callback(hObject, eventdata, handles)

    handles.g.error_cal = get(hObject,'Value');
    % Update handles structure
    guidata(hObject, handles);
    




function n_of_itterations_Callback(hObject, eventdata, handles)
    handles.g.N = str2double(get(hObject,'String'));
    % Update handles structure
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function n_of_itterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_of_itterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ip = findPeakCurr(handles,t,x,y,z,show_res)

%sns = handles.sen_num;
%maxdV = handles.maxdV;

% % %slope inrercept for dE peak method
% % m = [ 4.8930 3.9805 4.1938  NaN     2.1098 ...
% %       4.2255 NaN    2.0474  1.6661  2.0610  NaN];
% % c = [0.0382 0.0638  0.2847  NaN     1.2001 ...
% %      0.5778 NaN     1.8081  1.0528  1.1174 NaN];
% 
% % Slope intercept for Hill-Tra peak method
% m =  [2.6509 2.1414 2.3400 NaN NaN 2.2779 3.7111 1.1101 0.9198 1.1545 NaN];
% c = -[0.5544 0.1830 0.3358 NaN NaN 0.2436 0.1679 0.2617 0.2150 0.2912 NaN];
% 
% % Slope intercept for Hill-Tra peak method Using 1.13 power only >30km
% % sensors and y = mx fit
% 
% 
% r = sqrt((handles.sen_set.x - x).^2 + ...
%          (handles.sen_set.y - y).^2 + ...
%          (handles.sen_set.z - z).^2)/1000;
%      
%     sns
% r = r(sns)
% m = m(sns);
% c = c(sns);
% 
% % remove sensors less than 30km
% 
% % indx = find(r > 30);
% % 
% % r = r(indx);
% % c = c(indx);
% % m = m(indx);
% % maxdV = maxdV(indx);

% Ip = ( m.* maxdV .* (r.^1.13) + c)

% if show_res
%     fprintf('\nPeak current calculations:\n')
%     for i=1:length(Ip)
%         fprintf('\t%s\t%0.1f\n',handles.sen_set.sen_IDs{sns(indx(i))},-Ip(i))
%     end
%     
%     fprintf('\tAVG\t%0.1f\n',-nanmean(Ip))
% end

% Ip = -nanmean(Ip);
%Ip = 0;

%% New IP calculations
% Slope intercept for Hill-Tra peak method Using 1.13 power only >30km
% sensors and y = mx fit

% m =  [472.4 387.1 418.7 NaN NaN 411.9 676.4  204.2 165.6 209.0 NaN];
% r = sqrt((handles.sen_set.x - x).^2 + ...
%          (handles.sen_set.y - y).^2 + ...
%          (handles.sen_set.z - z).^2)/1000;
% 
% m = m(sns);
% r = r(sns);
% Ip = m.*maxdV.*((r/100).^1.13);
% indx = find(r > 30);
% Ip = -nanmean(Ip(indx));

%% NEW Ip calculations usin ch1 data
m   = [265.0 219.7 277.0 NaN NaN 247.5 336.5 105.1 84.3  205.1 NaN];
c   = [ -0.5  -0.5  -.7  NaN NaN  -0.5  -0.3  -0.8 -0.5   -0.9 NaN];

r = sqrt((handles.sen_set.x - x).^2 + ...
         (handles.sen_set.y - y).^2 + ...
         (handles.sen_set.z - z).^2);
  
% Let's find the arrival times at different locations
t_ar = t + r/ 2.988957706879363e+08;

% Lets find the HT peak arround this time +/- 2 data points
maxdV = nan(1,11);

% Calibration factors
factor = handles.g.factor(1:3:33);


for i=1:11 % for 11 sensors
    try
        ts = handles.ch1{i}.t;
        
        if ~isempty(ts)
            
            plot(ts,handles.ch1{i}.yH)
            
            lol = sum(ts < t_ar(i)-1e-6);
            ul = sum(ts < t_ar(i)+1e-6);
            
            dVm = max(handles.ch1{i}.yH(lol:ul));
            
            if abs(dVm*factor(i)) > 1
                maxdV(i) = dVm;
                plot(ts(lol),dVm,'o')
            end            
        end
    end
end

Ip = maxdV.*((r(1:11)/100e3).^1.13).*m+c;

indx = r(1:11) > 30;
Ip = -nanmean(Ip(indx));


function handles = cal_to_add(handles)

     cln = MFileLineNr();
    % Who are you? Uncomment the apropriate and comment others
     Iam =  1; % For Sumedhe
%      Iam =  2; % For Nadee
     %Iam =  3;  % For Dr. Marshall
     %Iam =  4; % clay
     
     team = {'Sumedhe' 'Nadee' 'Dr. Marshall'};
     fprintf('\n****************************************************')
     fprintf('\nCurrent user is sat to %s. \nTo change this, change variable ''Iam'' around\n',team{Iam}) 
     fprintf('line number %i in peak_modifier.m.',cln)
     fprintf('\n****************************************************\n\n')
     
     
     %%%  Input Arguments for pbfa_finder %%%

    arg.inner = handles.sen_num;
    arg.outer = [];
    arg.t_out = [];
    arg.t_in  = handles.t + handles.t_shifts1(handles.sen_num);
    arg.method = 3;
    
    % Removing NaNs
    indx = isnan(arg.t_in);
    arg.inner(indx) = [];
    arg.t_in(indx) = [];
    
    
    if length([arg.t_in arg.t_out]) < 5
        errordlg('Atleast 5 sensors need to calculate a PBFA point','PBFA error')
    else
        
       [x,y,z,t]=pbfa_finder(arg);
       
       N = length([arg.t_in arg.t_out])+20+Iam*100;
       % +10 is to understand manual calculation (not automatic)
       
       Ip = findPeakCurr(handles,t,x,y,z,0);
   
       [dx,dy,dz,dt]=pbfa_error(x,y,z,[arg.outer arg.inner],...
           handles.g.t_accu,500,arg.method,2);
       
       str = sprintf('x = (%.1f +/- %0.1f) m\ny = (%.1f  +/- %0.1f) m \nz = (%.1f +/- %0.1f) m \nt = (%6.7f +/- %0.6f) s\nN = %i\nIp = %0.1fkA\nAdd to the file?',...
           x,dx,y,dy,z,dz,t,dt,N,-Ip);
       fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\t%0.0f\t%0.0f\t%0.0f\n',...
           t,x,y,z,N,Ip,dx,dy,dz)
          
        handles.pbfa.x = x;
        handles.pbfa.y = y;
        handles.pbfa.z = real(z);
        handles.pbfa.t = t;
        handles.pbfa.N = N;
        handles.pbfa.Ip = Ip;
        handles.pbfa.dx = dx;
        handles.pbfa.dy = dy;
        handles.pbfa.dz = dz;
        handles.pbfa.dt = dt;
        handles.arg = arg;        
          
    end
  
    
function handles = cal_to_add2(handles)
% This function will calculate pbfa point just before calculate

     cln = MFileLineNr();
    % Who are you? Uncomment the apropriate and comment others
     Iam =  1; % For Sumedhe
    %  Iam =  2; % For Nadee
     %Iam =  3;  % For Dr. Marshall
     %Iam =  4; % clay
     
     team = {'Sumedhe' 'Nadee' 'Dr. Marshall'};
     fprintf('\n****************************************************')
     fprintf('\nCurrent user is sat to %s. \nTo change this, change variable ''Iam'' around\n',team{Iam}) 
     fprintf('line number %i in peak_modifier.m.',cln)
     fprintf('\n****************************************************\n\n')
     
     
     %%%  Input Arguments for pbfa_finder %%%

    arg.inner = handles.sen_num;
    arg.outer = [];
    arg.t_out = [];
    arg.t_in  = handles.t + handles.t_shifts1(handles.sen_num);
        
    % Removing NaNs
    indx = isnan(arg.t_in);
    arg.inner(indx) = [];
    arg.t_in(indx) = [];
    arg.method = 4;
    
   
    if length([arg.t_in arg.t_out]) < 5
        errordlg('Atleast 5 sensors need to calculate a PBFA point','PBFA error')
    else
              
        wbh = waitbar(0,'Ready for calculations','name','PBFA error propagation');
        N0 =  500;
        
        a2.inner = arg.inner;
        a2.outer = [];
        a2.t_in = arg.t_in;
        a2.t_out = [];
        
        a2.method = 3;
        [x03,y03,z03,t03]=pbfa_finder(a2);
        
        a2.method = 4;
        [x04,y04,z04,t04]=pbfa_finder(a2);
        
        dx_temp1 = zeros(1,N0);
        dy_temp1 = dx_temp1;
        dz_temp1 = dx_temp1;
        dt_temp1 = dx_temp1;
        
        dx_temp2 = dx_temp1;
        dy_temp2 = dx_temp1;
        dz_temp2 = dx_temp1;
        dt_temp2 = dx_temp1;
        
        for k = 1:N0

            try
                str = sprintf('Working on %i of %i',k,N0);
                waitbar(k/N0,wbh,str)
            catch
                return
            end
            
            
            % Add randowm numbers between +/- 1/2 time_accuracy to the time
            [raw1 col1] = size(arg.t_in);
            
            %r1 = -t_accu + 2*t_accu.*rand(1,n1);
            r1 =  normrnd(0,handles.g.t_accu,raw1,col1);
            %r2 = -t_accu/2 + t_accu.*rand(1,n2)
            
            
            a2.t_in  = arg.t_in + r1;
            a2.t_out = [];
            
            a2.method = 3;
            [x1,y1,z1,t1]=pbfa_finder(a2);
            
            dx_temp1(k) = (x1-x03);
            dy_temp1(k) = (y1-y03);
            dz_temp1(k) = (z1-z03);
            dt_temp1(k) = (t1);
            
            a2.method = 4;
            [x1,y1,z1,t1]=pbfa_finder(a2);
            
            dx_temp2(k) = (x1-x04);
            dy_temp2(k) = (y1-y04);
            dz_temp2(k) = (z1-z04);
            dt_temp2(k) = (t1);
                       
        end
        
        delete(wbh)
        
        if std(dz_temp1) < std(dz_temp2)
            fprintf('Better when use 6 - 30 km sensors for Z\n')
            x = x03;    y = y03;    z = z03;    t = t03;
            dx = std(dx_temp1);
            dy = std(dy_temp1);
            dz = std(dz_temp1);
            dt = std(dt_temp1);
            N = length([arg.t_in arg.t_out])+30+Iam*100;
            % +10 is to understand manual calculation (not automatic)
        else
            fprintf('Better when use > 6 km sensors for Z\n')
            x = x04;    y = y04;    z = z04;    t = t04;
            dx = std(dx_temp2);
            dy = std(dy_temp2);
            dz = std(dz_temp2);
            dt = std(dt_temp2);
            N = length([arg.t_in arg.t_out])+40+Iam*100;
            % +10 is to understand manual calculation (not automatic)
        end
      
        % Peak Current
        Ip = findPeakCurr(handles,t,x,y,z,0);
        
        
        % Ki -sqrd
        tt1 = arg.t_in';
        
        tt2 = t+sqrt((x - handles.sen_set.x(arg.inner)).^2 + ...
                                         (y - handles.sen_set.y(arg.inner)).^2 + ...
                                         (z - handles.sen_set.z(arg.inner)).^2)'/3e8;
        deg_free = length(arg.t_in)-4;                            
        ki_sqrd = 1/deg_free*sum(((tt1 - tt2)/handles.g.t_accu).^2);
                
      
        
%        str = sprintf('x = (%.1f +/- %0.1f) m\ny = (%.1f  +/- %0.1f) m \nz = (%.1f +/- %0.1f) m \nt = (%6.7f +/- %0.6f) s\nN = %i\nIp = %0.1fkA\nAdd to the file?',...
%            x,dx,y,dy,z,dz,t,dt,N,-Ip);
       fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\t%0.0f\t%0.0f\t%0.0f\t%0.1f\n',...
           t,x,y,z,N,Ip,dx,dy,dz,ki_sqrd)
          
        handles.pbfa.x = x;
        handles.pbfa.y = y;
        handles.pbfa.z = real(z);
        handles.pbfa.t = t;
        handles.pbfa.N = N;
        handles.pbfa.Ip = Ip;
        handles.pbfa.dx = dx;
        handles.pbfa.dy = dy;
        handles.pbfa.dz = dz;
        handles.pbfa.dt = dt;
        handles.pbfa.ki_sqrd = ki_sqrd;
        handles.arg = arg;     
    end
          
function handles = optimize_results(handles)

        x0 = handles.pbfa.x;
        y0 = handles.pbfa.y;
        z0 = handles.pbfa.z;
        t0 = handles.pbfa.t;
        dx = 10;
        dy = 10;
        dz = 10;
        dt = 1e-7;
        
        xs = x0 - 10*dx: dx : x0 + 10*dx;
        ys = y0 - 10*dy: dy : y0 + 10*dy;
        zs = z0 - 20*dz: dz : z0 + 20*dz;
        ts = t0 - 20*dt: dt : t0 + 20*dt;
        
        L1 = length(xs); L2 = length(ys); L3 = length(zs); L4 = length(ts);
        
        kis = NaN(L1,L2,L3,L4);
        
        tt1 = handles.t + handles.t_shifts1(handles.sen_num);
        deg_free = length(tt1)-4;
        
        % disp('working on')
        
        wbh = waitbar(0,'PBFA grid search ....','Name','Optimizing PBFA');
        i = 0;
        tot = L1*L2*L3*L4;
        for m = 1:L1
            x = xs(m);
         
            for n = 1:L2
                y = ys(m);
                for o = 1:L3
                    z = zs(o);
                      waitbar(i/tot,wbh)
                      i = i + L4;
                      
                    for p = 1:L4
                                   
   
            
                        t = ts(p);
                        
                        tt2 = t+sqrt((x - handles.sen_set.x(handles.sen_num)).^2 + ...
                            (y - handles.sen_set.y(handles.sen_num)).^2 + ...
                            (z - handles.sen_set.z(handles.sen_num)).^2)/3e8;                       
                        
                        ki_sqrd = 1/deg_free*sum(((tt1 - tt2)/handles.g.t_accu).^2);
                        
                        kis(m,n,o,p) = ki_sqrd;
                        
                    end
                end
            end
        end
        
        [min_ki, position]= min(kis(:));
        
        [i,j,k,l] = ind2sub(size(kis),position);
        
        fprintf('Optimized: \t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\n',...
                ts(l),xs(i),ys(j),zs(k),min_ki)
            
        delete(wbh)
        
        
function ki_sqrd = cal_ki_sqrd(t,x,y,z,arg,handles)

        % Ki -sqrd
        tt1 = arg.t_in';
        
        deg_free = length(arg.t_in);
        
        tt2 = t+sqrt((x - handles.sen_set.x(arg.inner)).^2 + ...
                                         (y - handles.sen_set.y(arg.inner)).^2 + ...
                                         (z - handles.sen_set.z(arg.inner)).^2)'/3e8;
                          
        ki_sqrd = 1/deg_free*sum(((tt1 - tt2)/handles.g.t_accu).^2);
            

function errors = pbfa_error_new(handles,N)

    a2.inner = handles.arg.inner;
    t_in = handles.arg.t_in;
   
    a2.sen_set = handles.sen_set;
    
    a2.outer = [];
          
    dx_temp1 = NaN(1,N);
    dy_temp1 = dx_temp1;
    dz_temp1 = dx_temp1;
    dt_temp1 = dx_temp1;
    
    dx_temp2 = NaN(1,N);
    dy_temp2 = dx_temp2;
    dz_temp2 = dx_temp2;
    dt_temp2 = dx_temp2;
    
    dx_temp3 = NaN(1,N);
    dy_temp3 = dx_temp3;
    dz_temp3 = dx_temp3;
    dt_temp3 = dx_temp3;
    
    errors = [];
    
    wbh = waitbar(0,'Ready for calculations','name','PBFA error propagation');
    tic
    for k = 1:N
        
        try
            str = sprintf('Working on %i of %i',k,N);
            waitbar(k/N,wbh,str)
        catch
            return
        end
        
        
        % Add randowm numbers between +/- 1/2 time_accuracy to the time
        [raw1 col1] = size(t_in);
        
        r1 =  normrnd(0,handles.g.t_accu,raw1,col1);
        

        
        a2.t_in  = t_in + r1;
        a2.t_out = [];
        
        a2.method = 3;
        
        [x1,y1,z1,t1]=pbfa_finder(a2);
        
        dx_temp1(k) = (x1);
        dy_temp1(k) = (y1);
        dz_temp1(k) = (z1);
        dt_temp1(k) = (t1);
        
        a2.method = 4;
        [x1,y1,z1,t1]=pbfa_finder(a2);
        
        dx_temp2(k) = (x1);
        dy_temp2(k) = (y1);
        dz_temp2(k) = (z1);
        dt_temp2(k) = (t1);
        
        a2.method = 5;
        [x1,y1,z1,t1]=pbfa_finder(a2);
        
        dx_temp3(k) = (x1);
        dy_temp3(k) = (y1);
        dz_temp3(k) = (z1);
        dt_temp3(k) = (t1);
        
        
        a2.method = 6;
        [x1,y1,z1,t1]=pbfa_finder(a2);
        
        dx_temp4(k) = (x1);
        dy_temp4(k) = (y1);
        dz_temp4(k) = (z1);
        dt_temp4(k) = (t1);
        
    end
            
    errors.dx = [std(dx_temp1) std(dx_temp2) std(dx_temp3) std(dx_temp4)];
    errors.dy = [std(dy_temp1) std(dy_temp2) std(dy_temp3) std(dy_temp4)];
    errors.dz = [std(dz_temp1) std(dz_temp2) std(dz_temp3) std(dz_temp4)];
    errors.dt = [std(dt_temp1) std(dt_temp2) std(dt_temp3) std(dt_temp4)]; 
 
    delete(wbh)
           
function removes = remove_check(arg,handles)

    L = length(arg.t_in);
    
    removes.s_ki = NaN(1,L);
    
    if L > 5
        
        arg.method = 5;
        ts = arg.t_in;
        sns = arg.inner;
        
        for i=1:L
            
            arg.t_in = ts;
            arg.inner = sns;
            
            arg.t_in(i) = [];
            arg.inner(i) = [];
            
            [x,y,z,t]=pbfa_finder(arg);
            removes.s_ki(i) = cal_ki_sqrd(t,x,y,z,arg,handles);
        end

    end
    
   

% --- Executes on selection change in user.
function user_Callback(hObject, eventdata, handles)

    sen_set = handles.sen_set;
    sen_set.userN = get(hObject,'Value');
    save('sensor_setting.mat','-Struct','sen_set')
    handles.sen_set = sen_set;
    guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function user_CreateFcn(hObject, eventdata, handles)
% hObject    handle to user (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function add_to_file2(handles,TOAM)
    
    % calculate pbfa point with ths suggested method
    arg = handles.arg;
    arg.method = TOAM;
    [x,y,z,t]=pbfa_finder(arg);
       
    handles.pbfa.N = length([arg.t_in arg.t_out])+10*TOAM+100*handles.sen_set.userN;
       % +10 is to understand manual calculation (not automatic)
       
    handles.pbfa.Ip = findPeakCurr(handles,t,x,y,z,0);
   
    % Ki -sqrd
    handles.pbfa.ki_sqrd = cal_ki_sqrd(t,x,y,z,arg,handles);
    
    handles.arg = arg;
    errors =pbfa_error_new(handles,250);
    
    handles.pbfa.t = t; handles.pbfa.x = x; handles.pbfa.y = y; handles.pbfa.z = z;
    
    handles.pbfa.dt = errors.dt(TOAM - 2);
    handles.pbfa.dx = errors.dx(TOAM - 2); 
    handles.pbfa.dy = errors.dy(TOAM - 2); 
    handles.pbfa.dz = errors.dz(TOAM - 2); 
      
    g=handles.g;
    
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
    %Derectory Name
    dir=sprintf('%s/PBFA/%s/%s/%s/',...
        handles.sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:});
    
    %file name
    fn=sprintf('pbfa_%s%s%s_%2.2i%2.2i.txt',...
        g.YYYY{:},g.MM{:},g.DD{:},g.hh,ext);
    
    %file name to keep peak informations
    fn2=sprintf('pbfa_%s%s%s_%2.2i%2.2i_Peaks.txt',...
        g.YYYY{:},g.MM{:},g.DD{:},g.hh,ext);
    
    % backup file name
    bfn=sprintf('pbfa_%s%s%s_%2.2i%2.2i_backup_%s.txt',...
        g.YYYY{:},g.MM{:},g.DD{:},g.hh,ext,datestr(now,'yyyymmdd_HHMMSS'));
    
 
    % create the folder if it is not excist
    if ~exist(dir,'dir')
        mkdir(dir)
    end
    
    % If not exsist let's write the header
    if ~exist([dir fn],'file')
        fid = fopen([dir fn],'a+');
        fprintf(fid,'ID\tt\t\tx\t\ty\t\tz');
        fprintf(fid,'\n===============================================================');
        fclose(fid);
    else
        % Let's make a backup so that we don't want to rgret
        copyfile([dir fn],[dir bfn])
    end
    
    % Lets read all the data in the file
    fid = fopen([dir fn]);
    data=textscan(fid,'%f %f %f %f %f %f %f %s %f %f %f %f %f','HeaderLines',2);
    fclose(fid);
    
    % make the last column to same as number of sensors
    data{8}=data{6};
    
    data = cell2mat(data);
    
    % find the index
    index = max(data(:,1));
    if isempty(index); index = 0; end
    index = index + 1;
     
    % Check whether user trying to add duplicate point
    index1 = find(data(:,2) == round(handles.pbfa.t*1e8)/1e8);
    
    if length(index1 > 1)
        index1 = index1(end); % Reduce confusions if there is more than 1 duplicated data points 
    end
    
    if ~isempty(index1)
        if (data(index1,3) == round(handles.pbfa.x*10)/10 && ...
                data(index1,4) == round(handles.pbfa.y*10)/10 && ...
                data(index1,5) == round(handles.pbfa.z*10)/10)
            btn = questdlg('Do you really want to add this (DUPLICATE) point?',...
                'Duplicate PBFA Entry','Yes','No','No') ;            
            if ~strcmp(btn,'Yes')
                return
            end
        end
    end
    
    sen_str = sprintf('%i',handles.arg.inner(1));
    for i=2:length(handles.arg.inner)
        sen_str = sprintf('%s,%i',sen_str,handles.arg.inner(i));
    end

    
    
    % Let's add the current data point to the end
    fid = fopen([dir fn],'a+');
    fprintf(fid,'\n%i\t%5.8f\t%.1f\t%.1f\t%.1f\t%i\t%0.1f\t%s\t%0.1f\t%i\t%i\t%i\t%0.1f',...
        index,handles.pbfa.t,handles.pbfa.x,handles.pbfa.y,handles.pbfa.z,...
        handles.pbfa.N,handles.pbfa.Ip,sen_str,handles.pbfa.dt*1e6,round(handles.pbfa.dx), ...
        round(handles.pbfa.dy),round(handles.pbfa.dz),handles.pbfa.ki_sqrd);
    fclose(fid);
    
    % Store peak info
    fid = fopen([dir fn2],'a+');
    
    fprintf(fid,'\n%i',index);    
    t_in = handles.t + handles.t_shifts1(handles.sen_num);
    j = 1;
    
    for i = 1:20
        try
            if i == handles.sen_num(j)
                fprintf(fid,'\t%0.8f',t_in(j));
                j = j + 1;
            else
                fprintf(fid,'\tNaN');
            end
        catch
            fprintf(fid,'\tNaN');
        end
    end
    
    fclose(fid);


% --- Executes on button press in del_all_peaks.
function del_all_peaks_Callback(hObject, eventdata, handles)

try
temp = nan(size(handles.t));

handles.t       = temp;
handles.y       = temp;
handles=plot_peaks(handles);
set(handles.info,'String','All points were delteted')
guidata(hObject, handles);
catch
    set(handles.info,'String','Something went wrong when deleting all points.')
end


% --- Executes on button press in change_peaks.
function change_peaks_Callback(hObject, eventdata, handles)

    if handles.add_btn || handles.del_btn
        datacursormode off;
        set(handles.info,'String','Please stop current operation by pressing Reset')
        blink_reset(handles)
        return
    end
    
    handles.add_btn = 1;
    guidata(hObject, handles)
   
    set(handles.info,'String','Click on plotted data to change a peak')

    datacursormode off;
    datacursormode on;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'UpdateFcn',@cur,'DisplayStyle','window')
    
    x=getCursorInfo(dcm_obj);
        
    while isempty(x) ...
            || ( isempty(find(handles.h1 == x.Target,1))...
            && isempty(find(handles.h2 == x.Target,1))...
            && isempty(find(handles.h3 == x.Target,1))...
            && isempty(find(handles.h4 == x.Target,1)))
         x=getCursorInfo(dcm_obj);
         if strcmp(get(dcm_obj,'Enable'),'off')
             %disp('Data Curser off - sumedhe')
             return
         end
         pause(0.01)
    end
    
    %Lets check whether user add the peak for a graph that we don't already
    % have a peak.
    index1 = find(handles.h4 == x.Target);
    
    if isempty(index1);     index1 = find(handles.h1 == x.Target);  end    
    if isempty(index1);     index1 = find(handles.h3 == x.Target);  end
    if isempty(index1);     index1 = find(handles.h2 == x.Target);  end
    
%     if ~isnan(handles.t(index1))
%         set(handles.info,'String',...
%             'A peak already exsist for the selected plot. Please delete that peak and try again.')
%         datacursormode off
%         return
%     end
    
    
    handles.t(index1) = x.Position(1);
    
    y = get(handles.h1(index1),'Ydata');
    handles.y(index1)       = y(x.DataIndex);
    
    guidata(hObject, handles)
    handles=plot_peaks(handles);
    
    datacursormode off;
    set(handles.info,'String','A peak added successfuly')
    
    handles.add_btn = 0;    
    guidata(hObject, handles)
    
    temp_cal_pbfa(handles)
    % if we have more points to add, let's continue
    change_peaks_Callback(hObject, eventdata, handles)
    
    
%     if sum(isnan(handles.t))>=1
%         add_point_Callback(hObject, eventdata, handles)
%     end

function temp_cal_pbfa(handles)

    arg.inner = handles.sen_num;
    arg.outer = [];
    arg.t_out = [];
    arg.t_in  = handles.t + handles.t_shifts1(handles.sen_num);
    
    arg.sen_set = handles.sen_set;
    
    % Removing NaNs
    indx = isnan(arg.t_in);
    arg.inner(indx) = [];
    arg.t_in(indx) = [];
    
      
    if length([arg.t_in arg.t_out]) < 5
        set(handles.live_cal_text,'String','Atleast 5 peaks needs to estimate position')
    else
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Method 3 (calculate z using sensors within 6 - 30 km)
       arg.method = 3;
       [x1,y1,z1,t1]=pbfa_finder(arg);
       
       %N1 = length([arg.t_in arg.t_out])+230;
       % +10 is to understand manual calculation (not automatic)
       
       Ip1 = findPeakCurr(handles,t1,x1,y1,z1,1);
   
       % Ki -sqrd
       ki_sqrd1 = cal_ki_sqrd(t1,x1,y1,z1,arg,handles);
       
     
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Method 5 (Optimize using levenberg-marquardt algorythm)
       arg.method = 5;
       
       [x3,y3,z3,t3]=pbfa_finder(arg);
       
       %N3 = length([arg.t_in arg.t_out])+250;
       % +10 is to understand manual calculation (not automatic)
       
       Ip3 = findPeakCurr(handles,t3,x3,y3,z3,0);
   
       % Ki -sqrd
       ki_sqrd3 = cal_ki_sqrd(t3,x3,y3,z3,arg,handles);
       
       str = sprintf('         Method3  Method5\n========================\n');
       str = sprintf('%st (s): %0.7f    %0.7f\n',str,t1,t3);       
       str = sprintf('%sx (m): %0.0f    %0.0f\n',str,x1,x3);
       str = sprintf('%sy (m): %0.0f    %0.0f\n',str,y1,y3);
       str = sprintf('%sz (m): %0.0f    %0.0f\n',str,z1,z3);              
       str = sprintf('%ski-sqrd: %0.2f    %0.2f\n',str,ki_sqrd1,ki_sqrd3);
       str = sprintf('%sIp (kA): %0.2f    %0.2f\n\n',str,Ip1,Ip3);
     
       % Check if user want to remove any sensors to improve the results
       removes = remove_check(arg,handles);
       
       str = sprintf('%sImprove ki-sqrd by removing (M5):\n',str);
       for i = 1:length(arg.t_in);
           str = sprintf('%s(%s) %0.2f    ',...
               str,handles.sen_set.sen_IDs{arg.inner(i)},removes.s_ki(i));
       end
       
       set(handles.live_cal_text,'String',str)
           
    end


% --- Executes on button press in cal_PBFA_RS.
function cal_PBFA_RS_Callback(hObject, eventdata, handles)
%%%  Input Arguments for pbfa_finder %%%

    arg.inner = handles.sen_num;
    arg.outer = [];
    arg.t_out = [];
    arg.t_in  = handles.t + handles.t_shifts1(handles.sen_num);
    arg.sen_set = handles.sen_set;
    
    % Removing NaNs
    indx = isnan(arg.t_in);
    arg.inner(indx) = [];
    arg.t_in(indx) = [];
    
      
    if length([arg.t_in arg.t_out]) < 5
        errordlg('Atleast 5 sensors need to calculate a PBFA point','PBFA error')
    else
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Method 3 (calculate z using sensors within 6 - 30 km)
       arg.method = 3;
       [x1,y1,z1,t1]=pbfa_finder(arg);
       
       %N1 = length([arg.t_in arg.t_out])+230;
       % +10 is to understand manual calculation (not automatic)
       
       Ip1 = findPeakCurr(handles,t1,x1,y1,z1,1);
   
       % Ki -sqrd
       ki_sqrd1 = cal_ki_sqrd(t1,x1,y1,z1,arg,handles);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Method 4 (calculate z using all the sensors at > 6 km)
       arg.method = 5;
       
       [x2,y2,z2,t2]=pbfa_finder(arg);
       
       %N2 = length([arg.t_in arg.t_out])+240;
       % +10 is to understand manual calculation (not automatic)
       
       Ip2 = findPeakCurr(handles,t2,x2,y2,z2,0);
   
       % Ki -sqrd
       ki_sqrd2 = cal_ki_sqrd(t2,x2,y2,z2,arg,handles);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Method 5 (Optimize using levenberg-marquardt algorythm)
       arg.method = 6;
       
       [x3,y3,z3,t3]=pbfa_finder(arg);
       
       %N3 = length([arg.t_in arg.t_out])+250;
       % +10 is to understand manual calculation (not automatic)
       
       Ip3 = findPeakCurr(handles,t3,x3,y3,z3,0);
   
       % Ki -sqrd
       ki_sqrd3 = cal_ki_sqrd(t3,x3,y3,z3,arg,handles);
       
       str = sprintf('         Method3  Mothod5  Method6\n=====================================\n');
       str = sprintf('%st (s):   %0.7f    %0.7f    %0.7f\n\n',str,t1,t2,t3);
       
       str = sprintf('%sx (m):   %0.0f    %0.0f    %0.0f\n',str,x1,x2,x3);
       str = sprintf('%sy (m):   %0.0f    %0.0f    %0.0f\n',str,y1,y2,y3);
       str = sprintf('%sz (m):   %0.0f    %0.0f    %0.0f\n\n',str,z1,z2,z3);
       
       str = sprintf('%ski-sqrd:   %0.2f    %0.2f    %0.2f\n',str,ki_sqrd1,ki_sqrd2,ki_sqrd3);
       str = sprintf('%sIp (kA):   %0.2f    %0.2f    %0.2f\n\n',str,Ip1,Ip2,Ip3);
    
       handles.arg = arg;
       guidata(hObject, handles)
       
       if handles.g.error_cal
                     
           errors =pbfa_error_new(handles,handles.g.N);
           
           str = sprintf('%sdt(us):   %0.2f    %0.2f    %0.2f\n',str,errors.dt*1e6);           
           str = sprintf('%sdx(m) :   %0.0f    %0.0f    %0.0f\n',str,errors.dx(1),errors.dx(2),errors.dx(3));
           str = sprintf('%sdy(m) :   %0.0f    %0.0f    %0.0f\n',str,errors.dy(1),errors.dy(2),errors.dy(3));
           str = sprintf('%sdz(m) :   %0.0f    %0.0f    %0.0f\n\n',str,errors.dz(1),errors.dz(2),errors.dz(3));
          
       end
       
       % Check if user want to remove any sensors to improve the results
       removes = remove_check(arg,handles);
       
       str = sprintf('%sImprove ki-sqrd by removing (M5):\n',str);
       for i = 1:length(arg.t_in);
           str = sprintf('%s(%s) %0.2f    ',...
               str,handles.sen_set.sen_IDs{arg.inner(i)},removes.s_ki(i));
       end
       
        %handles = optimize_results(handles);
        %return
       
        button = questdlg(str,'PBFA Finder','Add M3','Add M5','Add M6', 'Add M3');
        
        switch button
            case ''
                % do nothing
            case 'Add M3'
                  add_to_file2(handles,3)
           
            case 'Add M5'
                  add_to_file2(handles,5)
                
            case 'Add M6'
                  add_to_file2(handles,6)
                
        end      
           
    end
    
    
function handles = load_ch1_data(handles)
%% Load ch1 data for peak currrent calculations

settings = handles.settings;
g = handles.g;

% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

ch_on = [1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 ...
             1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 ...
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0] ;
         
         
%Is manual Time shift for each sensor activated
if settings.man_tshiftOn==1
    tshift=settings.t_shift;
else
    tshift=zeros(1,60);
end

ext = 1;
a.ch_fn = [];
a.h_fn = [];

for i=1:60
    if ch_on(i)==1
        
        % Finding the stattion ID
        sid=settings.sen_IDs{ceil(i/3)};
        
        % finding the file name
        filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        
        hfilename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);        
       
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

% Save ch1 data
handles.ch1 = {};

for i=1:3:60        
   
    sn = ceil(i/3);
        
    if strcmp(a.ch_fn{i},'')==0

        [t y] = FA_Extract1(a.ch_fn{i},a.h_fn{i},g.t1,g.t2,tshift(i),settings,i);
        
        if ~isempty(t)          
            
            % Get the Hilbert transform
            [yF yH] = hill_tra(t,y,2000);
            
            handles.ch1{sn}.t = t;
            handles.ch1{sn}.y = y;
            handles.ch1{sn}.yH = abs(yH);
            
        end 
    end
end
    
