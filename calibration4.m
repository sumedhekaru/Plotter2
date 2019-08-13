function varargout = calibration4(varargin)
% CALIBRATION4 M-file for calibration4.fig
%      CALIBRATION4, by itself, creates a new CALIBRATION4 or raises the existing
%      singleton*.
%
%      H = CALIBRATION4 returns the handle to a new CALIBRATION4 or the handle to
%      the existing singleton*.
%
%      CALIBRATION4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALIBRATION4.M with the given input arguments.
%
%      CALIBRATION4('Property','Value',...) creates a new CALIBRATION4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calibration4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calibration4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calibration4

% Last Modified by GUIDE v2.5 05-Apr-2011 22:26:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calibration4_OpeningFcn, ...
                   'gui_OutputFcn',  @calibration4_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before calibration4 is made visible.
function calibration4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calibration4 (see VARARGIN)

g=varargin{1};
handles.g=g;

% Choose default command line output for calibration4
handles.output = hObject;

% Add the standard tool set
set(hObject,'toolbar','figure');

% Disable save figure button in the tool box
set(findall(findall(gcf),'ToolTipString','Save Figure'),'Visible','Off');

% set default values (time shift in micro seconds)
handles.a.sa_gain = 1000;
handles.a.sa_offset = 0;
handles.a.sa_tshift = 0;

handles.a.sa_gain_fine = 0;
handles.a.sa_offset_fine = 0;
handles.a.sa_tshift_fine = -10;

handles.a.sa_gain_fine_max = 100;
handles.a.sa_offset_fine_max = 100;
handles.a.sa_tshift_fine_max = 20;

handles.a.fm_id = '';
handles.a.fm_tshift = -1000;
handles.a.sa_summing = 200;

handles.a.show_points = 0;
handles.a.title ='';
handles.a.lg='';

% data storage
handles.a.sa_t=[];
handles.a.sa_v=[];
handles.a.fm_v=[];
handles.a.fm_t=[];

% set default values
set(handles.sa_gain,'String',handles.a.sa_gain)
set(handles.sa_offset,'String',handles.a.sa_offset)
set(handles.sa_tshift,'String',handles.a.sa_tshift)

set(handles.sa_gain_fine,'String',handles.a.sa_gain_fine)
set(handles.sa_offset_fine,'String',handles.a.sa_offset_fine)
set(handles.sa_tshift_fine,'String',handles.a.sa_tshift_fine)

set(handles.fm_tshift,'String',handles.a.fm_tshift)
set(handles.sum,'String',handles.a.sa_summing)

set(handles.gain_slider,'Value',0.5)
set(handles.offset_slider,'Value',0.5)
set(handles.tshift_slider,'Value',0.25)


% Choose current axis to plot things
axes(handles.axes1);
hold all

[hObject, handles]=load_data(hObject,handles);
set(gca,'xlim',[0 max(handles.a.sa_t)-handles.a.t_min])


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes calibration4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = calibration4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function fm_id_Callback(hObject, eventdata, handles)
% hObject    handle to fm_id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fm_id as text
%        str2double(get(hObject,'String')) returns contents of fm_id as a double
handles.a.fm_id=get(hObject,'String');
guidata(hObject, handles);
[hObject,handles]=load_data(hObject,handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function fm_id_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fm_id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sum_Callback(hObject, eventdata, handles)
% hObject    handle to sum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sum as text
%        str2double(get(hObject,'String')) returns contents of sum as a double
handles.a.sa_summing=str2double(get(hObject,'String'));
guidata(hObject, handles);
[hObject,handles]=load_data(hObject,handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fm_tshift_Callback(hObject, eventdata, handles)
% hObject    handle to fm_tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fm_tshift as text
%        str2double(get(hObject,'String')) returns contents of fm_tshift as a double
handles.a.fm_tshift=str2double(get(hObject,'String'));
guidata(hObject, handles);
[hObject,handles]=load_data(hObject,handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function fm_tshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fm_tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sa_gain_Callback(hObject, eventdata, handles)
% hObject    handle to sa_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sa_gain as text
%        str2double(get(hObject,'String')) returns contents of sa_gain as a double
handles.a.sa_gain=str2double(get(hObject,'String'));
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function sa_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sa_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sa_tshift_Callback(hObject, eventdata, handles)
% hObject    handle to sa_tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sa_tshift as text
%        str2double(get(hObject,'String')) returns contents of sa_tshift as a double
handles.a.sa_tshift=str2double(get(hObject,'String'));
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function sa_tshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sa_tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sa_offset_Callback(hObject, eventdata, handles)
% hObject    handle to sa_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sa_offset as text
%        str2double(get(hObject,'String')) returns contents of sa_offset as a double
handles.a.sa_offset=str2double(get(hObject,'String'));
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function sa_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sa_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sa_gain_fine_Callback(hObject, eventdata, handles)
% hObject    handle to sa_gain_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sa_gain_fine as text
%        str2double(get(hObject,'String')) returns contents of sa_gain_fine as a double
val=str2double(get(hObject,'String'));

if val < 0
    handles.a.sa_gain_fine_max = -val;
    handles.a.sa_gain_fine = val;
    set(handles.gain_slider,'Value',0)
elseif val > 0
    handles.a.sa_gain_fine_max = val;
    handles.a.sa_gain_fine = val;
    set(handles.gain_slider,'Value',1)
else
    handles.a.sa_gain_fine_max = val;
    handles.a.sa_gain_fine = val;
    set(handles.gain_slider,'Value',0.5)   
end

set(handles.auto,'Value',0)

guidata(hObject, handles);
plot_now(hObject,handles)



% --- Executes during object creation, after setting all properties.
function sa_gain_fine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sa_gain_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sa_tshift_fine_Callback(hObject, eventdata, handles)
% hObject    handle to sa_tshift_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sa_tshift_fine as text
%        str2double(get(hObject,'String')) returns contents of sa_tshift_fine as a double
val=str2double(get(hObject,'String'));

if val < 0
    handles.a.sa_tshift_fine_max = -val;
    handles.a.sa_tshift_fine = val;
    set(handles.tshift_slider,'Value',0)
elseif val > 0
    handles.a.sa_tshift_fine_max = val;
    handles.a.sa_tshift_fine = val;
    set(handles.tshift_slider,'Value',1)
else
    handles.a.sa_tshift_fine_max = val;
    handles.a.sa_tshift_fine = val;
    set(handles.tshift_slider,'Value',0.5)   
end

guidata(hObject, handles);
plot_now(hObject,handles)




% --- Executes during object creation, after setting all properties.
function sa_tshift_fine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sa_tshift_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sa_offset_fine_Callback(hObject, eventdata, handles)
% hObject    handle to sa_offset_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sa_offset_fine as text
%        str2double(get(hObject,'String')) returns contents of sa_offset_fine as a double

val=str2double(get(hObject,'String'));

if val < 0
    handles.a.sa_offset_fine_max = -val;
    handles.a.sa_offset_fine = val;
    set(handles.offset_slider,'Value',0)
elseif val > 0
    handles.a.sa_offset_fine_max = val;
    handles.a.sa_offset_fine = val;
    set(handles.offset_slider,'Value',1)
else
    handles.a.sa_offset_fine_max = val;
    handles.a.sa_offset_fine = val;
    set(handles.offset_slider,'Value',0.5)   
end

set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)



% --- Executes during object creation, after setting all properties.
function sa_offset_fine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sa_offset_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_points.
function show_points_Callback(hObject, eventdata, handles)
% hObject    handle to show_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.a.show_points == 0
    handles.a.show_points = 1;
    set(handles.show_points,'String','Hide Points')
else
    handles.a.show_points = 0;
    set(handles.show_points,'String','Show Points')
end
guidata(hObject, handles);
plot_now(hObject,handles)





% --- Executes on slider movement.
function gain_slider_Callback(hObject, eventdata, handles)
% hObject    handle to gain_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=-handles.a.sa_gain_fine_max+2*get(hObject,'Value')*handles.a.sa_gain_fine_max;
handles.a.sa_gain_fine = val;
set(handles.sa_gain_fine,'String',round(val*100)/100)
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function gain_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gain_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function offset_slider_Callback(hObject, eventdata, handles)
% hObject    handle to offset_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=-handles.a.sa_offset_fine_max+2*get(hObject,'Value')*handles.a.sa_offset_fine_max;
handles.a.sa_offset_fine = val;
set(handles.sa_offset_fine,'String',round(val*100)/100)
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function offset_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function tshift_slider_Callback(hObject, eventdata, handles)
% hObject    handle to tshift_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=-handles.a.sa_tshift_fine_max+2*get(hObject,'Value')*handles.a.sa_tshift_fine_max;
handles.a.sa_tshift_fine = val;
set(handles.sa_tshift_fine,'String',round(val*100)/100)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function tshift_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tshift_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function [hObject,handles]=load_data(hObject,handles)

settings=open('sensor_setting.mat');

g=handles.g;
a=g;

bfolder=sprintf('%s/%s/%s/%s/',settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

%% Generating file names for lp graphs and check those files are availble

% checking witch graphs are on
lp_on=g.lpgraphs;

% Save absant file extentions
a.absant_fn={};

% just need to find a single SA file
filename='';

% legend
lg= {};

for i=1:18
    if lp_on(i)==1
        % finding the file extention number
        ext=mod(i,3);
        if ext==0
            ext=3;
        end
        
        sid = settings.sen_IDs{ceil(i/3)};
    
        handles.a.sa_id = sid;
        
        % finding the file name
        filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.lp%1.1i', ...
                bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        % Check whether the file is exists
        if exist(filename,'file')~=0
            % Check whether the file is exists if so load data
            sa_tshift = (handles.a.sa_tshift + handles.a.sa_tshift_fine)/1e3 + settings.t_shift(i);
            
            % Load data
            [sa_t,sa_v]=SA_Extract1(filename,g.t1,g.t2,sa_tshift);
            
            % Number of points that should be summed
            summing = round(handles.a.sa_summing);
            
            ext = sprintf('%1.1i',ext);
            lg =[lg [handles.a.sa_id ':lp' ext]];               
            
            if summing == 1
                handles.a.sa_t = sa_t;
                handles.a.sa_v = sa_v;
                
            elseif summing > 1
                
                % Reduce Frequency of data( (lp3)
                handles.a.sa_t=nan(1,ceil(length(sa_t)/summing));
                handles.a.sa_v=handles.a.sa_t;                               
                
                for k=1:length(handles.a.sa_v)-1
                    handles.a.sa_t(k)=mean(sa_t((k-1)*summing+1:k*summing));
                    handles.a.sa_v(k)= mean(sa_v((k-1)*summing+1:k*summing));
                end
            else
                handles.a.sa_t = sa_t;
                handles.a.sa_v = sa_v;
                fprintf('Data summing failed. Please see the "SA Summing" entry/n')
            end      
            
   
        else
            % If file is not exist don't store the file name
            a.lp_fn{i}='';
            % Absant file
            a.absant_fn=[a.absant_fn filename];
        end
        
        % Exit if you find first SA file name
        
        break
        
    else
        a.lp_fn{i}='';
    end        
end

%% finding the file name for FMs data
if strcmp(handles.a.fm_id,'')
    fm_ID = sprintf('KSC%s',handles.a.sa_id(2:3));
    handles.a.fm_id= fm_ID;
else
    fm_ID = handles.a.fm_id;
end

if a.mm < 30
    ext1=0;
else
    ext1=30;
end

fm_fn=sprintf('%s/FM/%s/%s/%s/%s%s%s%s_%2.2i%2.2i.txt',...
    settings.base_dir,a.YYYY{:},a.MM{:},a.DD{:},fm_ID,...
    a.YYYY{:},a.MM{:},a.DD{:},a.hh,ext1);

fm_tshift = handles.a.fm_tshift/1e3;

if (exist(fm_fn,'file'))
    % if file exist load data
    [handles.a.fm_t,handles.a.fm_v]=fieldMillExtract3(fm_fn,a.t1,a.t2,fm_tshift);
    set(handles.fm_id,'String',fm_ID)
    lg = [lg fm_ID];
    
else
    a.absant_fn=[a.absant_fn fm_fn];    
end
    
%% Letting user know about missing file

a.absant_fn=sort(a.absant_fn);

if isempty (a.absant_fn)==0
    errordlg(a.absant_fn,'Files not found!','modal')
    return
end

handles.a.title=sprintf('SA Calibration  %s-%s-%s    UT: %2.2i:%2.2i:%2.2i',...
    g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss);
handles.a.lg = lg;

plot_now(hObject,handles)
handles.a.t_min = min(handles.a.sa_t);
legend(lg)
xlabel('Time (s)')
ylabel('E (V/m)')
box on
grid on




function plot_now(hObject,handles)
% clear axis
cla
sa_gain = handles.a.sa_gain + handles.a.sa_gain_fine;
sa_offset = handles.a.sa_offset + handles.a.sa_offset_fine;
sa_tshift = (handles.a.sa_tshift + handles.a.sa_tshift_fine)/1e3;

tmin = min(handles.a.sa_t);

if handles.a.show_points == 1
    str1='b-o';
    str2='r-o';
else
    str1='b-';
    str2='r-';
end

plot(handles.a.sa_t-tmin+sa_tshift,handles.a.sa_v*sa_gain+sa_offset,str1)
plot(handles.a.fm_t-tmin,-handles.a.fm_v,str2)

tstr=handles.a.title;
lg = handles.a.lg;

tstr=sprintf('%s\n%s Gain = %.1f m^{-1}    Offset = %.1f Vm^{-1}    T-shift = %.1f ms     %s : T-Shift = %.1f ms',...
    tstr(1,1:42),lg{1},sa_gain,sa_offset,sa_tshift*1e3,lg{2},handles.a.fm_tshift);

title(tstr)


% --- Executes on button press in update_time.
function update_time_Callback(hObject, eventdata, handles)
% hObject    handle to update_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AX=findall(gcf,'Type','axes');
clc
try
    h=guidata(findall(0,'Tag','plotter2'));

    %handle to t1 and t2 in plotter 
    ht1=findall(0,'Tag','t1');
    ht2=findall(0,'Tag','t2');

    xlimits = get(AX(end),'XLim');

    tmin = min(handles.a.sa_t);
    
    h.g.t1=xlimits(1)+tmin;
    h.g.t2=xlimits(2)+tmin;
    
    t1=sprintf('%.6f',h.g.t1);
    t2=sprintf('%.6f',h.g.t2);
    
    set(ht1,'String',t1)
    set(ht2,'String',t2)
    
    %Update plotter 2 handles
    guidata(findall(0,'Tag','plotter2'), h)

catch
    errordlg('Unknown error occured. May be wrong reference for a plot!',...
        'Time Range Error','modal')
end



% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% clear axis
figure
hold on
box on
grid on

sa_gain = handles.a.sa_gain + handles.a.sa_gain_fine;
sa_offset = handles.a.sa_offset + handles.a.sa_offset_fine;
sa_tshift = (handles.a.sa_tshift + handles.a.sa_tshift_fine)/1e3;

tmin = min(handles.a.sa_t);

if handles.a.show_points == 1
    str1='b-o';
    str2='r-o';
else
    str1='b-';
    str2='r-';
end

plot(handles.a.sa_t+sa_tshift,handles.a.sa_v*sa_gain+sa_offset,str1)
plot(handles.a.fm_t,-handles.a.fm_v,str2)

tstr=handles.a.title;
lg = handles.a.lg;

tstr=sprintf('%s\n%s Gain = %.1f m^{-1}    Offset = %.1f Vm^{-1}    T-shift = %.1f ms     %s : T-Shift = %.1f ms',...
    tstr(1,1:42),lg{1},sa_gain,sa_offset,sa_tshift*1e3,lg{2},handles.a.fm_tshift);

set(gca,'xlim',[min(handles.a.sa_t) max(handles.a.sa_t)])

title(tstr)
xlabel('Time (s)')
ylabel('E (v/m)')
legend(lg)


% --- Executes on button press in auto.
function auto_Callback(hObject, eventdata, handles)

g1 = range(handles.a.sa_v);
g2 = range(handles.a.fm_v);
m1 = max(handles.a.sa_v)*g2/g1;
m2 = -min(handles.a.fm_v);

handles.a.sa_gain = round(g2/g1);
handles.a.sa_offset = round(m2-m1);
handles.a.sa_gain_fine = 0;
handles.a.sa_offset_fine = 0;


guidata(hObject, handles);
plot_now(hObject,handles)



set(handles.sa_gain,'String',handles.a.sa_gain)
set(handles.sa_offset,'String',handles.a.sa_offset)
set(handles.sa_gain_fine,'String','0')
set(handles.sa_offset_fine,'String','0')
set(handles.gain_slider,'Value',0.5) 
set(handles.offset_slider,'Value',0.5)

guidata(hObject, handles);
plot_now(hObject,handles)


