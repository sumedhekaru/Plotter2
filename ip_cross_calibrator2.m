function varargout = ip_cross_calibrator2(varargin)
% IP_CROSS_CALIBRATOR2 MATLAB code for ip_cross_calibrator2.fig
%      IP_CROSS_CALIBRATOR2, by itself, creates a new IP_CROSS_CALIBRATOR2 or raises the existing
%      singleton*.
%
%      H = IP_CROSS_CALIBRATOR2 returns the handle to a new IP_CROSS_CALIBRATOR2 or the handle to
%      the existing singleton*.
%
%      IP_CROSS_CALIBRATOR2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IP_CROSS_CALIBRATOR2.M with the given input arguments.
%
%      IP_CROSS_CALIBRATOR2('Property','Value',...) creates a new IP_CROSS_CALIBRATOR2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ip_cross_calibrator2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ip_cross_calibrator2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ip_cross_calibrator2

% Last Modified by GUIDE v2.5 26-Sep-2012 22:52:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ip_cross_calibrator2_OpeningFcn, ...
                   'gui_OutputFcn',  @ip_cross_calibrator2_OutputFcn, ...
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


% --- Executes just before ip_cross_calibrator2 is made visible.
function ip_cross_calibrator2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ip_cross_calibrator2 (see VARARGIN)
clc
g=varargin{1};
handles.g=g;

% Choose default command line output for ip_cross_calibrator2
handles.output = hObject;
% Add the standard tool set
set(hObject,'toolbar','figure');

% Disable save figure button in the tool box
set(findall(findall(gcf),'ToolTipString','Save Figure'),'Visible','Off');

% set default values (time shift in micro seconds)
handles.a.sa_gain = 1;
handles.a.sa_offset = 0;
handles.a.sa_tshift = 0;

handles.a.sa_gain_fine = 0;
handles.a.sa_offset_fine = 0;
handles.a.sa_tshift_fine = 0;

handles.a.sa_gain_fine_max = 1;
handles.a.sa_offset_fine_max = .1;
handles.a.sa_tshift_fine_max = 20;

handles.a.fm_id = '';
handles.a.fm_tshift = 0;

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

% UIWAIT makes ip_cross_calibrator2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ip_cross_calibrator2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function sa_gain_Callback(hObject, eventdata, handles)
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



function sa_gain_fine_Callback(hObject, eventdata, handles)
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






function sa_offset_Callback(hObject, eventdata, handles)
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



function sa_offset_fine_Callback(hObject, eventdata, handles)
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



function sa_tshift_Callback(hObject, eventdata, handles)
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



function sa_tshift_fine_Callback(hObject, eventdata, handles)
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


% --- Executes on slider movement.
function gain_slider_Callback(hObject, eventdata, handles)
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




% --- Executes on button press in show_points.
function show_points_Callback(hObject, eventdata, handles)
if handles.a.show_points == 0
    handles.a.show_points = 1;
    set(handles.show_points,'String','Hide Points')
else
    handles.a.show_points = 0;
    set(handles.show_points,'String','Show Points')
end
guidata(hObject, handles);
plot_now(hObject,handles)



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

lps = find(lp_on == 1);
sns = ceil(lps/3);
ext = mod(lps,3);
ext(find(ext == 0)) = 3;
sids = settings.sen_IDs(sns);

% legend
lg = sids;

% Load data for first lp file (consider as the reference)
filename = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.lp%1.1i', ...
                bfolder,sids{1},g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext(1));
[handles.a.fm_t,handles.a.fm_v]=SA_Extract1(filename,g.t1,g.t2,settings.t_shift(lps(1)));


% Load data for second lp file (which is going to change)
filename = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.lp%1.1i', ...
                bfolder,sids{2},g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext(2));
sa_tshift = (handles.a.sa_tshift + handles.a.sa_tshift_fine)/1e3 + settings.t_shift(lps(2));
[handles.a.sa_t,handles.a.sa_v]=SA_Extract1(filename,g.t1,g.t2,sa_tshift);


set(handles.uipanel2,'Title',['Controls for ' sids{2}])
    
%% Letting user know about missing file

a.absant_fn=sort(a.absant_fn);

if isempty (a.absant_fn)==0
    errordlg(a.absant_fn,'Files not found!','modal')
    return
end

handles.a.title=sprintf('LP Cross Calibration  %s-%s-%s    UT: %2.2i:%2.2i:%2.2i',...
    g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss);
handles.a.lg = lg;

% do a auto calibration
handles = auto_Callback(hObject, [], handles);

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

plot(handles.a.fm_t-tmin,handles.a.fm_v,str2)
plot(handles.a.sa_t-tmin+sa_tshift,handles.a.sa_v*sa_gain+sa_offset,str1)


tstr=handles.a.title;
lg = handles.a.lg;

 tstr=sprintf('%s\n%s Gain = %.3f m^{-1}    Offset = %.3f Vm^{-1}    T-shift = %.1f ms     %s : T-Shift = %.1f ms',...
     tstr(1,1:42),lg{2},sa_gain,sa_offset,sa_tshift*1e3,lg{2},handles.a.fm_tshift);

title(tstr)

% --- Executes on button press in update_time.
function update_time_Callback(hObject, eventdata, handles)
% hObject    handle to update_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AX=findall(gcf,'Type','axes');

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
plot(handles.a.fm_t,handles.a.fm_v,str2)

tstr=handles.a.title;
lg = handles.a.lg;

tstr=sprintf('%s\n%s Gain = %.3f m^{-1}    Offset = %.3f Vm^{-1}    T-shift = %.1f ms     %s : T-Shift = %.1f ms',...
    tstr(1,1:42),lg{1},sa_gain,sa_offset,sa_tshift*1e3,lg{2},handles.a.fm_tshift);

set(gca,'xlim',[min(handles.a.sa_t) max(handles.a.sa_t)])

title(tstr)
xlabel('Time (s)')
ylabel('E (v/m)')
legend(lg)


% --- Executes on button press in auto.
function handles = auto_Callback(hObject, eventdata, handles)

g1 = range(handles.a.sa_v);
g2 = range(handles.a.fm_v);
m1 = max(handles.a.sa_v)*g2/g1;
m2 = max(handles.a.fm_v);

handles.a.sa_gain = round((g2/g1)*1000)/1000;
handles.a.sa_offset = round((m2-m1)*1000)/1000;
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
