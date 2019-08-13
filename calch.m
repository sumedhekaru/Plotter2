function varargout = calch(varargin)
% CALCH M-file for calch.fig
%      CALCH, by itself, creates a new CALCH or raises the existing
%      singleton*.
%
%      H = CALCH returns the handle to a new CALCH or the handle to
%      the existing singleton*.
%
%      CALCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALCH.M with the given input arguments.
%
%      CALCH('Property','Value',...) creates a new CALCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calch

% Last Modified by GUIDE v2.5 18-Apr-2011 11:34:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calch_OpeningFcn, ...
                   'gui_OutputFcn',  @calch_OutputFcn, ...
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


% --- Executes just before calch is made visible.
function calch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calch (see VARARGIN)
%clc

handles.g=varargin{1};

% Choose default command line output for calch
handles.output = hObject;

% Store the info about lines drew
handles.a.line_h = 0;

% Store data for 3 channels
handles.a.ch1_v=[];
handles.a.ch1_t=[];
handles.a.ch2_v=[];
handles.a.ch2_t=[];
handles.a.ch3_v=[];
handles.a.ch3_t=[];

% Othere important data (sumedhe)
handles.a.ch1_gain = 120;
handles.a.ch2_gain = 1000;
handles.a.ch3_gain = 120;

handles.a.ch1_g_fine = 0;
%handles.a.ch2_g_fine = 10;
handles.a.ch3_g_fine = 0;

handles.a.ch1_g_fine_max = 10;
%handles.a.ch2_g_fine = 10;
handles.a.ch3_g_fine_max = 10;


handles.a.ch1_offset = 0;
handles.a.ch2_offset = 0;
handles.a.ch3_offset = 0;

handles.a.ch1_offset_fine = 0;
%handles.a.ch2_offset_fine = 10;
handles.a.ch3_offset_fine = 0;


handles.a.ch1_offset_fine_max = 10;
%handles.a.ch2_offset_fine = 10;
handles.a.ch3_offset_fine_max = 10;

% Time shifts in micro seconds
handles.a.ch1_tshift = 0;
handles.a.ch1_tshift_fine = 0;
handles.a.ch3_tshift = 0;
handles.a.ch3_tshift_fine = 0;
handles.a.ch1_tshift_fine_max = 10;
handles.a.ch3_tshift_fine_max = 10;

handles.a.show_points = 0;


handles.a.is_norm  = 0;
handles.a.x_norm = NaN;
handles.a.y_norm = NaN;
handles.a.z_norm = NaN;


set(handles.ch1_tshift_slider,'Value',0.5)
set(handles.ch3_tshift_slider,'Value',0.5)



% set Default values
set(handles.ch1_gain,'String',handles.a.ch1_gain)
set(handles.ch2_gain,'String',handles.a.ch2_gain)
set(handles.ch3_gain,'String',handles.a.ch3_gain)

set(handles.ch1_offset,'String',handles.a.ch1_offset)
set(handles.ch2_offset,'String',handles.a.ch2_offset)
set(handles.ch3_offset,'String',handles.a.ch3_offset)

set(handles.ch1_offset_fine,'String',handles.a.ch1_offset_fine)
set(handles.ch3_offset_fine,'String',handles.a.ch3_offset_fine)

set(handles.ch1_g_fine,'String',handles.a.ch1_g_fine)
set(handles.ch3_g_fine,'String',handles.a.ch1_g_fine)

set(handles.ch1_tshift,'String',handles.a.ch1_tshift)
set(handles.ch1_tshift_fine,'String',handles.a.ch1_tshift_fine)

set(handles.ch3_tshift,'String',handles.a.ch3_tshift)
set(handles.ch3_tshift_fine,'String',handles.a.ch3_tshift_fine)

%sliders
set(handles.ch1_offset_slider,'Value',0.5)
set(handles.ch3_offset_slider,'Value',0.5)
set(handles.ch1_gain_slider,'Value',0.5)
set(handles.ch3_gain_slider,'Value',0.5)

% legend
handles.a.legend ='';

% Add the standard tool set
set(hObject,'toolbar','figure');

% Disable save figure button in the tool box
set(findall(findall(gcf),'ToolTipString','Save Figure'),'Visible','Off');


% Load data
[hObject,handles]=load_data(hObject, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes calch wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = calch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function ch1_gain_slider_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_gain_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val=-handles.a.ch1_g_fine_max+2*get(hObject,'Value')*handles.a.ch1_g_fine_max;
handles.a.ch1_g_fine = val;
set(handles.ch1_g_fine,'String',round(val*100)/100)
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ch1_gain_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1_gain_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function ch1_gain_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch1_gain as text
%        str2double(get(hObject,'String')) returns contents of ch1_gain as a double\

handles.a.ch1_gain=str2double(get(hObject,'String'));
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)




% --- Executes during object creation, after setting all properties.
function ch1_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch3_gain_Callback(hObject, eventdata, handles)
% hObject    handle to ch3_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch3_gain as text
%        str2double(get(hObject,'String')) returns contents of ch3_gain as a double
handles.a.ch3_gain=str2double(get(hObject,'String'));
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ch3_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch3_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function ch3_gain_slider_Callback(hObject, eventdata, handles)
% hObject    handle to ch3_gain_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=-handles.a.ch3_g_fine_max+2*get(hObject,'Value')*handles.a.ch3_g_fine_max;
handles.a.ch3_g_fine = val;
set(handles.ch3_g_fine,'String',round(val*100)/100)
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ch3_gain_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch3_gain_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function ch1_g_fine_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_g_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch1_g_fine as text
%        str2double(get(hObject,'String')) returns contents of ch1_g_fine as a double
val=str2double(get(hObject,'String'));

if val < 0
    handles.a.ch1_g_fine_max = -val;
    handles.a.ch1_g_fine = val;
    set(handles.ch1_gain_slider,'Value',0)
elseif val > 0
    handles.a.ch1_g_fine_max = val;
    handles.a.ch1_g_fine = val;
    set(handles.ch1_gain_slider,'Value',1)
else
    handles.a.ch1_g_fine_max = val;
    handles.a.ch1_g_fine = val;
    set(handles.ch1_gain_slider,'Value',0.5)   
end

set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ch1_g_fine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1_g_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch3_g_fine_Callback(hObject, eventdata, handles)
% hObject    handle to ch3_g_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch3_g_fine as text
%        str2double(get(hObject,'String')) returns contents of ch3_g_fine as a double

val=str2double(get(hObject,'String'));

if val < 0
    handles.a.ch3_g_fine_max = -val;
    handles.a.ch3_g_fine = val;
    set(handles.ch3_gain_slider,'Value',0)
elseif val > 0
    handles.a.ch3_g_fine_max = val;
    handles.a.ch3_g_fine = val;
    set(handles.ch3_gain_slider,'Value',1)
else
    handles.a.ch3_g_fine_max = val;
    handles.a.ch3_g_fine = val;
    set(handles.ch3_gain_slider,'Value',0.5)   
end

set(handles.auto,'Value',0)

guidata(hObject, handles);
plot_now(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ch3_g_fine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch3_g_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch1_offset_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch1_offset as text
%        str2double(get(hObject,'String')) returns contents of ch1_offset as a double
handles.a.ch1_offset=str2double(get(hObject,'String'));
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ch1_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch3_offset_Callback(hObject, eventdata, handles)
% hObject    handle to ch3_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch3_offset as text
%        str2double(get(hObject,'String')) returns contents of ch3_offset as a double
handles.a.ch3_offset=str2double(get(hObject,'String'));
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ch3_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch3_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch1_offset_fine_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_offset_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch1_offset_fine as text
%        str2double(get(hObject,'String')) returns contents of ch1_offset_fine as a double

val=str2double(get(hObject,'String'));

if val < 0
    handles.a.ch1_offset_fine_max = -val;
    handles.a.ch1_offset_fine = val;
    set(handles.ch1_offset_slider,'Value',0)
elseif val > 0
    handles.a.ch1_offset_fine_max = val;
    handles.a.ch1_offset_fine = val;
    set(handles.ch1_offset_slider,'Value',1)
else
    handles.a.ch1_offset_fine_max = val;
    handles.a.ch1_offset_fine = val;
    set(handles.ch1_offset_slider,'Value',0.5)   
end

set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ch1_offset_fine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1_offset_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch3_offset_fine_Callback(hObject, eventdata, handles)
% hObject    handle to ch3_offset_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch3_offset_fine as text
%        str2double(get(hObject,'String')) returns contents of ch3_offset_fine as a double

val=str2double(get(hObject,'String'));

if val < 0
    handles.a.ch3_offset_fine_max = -val;
    handles.a.ch3_offset_fine = val;
    set(handles.ch3_offset_slider,'Value',0)
elseif val > 0
    handles.a.ch3_offset_fine_max = val;
    handles.a.ch3_offset_fine = val;
    set(handles.ch3_offset_slider,'Value',1)
else
    handles.a.ch3_offset_fine_max = val;
    handles.a.ch3_offset_fine = val;
    set(handles.ch3_offset_slider,'Value',0.5)   
end
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ch3_offset_fine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch3_offset_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function ch1_offset_slider_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_offset_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=-handles.a.ch1_offset_fine_max+2*get(hObject,'Value')*handles.a.ch1_offset_fine_max;
handles.a.ch1_offset_fine = val;
set(handles.ch1_offset_fine,'String',round(val*100)/100)
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ch1_offset_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1_offset_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function ch3_offset_slider_Callback(hObject, eventdata, handles)
% hObject    handle to ch3_offset_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=-handles.a.ch3_offset_fine_max+2*get(hObject,'Value')*handles.a.ch3_offset_fine_max;
handles.a.ch3_offset_fine = val;
set(handles.ch3_offset_fine,'String',round(val*100)/100)
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ch3_offset_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch3_offset_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function ch2_gain_Callback(hObject, eventdata, handles)
% hObject    handle to ch2_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch2_gain as text
%        str2double(get(hObject,'String')) returns contents of ch2_gain as a double
handles.a.ch2_gain=str2double(get(hObject,'String'));
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ch2_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch2_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch2_offset_Callback(hObject, eventdata, handles)
% hObject    handle to ch2_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch2_offset as text
%        str2double(get(hObject,'String')) returns contents of ch2_offset as a double
handles.a.ch2_offset=str2double(get(hObject,'String'));
set(handles.auto,'Value',0)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ch2_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch2_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


figure

gain1 = handles.a.ch1_gain+handles.a.ch1_g_fine;
gain2 = handles.a.ch2_gain;
gain3 = handles.a.ch3_gain+handles.a.ch3_g_fine;

offset1 = handles.a.ch1_offset+handles.a.ch1_offset_fine;
offset3 = handles.a.ch3_offset+handles.a.ch3_offset_fine;

tshift1 = (handles.a.ch1_tshift + handles.a.ch1_tshift_fine)/1e6;
tshift3 = (handles.a.ch3_tshift + handles.a.ch3_tshift_fine)/1e6;

t_min =min(handles.a.ch1_t);


if handles.a.show_points == 1
    str1='b-o';
    str2='k-o';
    str3='r-o';
else
    str1='b-';
    str2='k-';
    str3='r-';
end


handles.a.line_h(1)=plot(handles.a.ch1_t+tshift1, handles.a.ch1_v*(gain1)+...
                          offset1,str1);

hold all
handles.a.line_h(2)=plot(handles.a.ch2_t,handles.a.ch2_v*(gain2)+...
                          handles.a.ch2_offset,str2);
handles.a.line_h(3)=plot(handles.a.ch3_t+tshift3,handles.a.ch3_v*(gain3)+...
                          offset3,str3);



%2010-08-06 1720 62476.912679 62476.912801

% wbh= waitbar(0,'Please wait...','Name','Plotter Busy');

lg=handles.a.legend;
tstr=handles.a.title;

tstr=sprintf('%s\n%s Gain = %.1f m^{-1}    %s Gain = %.1f m^{-1}    %s Gain = %.1f m^{-1}',...
    tstr(1,1:41),lg{1},gain1,lg{2},gain2,lg{3},gain3);




box on
grid on
xlabel('Time (s)')
xlim([handles.g.t1 handles.g.t2])
ylabel('E (V/m)')
legend(lg)
title(tstr)

% Text box to print data
    h=annotation('textbox',[ 0.75, 0.5, 0.20,0.1]);
    str=sprintf('Ch1 offset = %0.1f v/m\nCh2 Offset = %.1f v/m\nCh3 offset = %.1f v/m\nCh1 T-shift = %.1f us \nCh3 T-shift = %0.1f us',...
        offset1,handles.a.ch2_offset,offset3,tshift1*1e6,tshift3*1e6);
    set(h,'String',str)
    %set(h,'FontSize',a.)
    set(h,'FitBoxToText','on')
    set(h,'VerticalAlignment','middle')
    set(h,'HorizontalAlignment','Left')
   





function [hObject, handles] = load_data(hObject,handles)

settings=open('sensor_setting.mat');


handles.sen_set = settings ;

g=handles.g;

bfolder=sprintf('%s/%s/%s/%s/', ...
                settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

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

%Is manual Time shift for each sensor activated
if settings.man_tshiftOn==1
    tshift=settings.t_shift;
else
    tshift=zeros(1,60);
end

%Is manual gain should be included?
%gain = zeros(1,60) + 1;


% the variable for real legend
lg={};

% Generating Ch Legend names
ch_legend=cell(1,60);

for i=1:20;
%     lp_legend{i*3-2}=[settings.sen_IDs{i} ':lp1'];
%     lp_legend{i*3-1}=[settings.sen_IDs{i} ':lp2'];
%     lp_legend{(i*3)}=[settings.sen_IDs{i} ':lp3'];
    
    ch_legend{i*3-2}=[settings.sen_IDs{i} ':ch1'];
    ch_legend{i*3-1}=[settings.sen_IDs{i} ':ch2'];
    ch_legend{i*3}=[settings.sen_IDs{i} ':ch3'];
end

% Generating file name for LDAR data

if g.mm < 30
    ext=0;
else
    ext=30;
end

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);

% If the file exist let's load data
if exist(ldar_fn,'file')~=0
    % Load ldar data
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        0,0,0,0);
end

% Let's load CG first, if not let's try to load DLS
if ~isnan(CG(1,1))
        
        handles.a.x_norm = CG(1,6);
        handles.a.y_norm = CG(1,7);
        handles.a.z_norm = CG(1,8);
    
        set(handles.x_norm,'String',handles.a.x_norm)
        set(handles.y_norm,'String',handles.a.y_norm)
        set(handles.z_norm,'String',handles.a.z_norm)
        
elseif ~isnan(DLS(1,1))
    
        handles.a.x_norm = DLS(1,6);
        handles.a.y_norm = DLS(1,7);
        handles.a.z_norm = DLS(1,8);
    
        set(handles.x_norm,'String',handles.a.x_norm)
        set(handles.y_norm,'String',handles.a.y_norm)
        set(handles.z_norm,'String',handles.a.z_norm)
       
end


%% Loading and plotting ch data



t=[];   % Variable for time
y=[];   % Variable for voltage
empty_trigs=[]; % Variable for empty triggers

% fg=figure;
% set(fg,'visible','off')

axes(handles.axes1);

hold all

tit= '';

% Counter to save data as GUI data 
k=0;

for i=1:60
    
    %waitbar(0.3+i*0.02,wbh,'Loading Fast Antenna data','Name','Plotter Busy')
    
    t=[];   % Variable for time
    y=[];   % Variable for voltage
    
    if ~strcmp(a.ch_fn{i},'')
        [t,y,ch_freq_str]=FA_Extract1(a.ch_fn{i}, a.h_fn{i} ,g.t1,g.t2,tshift(i),settings,1);
    end
    
    if isempty(t)==0
        
        % Save up to 3 data sets as GUI data
        switch k
            case 0
                handles.a.ch1_v = y;
                handles.a.ch1_t = t;
                k = k+1;
                
                % Sensor number for later use
                handles.a.sn1 = ceil(i/3);
                
                factor = handles.a.ch1_gain + handles.a.ch1_g_fine;
                
            case 1
                handles.a.ch2_v = y;
                handles.a.ch2_t = t;
                k = k+1;
                
                % Sensor number for later use
                handles.a.sn2 = ceil(i/3);
                
                factor = handles.a.ch2_gain;
                
                
            case 2
                handles.a.ch3_v = y;
                handles.a.ch3_t = t;
                k = k+1;
                
%                 disp('test')
%                 size(y)
%                 size(t)
                
                % Sensor number for later use
                handles.a.sn3 = ceil(i/3);
                
                factor = handles.a.ch3_gain + handles.a.ch3_g_fine;
                
                
            otherwise
                % Do nothing
        end
        
        lg=[lg ch_legend{i}];
        
        tit = sprintf('%s %s Gain = %0.1f m^{-1}     ',tit,ch_legend{i},factor);
        
        empty_plot=false;
        
        % No more than 3 graphs can be selected
        if k == 3
            break
            fprintf('Done')
        end
        
    end
end











tit=sprintf('CH Calibration  %s-%s-%s    UT: %2.2i:%2.2i:%2.2i\n%s',...
    g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,tit);

handles.a.title = tit;

ylabel('E (V/m)')
xlabel('Time (s)')
title(tit)

%delete(wbh)
box on
grid on

handles.a.legend = lg;

legend(lg)

plot_now(hObject,handles)
guidata(hObject, handles);


function plot_now(hObject,handles)
%clc
% This script can use to calibrate ch can be used to calibrate ch channels
% to each other. All the inputs will take from the plotter. So run this
% program after the plotter.

% % Delete previosely drawn lines
% try
%     for ind = 1: length(handles.a.line_h)
%         if handles.a.line_h(ind) ~=0
%             delete(handles.a.line_h(ind))
%         end
%     end
% catch
%     fprintf('couldn''t delete last plots!\n')
% end

% Clear the axis
cla

%handles.a

% global g;

% gain1 = handles.a.ch1_gain+handles.a.ch1_g_fine;
% gain2 = handles.a.ch2_gain;
% gain3 = handles.a.ch3_gain+handles.a.ch3_g_fine;

tshift1 = (handles.a.ch1_tshift + handles.a.ch1_tshift_fine)/1e6;
tshift3 = (handles.a.ch3_tshift + handles.a.ch3_tshift_fine)/1e6;

t_min =min(handles.a.ch1_t);


if handles.a.show_points == 1
    str1='b-o';
    str2='k-o';
    str3='r-o';
else
    str1='b-';
    str2='k-';
    str3='r-';
end

          

gain1 = handles.a.ch1_gain+handles.a.ch1_g_fine;
gain2 = handles.a.ch2_gain;
gain3 = handles.a.ch3_gain+handles.a.ch3_g_fine;



handles.a.line_h(1)=plot(handles.a.ch1_t - t_min + tshift1, handles.a.ch1_v*(gain1)+...
                          handles.a.ch1_offset+handles.a.ch1_offset_fine,str1);

hold all
handles.a.line_h(2)=plot(handles.a.ch2_t - t_min,handles.a.ch2_v*(gain2)+...
                          handles.a.ch2_offset,str2);
                     

                      
handles.a.line_h(3)=plot(handles.a.ch3_t - t_min + tshift3,handles.a.ch3_v*(gain3)+...
                          handles.a.ch3_offset+handles.a.ch3_offset_fine,str3);



%2010-08-06 1720 62476.912679 62476.912801

% wbh= waitbar(0,'Please wait...','Name','Plotter Busy');

lg=handles.a.legend;
tstr=handles.a.title;

tstr=sprintf('%s\n%s Gain = %.1f m^{-1}    %s Gain = %.1f m^{-1}    %s Gain = %.1f m^{-1}',...
    tstr(1,1:42),lg{1},gain1,lg{2},gain2,lg{3},gain3);

title(tstr)
legend(lg)

% set(fg,'visible','on')

% Update handles structure
guidata(hObject, handles);


% --- Executes on slider movement.
function ch1_tshift_slider_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_tshift_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=-handles.a.ch1_tshift_fine_max+2*get(hObject,'Value')*handles.a.ch1_tshift_fine_max;
handles.a.ch1_tshift_fine = val;
set(handles.ch1_tshift_fine,'String',round(val*100)/100)
guidata(hObject, handles);
plot_now(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ch1_tshift_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1_tshift_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function ch1_tshift_Callback(hObject, eventdata, handles)
% hObject    handle to textdddd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textdddd as text
%        str2double(get(hObject,'String')) returns contents of textdddd as a double
handles.a.ch1_tshift=str2double(get(hObject,'String'));
guidata(hObject, handles);
plot_now(hObject,handles)


% --- Executes during object creation, after setting all properties.
function textdddd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textdddd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch1_tshift_fine_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_tshift_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch1_tshift_fine as text
%        str2double(get(hObject,'String')) returns contents of ch1_tshift_fine as a double
val=str2double(get(hObject,'String'));

if val < 0
    handles.a.ch1_tshift_fine_max = -val;
    handles.a.ch1_tshift_fine = val;
    set(handles.ch1_tshift_slider,'Value',0)
elseif val > 0
    handles.a.ch1_tshift_fine_max = val;
    handles.a.ch1_tshift_fine = val;
    set(handles.ch1_tshift_slider,'Value',1)
else
    handles.a.ch1_tshift_fine_max = val;
    handles.a.ch1_tshift_fine = val;
    set(handles.ch1_tshift_slider,'Value',0.5) 
end


guidata(hObject, handles);
plot_now(hObject,handles)



% --- Executes during object creation, after setting all properties.
function ch1_tshift_fine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1_tshift_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch3_tshift_Callback(hObject, eventdata, handles)
% hObject    handle to ch3_tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch3_tshift as text
%        str2double(get(hObject,'String')) returns contents of ch3_tshift as a double
handles.a.ch3_tshift=str2double(get(hObject,'String'));
guidata(hObject, handles);
plot_now(hObject,handles)


% --- Executes during object creation, after setting all properties.
function ch3_tshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch3_tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch3_tshift_fine_Callback(hObject, eventdata, handles)
% hObject    handle to ch3_tshift_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch3_tshift_fine as text
%        str2double(get(hObject,'String')) returns contents of ch3_tshift_fine as a double
val=str2double(get(hObject,'String'));

if val < 0
    handles.a.ch1_tshift_fine_max = -val;
    handles.a.ch1_tshift_fine = val;
    set(handles.ch1_tshift_slider,'Value',0)
elseif val > 0
    handles.a.ch1_tshift_fine_max = val;
    handles.a.ch1_tshift_fine = val;
    set(handles.ch1_tshift_slider,'Value',1)
else
    handles.a.ch1_tshift_fine_max = val;
    handles.a.ch1_tshift_fine = val;
    set(handles.ch1_tshift_slider,'Value',0.5) 
end


guidata(hObject, handles);
plot_now(hObject,handles)




% --- Executes during object creation, after setting all properties.
function ch3_tshift_fine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch3_tshift_fine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function ch3_tshift_slider_Callback(hObject, eventdata, handles)
% hObject    handle to ch3_tshift_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=-handles.a.ch3_tshift_fine_max+2*get(hObject,'Value')*handles.a.ch3_tshift_fine_max;
handles.a.ch3_tshift_fine = val;
set(handles.ch3_tshift_fine,'String',round(val*100)/100)
guidata(hObject, handles);
plot_now(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ch3_tshift_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch3_tshift_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
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


% --- Executes during object creation, after setting all properties.
function ch1_tshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1_tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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

    tmin = min(handles.a.ch1_t);
    
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


% --- Executes on button press in auto.
function auto_Callback(hObject, eventdata, handles)
g1 = range(handles.a.ch1_v);
g2 = range(handles.a.ch2_v)*handles.a.ch2_gain;
g3 = range(handles.a.ch3_v);

ch1_g = round(g2/g1);
ch2_g = handles.a.ch2_gain;
ch3_g = round(g2/g3);

m1 = max(handles.a.ch1_v)*ch1_g;
m2 = max(handles.a.ch2_v)*ch2_g;
m3 = max(handles.a.ch3_v)*ch3_g;

handles.a.ch1_gain= ch1_g;
handles.a.ch1_offset = round(m2-m1);
handles.a.ch3_gain = ch3_g;
handles.a.ch3_offset = round(m2-m3);

handles.a.ch1_gain_fine = 0;
handles.a.ch3_gain_fine =0;
handles.a.ch1_offset_fine = 0;
handles.a.ch3_offset_fine = 0;

guidata(hObject, handles);
plot_now(hObject,handles)

set(handles.ch1_gain,'String',handles.a.ch1_gain)
set(handles.ch3_gain,'String',handles.a.ch3_gain)

set(handles.ch1_offset,'String',handles.a.ch1_offset)
set(handles.ch3_offset,'String',handles.a.ch3_offset)

set(handles.ch1_g_fine,'String',handles.a.ch1_gain_fine)
set(handles.ch3_g_fine,'String',handles.a.ch3_gain_fine)

set(handles.ch1_offset_fine,'String',handles.a.ch1_offset_fine)
set(handles.ch3_offset_fine,'String',handles.a.ch3_offset_fine)

set(handles.ch1_gain_slider,'Value',0.5)
set(handles.ch3_gain_slider,'Value',0.5)
set(handles.ch1_offset_slider,'Value',0.5)
set(handles.ch3_offset_slider,'Value',0.5)


% --- Executes on button press in norm.
function norm_Callback(hObject, eventdata, handles)
    
    handles.a.is_norm = get(hObject,'Value');
    
    if handles.a.is_norm == 1
        set(handles.x_norm,'Enable','on')
        set(handles.y_norm,'Enable','on')
        set(handles.z_norm,'Enable','on')
    else
        set(handles.x_norm,'Enable','off')
        set(handles.y_norm,'Enable','off')
        set(handles.z_norm,'Enable','off')
    end

    
    % Plot 100km Normalized data?
    
    nd1 = sqrt((handles.sen_set.x(handles.a.sn1) - handles.a.x_norm)^2+...
        (handles.sen_set.y(handles.a.sn1) - handles.a.y_norm)^2+...
        (handles.sen_set.z(handles.a.sn1) - handles.a.z_norm)^2)/100000;
    nd2 = sqrt((handles.sen_set.x(handles.a.sn2) - handles.a.x_norm)^2+...
        (handles.sen_set.y(handles.a.sn2) - handles.a.y_norm)^2+...
        (handles.sen_set.z(handles.a.sn2) - handles.a.z_norm)^2)/100000;
    nd3 = sqrt((handles.sen_set.x(handles.a.sn3) - handles.a.x_norm)^2+...
        (handles.sen_set.y(handles.a.sn3) - handles.a.y_norm)^2+...
        (handles.sen_set.z(handles.a.sn3) - handles.a.z_norm)^2)/100000;
    

if ~isnan(nd1) && ~isnan(nd2) && ~isnan(nd3)
    
    if handles.a.is_norm
            handles.a.ch1_v = handles.a.ch1_v*nd1;
            handles.a.ch2_v = handles.a.ch2_v*nd2;
            handles.a.ch3_v = handles.a.ch3_v*nd3;
    else
            handles.a.ch1_v = handles.a.ch1_v/nd1;
            handles.a.ch2_v = handles.a.ch2_v/nd2;
            handles.a.ch3_v = handles.a.ch3_v/nd3;
    end
 end

   
    
    plot_now(hObject,handles)
    
    guidata(hObject, handles);


function x_norm_Callback(hObject, eventdata, handles)
handles.a.x_norm = str2double(get(hObject,'String')) ;
str2double(get(hObject,'String')) ;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function x_norm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_norm_Callback(hObject, eventdata, handles)
handles.a.y_norm = str2double(get(hObject,'String')) ;
str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function y_norm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_norm_Callback(hObject, eventdata, handles)
handles.a.z_norm = str2double(get(hObject,'String')) ;
str2double(get(hObject,'String')); 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function z_norm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
