function varargout = cal_60Hz_noise(varargin)
% CAL_60HZ_NOISE MATLAB code for cal_60Hz_noise.fig
%      CAL_60HZ_NOISE, by itself, creates a new CAL_60HZ_NOISE or raises the existing
%      singleton*.
%
%      H = CAL_60HZ_NOISE returns the handle to a new CAL_60HZ_NOISE or the handle to
%      the existing singleton*.
%
%      CAL_60HZ_NOISE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAL_60HZ_NOISE.M with the given input arguments.
%
%      CAL_60HZ_NOISE('Property','Value',...) creates a new CAL_60HZ_NOISE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cal_60Hz_noise_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cal_60Hz_noise_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cal_60Hz_noise

% Last Modified by GUIDE v2.5 12-Feb-2015 08:10:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cal_60Hz_noise_OpeningFcn, ...
                   'gui_OutputFcn',  @cal_60Hz_noise_OutputFcn, ...
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


% --- Executes just before cal_60Hz_noise is made visible.
function cal_60Hz_noise_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cal_60Hz_noise (see VARARGIN)

clc
% get Plotter2 data
try 
    h=guidata(findall(0,'Tag','plotter2'));
    handles.g = h.g;
    handles.sen_set = h.sen_set;
catch    
    disp('Run plotter2 first!'); 
end

try
    % Load last data
    handles.a = load('cal_60Hz_noize.mat');
catch
    % Default values;
    handles.a.show_original = 1;
    
    handles.a.amp = 2 ;
    handles.a.fine_amp = 1;
    handles.a.slide_amp = 0;
    
    handles.a.phase = 180;
    handles.a.fine_phase = 180;
    handles.a.slide_phase = 0;
    
    handles.a.freq = 60;
    handles.a.fine_freq = 0.2;
    handles.a.slide_freq = 0;
    
end

set_values(handles);

handles = loadData(handles);

% Choose default command line output for cal_60Hz_noise
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cal_60Hz_noise wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cal_60Hz_noise_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sl_amplitude_Callback(hObject, eventdata, handles)
% hObject    handle to sl_amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(hObject,'Value');
handles.a.slide_amp = val;
handles = plot_corrected_data(handles);

% Update handles structure
guidata(hObject, handles)



% --- Executes during object creation, after setting all properties.
function sl_amplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sl_phase_Callback(hObject, eventdata, handles)

val = get(hObject,'Value');
handles.a.slide_phase = val;
handles = plot_corrected_data(handles);

% Update handles structure
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function sl_phase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function amplitude_Callback(hObject, eventdata, handles)
% hObject    handle to amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of amplitude as text
%        str2double(get(hObject,'String')) returns contents of amplitude as a double


handles.a.amp = str2double(get(hObject,'String'));
handles = plot_corrected_data(handles);

% Update handles structure
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function amplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fn_amplitude_Callback(hObject, eventdata, handles)

val = str2double(get(hObject,'String')); 
handles.a.fine_amp = val;
set(handles.sl_amplitude,'Max',val)
set(handles.sl_amplitude,'Min',-val)
set(handles.sl_amplitude,'Value',0)
handles.a.slide_amp = 0;
handles = plot_corrected_data(handles);

% Update handles structure
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function fn_amplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fn_amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function phase_Callback(hObject, eventdata, handles)

handles.a.phase = str2double(get(hObject,'String'));
handles = plot_corrected_data(handles);

% Update handles structure
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function phase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fn_phase_Callback(hObject, eventdata, handles)

val = str2double(get(hObject,'String')); 
handles.a.fine_phase = val;
set(handles.sl_phase,'Max',val)
set(handles.sl_phase,'Min',-val)
set(handles.sl_phase,'Value',0)
handles.a.slide_phase = 0;
handles = plot_corrected_data(handles);
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function fn_phase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fn_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_original_data.
function show_original_data_Callback(hObject, eventdata, handles)
% hObject    handle to show_original_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_original_data
handles.show_original = get(hObject,'Value');

if handles.show_original  
    d = handles.d;
    axes(handles.axes1)
    hold all
    handles.h_ori_plot = plot(handles.axes1,...
        handles.t,handles.y);
else
    delete(handles.h_ori_plot)
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function handles = loadData(handles)

% If type = 1 then program will force to load NEW data according to
% plotter2 inputs. If not, it will attempt to load data set from hard
% drive.
g = handles.g;

% Turn off unncessary plots
g.linet = 0;
g.ldar = 0;
g.cglss = 0;
g.pbfa = 0;
g.nldn = 0;
g.pbfa_old = 0;
g.chgraphs = zeros(1,60);
lpgraphs = zeros(1,60);

sns = find(g.lpgraphs == 1);

if ~isempty(sns)
    channel = sns(1);
else
    msgbox('No plots were selected, Select a LP plot and re-load data','Noise finder','modal')
    return
end

g.lpgraphs = zeros(1,60);
g.lpgraphs(channel) = 1;
% Load data
d = load_data(g,handles.sen_set);
handles.channel = channel;
handles.t = d.lp_data(channel).t;
handles.y = (d.lp_data(channel).y - d.vshift(channel))/d.gain(channel)/d.factor(channel);

handles.d = d;

hold(handles.axes1,'all')

if handles.a.show_original   
    handles.h_ori_plot = plot(handles.axes1,handles.t,handles.y);
end

xlabel(handles.axes1,'Time (s)')
ylabel(handles.axes1,'Uncalibrated voltage (V)')

handles = plot_corrected_data(handles);



% --- Executes on button press in re_load.
function re_load_Callback(hObject, eventdata, handles)
handles = loadData(handles);

% Update handles structure
guidata(hObject, handles);

function handles = plot_corrected_data(handles)

% sine_wave
A = handles.a.amp + handles.a.slide_amp;
ph = handles.a.phase + handles.a.slide_phase;
f = handles.a.freq + handles.a.slide_freq;

sny = A*sin(2*pi*f*handles.t' + ph/180*pi);

try 
    delete(handles.h_cor_plot)
end

handles.h_cor_plot = plot(handles.axes1,handles.t,handles.y - sny,'r');
title(handles.axes1,sprintf('Frequency = %0.6f Hz     Amplitude = %0.5f V       Phase = %0.2f ',f,A,ph))

function set_values(handles)

set(handles.show_original_data,'Value',handles.a.show_original);
set(handles.amplitude,'String',handles.a.amp)
set(handles.fn_amplitude,'String',handles.a.fine_amp)
set(handles.sl_amplitude,'Value',handles.a.slide_amp)
set(handles.sl_amplitude,'Max',handles.a.fine_amp)
set(handles.sl_amplitude,'Min',-handles.a.fine_amp)

set(handles.phase,'String',handles.a.phase)
set(handles.fn_phase,'String',handles.a.fine_phase)
set(handles.sl_phase,'Value',handles.a.slide_phase)
set(handles.sl_phase,'Max',handles.a.fine_phase)
set(handles.sl_phase,'Min',-handles.a.fine_phase)

set(handles.freq,'String',handles.a.freq)
set(handles.fine_freq,'String',handles.a.fine_freq)
set(handles.sl_freq,'Value',handles.a.slide_freq)
set(handles.sl_freq,'Max',handles.a.fine_freq)
set(handles.sl_freq,'Min',-handles.a.fine_freq)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
a = handles.a;
save('cal_60Hz_noize.mat','-Struct','a')

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in auto.
function auto_Callback(hObject, eventdata, handles)
clc
p = polyfit(handles.t-handles.t(1),handles.y',1);
handles.p = p;
yfit =  p(1) * (handles.t - handles.t(1)) + p(2);
plot(handles.axes1,handles.t,yfit,'k')

figure; hold all
plot(handles.t,handles.y-yfit','k')
[yp, ind] = findpeaks(handles.y-yfit','MINPEAKDISTANCE',20);

range = [];
% get the amplitude
for i = 2:length(ind) - 2;
   range = [range, abs(yp(i+1) - yp(i))/2];
end

A = mean(range);

dxs = nan(1,360);
for i = 1:360
    sny = A*sin(2*pi*60*handles.t + i/180*pi);
    dxs(i) = sum(abs(handles.y-yfit'-sny'));
end

ph = min(dxs(i));

sny = A*sin(2*pi*60*handles.t + ph/180*pi);

plot(handles.t,handles.y-yfit' + sny')


% --- Executes on slider movement.
function sl_freq_Callback(hObject, eventdata, handles)


val = get(hObject,'Value');
handles.a.slide_freq = val;
handles = plot_corrected_data(handles);

% Update handles structure
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function sl_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function freq_Callback(hObject, eventdata, handles)
handles.a.freq = str2double(get(hObject,'String'));
handles = plot_corrected_data(handles);

% Update handles structure
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fine_freq_Callback(hObject, eventdata, handles)

val = str2double(get(hObject,'String')); 
handles.a.fine_freq = val;
set(handles.sl_freq,'Max',val)
set(handles.sl_freq,'Min',-val)
set(handles.sl_freq,'Value',0)
handles.a.slide_freq = 0;
handles = plot_corrected_data(handles);
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function fine_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fine_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
