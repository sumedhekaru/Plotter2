%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the updated version of plotter program to handle more sensors
%   Things removed : Notes, ch filters, lp filters
%
%   Things added :
%
%   Author: Sumedhe Karunarathne
%
%   Modification History:
%       06/26/2011  Start building GUI plotter 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function varargout = plotter2(varargin)
% PLOTTER2 MATLAB code for plotter2.fig
%      PLOTTER2, by itself, creates a new PLOTTER2 or raises the existing
%      singleton*.
%
%      H = PLOTTER2 returns the handle to a new PLOTTER2 or the handle to
%      the existing singleton*.
%
%      PLOTTER2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTTER2.M with the given input arguments.
%
%      PLOTTER2('Property','Value',...) creates a new PLOTTER2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plotter2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plotter2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)"
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plotter2

% Last Modified by GUIDE v2.5 10-Aug-2015 12:26:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @plotter2_OpeningFcn, ...
    'gui_OutputFcn',  @plotter2_OutputFcn, ...
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


% --- Executes just before plotter2 is made visible.
function plotter2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plotter2 (see VARARGIN)

g=struct( ...
    'Sid',{'Ks'},'Sidn',{1} ...
    ,'YYYY',{'yyyy'},'YYYYn',{1}...
    ,'MM',{'mm'},'MMn',{1}...
    ,'DD',{'dd'},'DDn',{1}...
    ,'hh',{NaN},'hhn',{1}...
    ,'mm',{NaN},'mmn',{1}...
    ,'ss',{NaN},'ssn',{1}...
    ,'fn',{0},'pn',{0}...
    ,'save_fn',{0},'save_pn',{0}...
    ,'lpgraphs',{zeros(1,60)},'chgraphs',{zeros(1,60)}...
    ,'linet',{0},'ldar',{0},'cglss',{0},'pbfa',{0},'nldn',{0},'pbfa_old',{0}...
    ,'dedt',{0},'dedt_trig',{0},'dedt_pps',{0},'dedt_opt',{0}...
    ,'t1',{0},'t2',{0}...
    ,'t1min',{0},'t2max',{0}...
    ,'isSaved',{1} ...
    ,'nldn2',{0} ...
    ,'KSC_FM',{zeros(1,34)} ...
    ,'version',{'3.0'},'modi_date','2013-July-06' ...
    );

handles.g=g;

handles.sen_set=open('sensor_setting.mat');

try
    % check whether user asking for new plotter
    if varargin{1}==1
        % Do nothing        
    end
catch
    try
        % If not try to load last GUI data
        g=open('plotter_last_gui_data.mat');
        g.isSaved=1;
        g.save_fn = 0;
        handles.g=g;
        set_values(handles);
    catch
        % Do nothing
    end
end

% High pass filter
handles.g.temp.ch_hpf = 0;

handles.sen_set.working_dir = pwd;

set(handles.plot_next,'String',['Next ' num2str(handles.sen_set.next_prev_dt) 's'])
set(handles.plot_prev,'String',['Prev ' num2str(handles.sen_set.next_prev_dt) 's'])


handles=calibration(handles);
set_sensor_ids(handles);
setup_buttons(handles);

% Turn on time shift warning
handles.tshift = csvread('time_shifts_2011.csv',1,0);
handles.tshift_war = 1;


% Undo Redo data
handles.temp.UR_data = [];
handles.temp.UR_index = 0;
handles.temp.UR_cindex =0;

% Computer name
[~, cname] = system('hostname');
handles.cname =  deblank(cname);

% set prefix according to the computer name
handles = set_prefix(handles);

% Show Current radar station on settings
set(handles.radar_station,'Label',['Radara station [' ...
    handles.sen_set.radarStations{handles.sen_set.radarStationID} ']']);

% Choose default command line output for plotter2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plotter2 wait for user response (see UIRESUME)
% uiwait(handles.plotter2);


% --- Outputs from this function are returned to the command line.
function varargout = plotter2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in hour.
function hour_Callback(hObject, eventdata, handles)

g=handles.g;

g.hhn=get(hObject,'Value');
l=get(hObject,'String');
g.hh=str2double(l(g.hhn));

t1minimum=0;
t2maximum=0;

if ~isnan(g.hh)
    t1minimum=3600*g.hh;
    t2maximum=t1minimum+300;
end

if ~isnan(g.mm)
    t1minimum=t1minimum+60*g.mm;
    t2maximum=t1minimum+(5-mod(g.mm,5))*60;
end

if ~isnan(g.ss)
    t1minimum=t1minimum+g.ss;
end

g.t1min=t1minimum;
g.t2max=t2maximum;
g.t1=t1minimum;
g.t2=t2maximum;

set(handles.t1,'String',t1minimum)
set(handles.t2,'String',t2maximum)


t1min_s=sprintf('t Min = %is',t1minimum);
t2max_s=sprintf('t Max = %is',t2maximum);
set(handles.t1min,'String',t1min_s)
set(handles.t2max,'String',t2max_s)

g.isSaved=0;

handles.g=g;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function hour_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in minute.
function minute_Callback(hObject, eventdata, handles)
g=handles.g;

g.mmn=get(hObject,'Value');
l=get(hObject,'String');
g.mm=str2double(l(g.mmn));

t1minimum=0;
t2maximum=0;

if ~isnan(g.hh)
    t1minimum=3600*g.hh;
    t2maximum=t1minimum+300;
end

if ~isnan(g.mm)
    t1minimum=t1minimum+60*g.mm;
    t2maximum=t1minimum+(5-mod(g.mm,5))*60;
end

if ~isnan(g.ss)
    t1minimum=t1minimum+g.ss;
end

g.t1min=t1minimum;
g.t2max=t2maximum;
g.t1=t1minimum;
g.t2=t2maximum;

set(handles.t1,'String',t1minimum)
set(handles.t2,'String',t2maximum)
t1min_s=sprintf('t Min = %is',t1minimum);
t2max_s=sprintf('t Max = %is',t2maximum);
set(handles.t1min,'String',t1min_s)
set(handles.t2max,'String',t2max_s)

g.isSaved=0;

handles.g=g;

guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function minute_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in sec.
function sec_Callback(hObject, eventdata, handles)

g=handles.g;

g.ssn=get(hObject,'Value');
l=get(hObject,'String');
g.ss=str2double(l(g.ssn));

t1minimum=0;
t2maximum=0;

if ~isnan(g.hh)
    t1minimum=3600*g.hh;
    t2maximum=t1minimum+300;
end

if ~isnan(g.mm)
    t1minimum=t1minimum+60*g.mm;
    t2maximum=t1minimum+(5-mod(g.mm,5))*60;
end

if ~isnan(g.ss)
    t1minimum=t1minimum+g.ss;
end

g.t1min=t1minimum;
g.t2max=t2maximum;

set(handles.t1,'String',t1minimum)
set(handles.t2,'String',t2maximum)
t1min_s=sprintf('t Min = %is',t1minimum);
t2max_s=sprintf('t Max = %is',t2maximum);
set(handles.t1min,'String',t1min_s)
set(handles.t2max,'String',t2max_s)

g.isSaved=0;

handles.g=g;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFc

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in year.
function year_Callback(hObject, eventdata, handles)

handles.g.YYYYn=get(hObject,'Value');
l=get(hObject,'String');
handles.g.YYYY=l(handles.g.YYYYn);
handles.g.isSaved=0;
handles.tshift_war = 1; % Turn on time shift warning
% calibration(handles)
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function year_CreateFcn(hObject, eventdata, handles)
% hObject    handle to year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in month.
function month_Callback(hObject, eventdata, handles)

handles.g.MMn=get(hObject,'Value');
l=get(hObject,'String');
handles.g.MM=l(handles.g.MMn);
handles.g.isSaved=0;
handles.tshift_war = 1; % Turn on time shift warning
%calibration(handles)
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function month_CreateFcn(hObject, eventdata, handles)
% hObject    handle to month (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in date.
function date_Callback(hObject, eventdata, handles)
handles.g.DDn=get(hObject,'Value');
l=get(hObject,'String');
handles.g.DD=l(handles.g.DDn);
handles.g.isSaved=0;
handles.tshift_war = 1; % Turn on time shift warning
%calibration(handles)
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function date_CreateFcn(hObject, eventdata, handles)
% hObject    handle to date (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t1_Callback(hObject, eventdata, handles)

value=str2double(get(hObject,'String'));

if value < handles.g.t1min
    value = handles.g.t1min;
    set(handles.t1,'String',value)
end

handles.g.t1=value;

handles.g.isSaved=0;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function t1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t2_Callback(hObject, eventdata, handles)

value=str2double(get(hObject,'String'));

if value > handles.g.t2max
    value = handles.g.t2max;
    set(handles.t2,'String',value)
end

handles.g.t2=value;

handles.g.isSaved=0;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function t2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ch011.
function ch011_Callback(hObject, eventdata, handles)

handles.g.chgraphs(1)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

%disp(handles.g)

% --- Executes on button press in lp011.
function lp011_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(1)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

%disp(handles)


% --- Executes on button press in ch012.
function ch012_Callback(hObject, eventdata, handles)
handles.g.chgraphs(2)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);
%disp(handles.g.chgraphs)

% --- Executes on button press in lp012.
function lp012_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(2)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

%disp(handles.g.lpgraphs)


% --- Executes on button press in ch013.
function ch013_Callback(hObject, eventdata, handles)
handles.g.chgraphs(3)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

% --- Executes on button press in lp013.
function lp013_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(3)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch021.
function ch021_Callback(hObject, eventdata, handles)
handles.g.chgraphs(4)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp021.
function lp021_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(4)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch022.
function ch022_Callback(hObject, eventdata, handles)

handles.g.chgraphs(5)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp022.
function lp022_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(5)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch023.
function ch023_Callback(hObject, eventdata, handles)

handles.g.chgraphs(6)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp023.
function lp023_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(6)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch031.
function ch031_Callback(hObject, eventdata, handles)

handles.g.chgraphs(7)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp031.
function lp031_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(7)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch032.
function ch032_Callback(hObject, eventdata, handles)
handles.g.chgraphs(8)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp032.
function lp032_Callback(hObject, eventdata, handles)
handles.g.lpgraphs(8)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

% --- Executes on button press in ch033.
function ch033_Callback(hObject, eventdata, handles)
handles.g.chgraphs(9)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp033.
function lp033_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(9)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch041.
function ch041_Callback(hObject, eventdata, handles)
handles.g.chgraphs(10)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp041.
function lp041_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(10)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch042.
function ch042_Callback(hObject, eventdata, handles)
handles.g.chgraphs(11)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp042.
function lp042_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(11)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch043.
function ch043_Callback(hObject, eventdata, handles)
handles.g.chgraphs(12)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp043.
function lp043_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(12)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch051.
function ch051_Callback(hObject, eventdata, handles)
handles.g.chgraphs(13)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp051.
function lp051_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(13)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch052.
function ch052_Callback(hObject, eventdata, handles)

handles.g.chgraphs(14)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp052.
function lp052_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(14)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch053.
function ch053_Callback(hObject, eventdata, handles)

handles.g.chgraphs(15)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp053.
function lp053_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(15)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch061.
function ch061_Callback(hObject, eventdata, handles)

handles.g.chgraphs(16)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp061.
function lp061_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(16)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch062.
function ch062_Callback(hObject, eventdata, handles)

handles.g.chgraphs(17)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp062.
function lp062_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(17)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

% --- Executes on button press in ch063.
function ch063_Callback(hObject, eventdata, handles)

handles.g.chgraphs(18)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp063.
function lp063_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(18)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch071.
function ch071_Callback(hObject, eventdata, handles)

handles.g.chgraphs(19)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp071.
function lp071_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(19)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch072.
function ch072_Callback(hObject, eventdata, handles)

handles.g.chgraphs(20)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp072.
function lp072_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(20)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch073.
function ch073_Callback(hObject, eventdata, handles)

handles.g.chgraphs(21)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp073.
function lp073_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(21)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch081.
function ch081_Callback(hObject, eventdata, handles)

handles.g.chgraphs(22)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp081.
function lp081_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(22)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch082.
function ch082_Callback(hObject, eventdata, handles)

handles.g.chgraphs(23)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp082.
function lp082_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(23)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch083.
function ch083_Callback(hObject, eventdata, handles)

handles.g.chgraphs(24)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp083.
function lp083_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(24)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch091.
function ch091_Callback(hObject, eventdata, handles)

handles.g.chgraphs(25)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp091.
function lp091_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(25)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch092.
function ch092_Callback(hObject, eventdata, handles)

handles.g.chgraphs(26)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp092.
function lp092_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(26)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch093.
function ch093_Callback(hObject, eventdata, handles)

handles.g.chgraphs(27)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp093.
function lp093_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(27)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch101.
function ch101_Callback(hObject, eventdata, handles)

handles.g.chgraphs(28)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp101.
function lp101_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(28)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch102.
function ch102_Callback(hObject, eventdata, handles)

handles.g.chgraphs(29)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp102.
function lp102_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(29)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch103.
function ch103_Callback(hObject, eventdata, handles)

handles.g.chgraphs(30)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp103.
function lp103_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(30)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

% --- Executes on button press in ch111.
function ch111_Callback(hObject, eventdata, handles)
handles.g.chgraphs(31)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp111.
function lp111_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(31)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch112.
function ch112_Callback(hObject, eventdata, handles)
handles.g.chgraphs(32)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp112.
function lp112_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(32)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch113.
function ch113_Callback(hObject, eventdata, handles)
handles.g.chgraphs(33)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);



% --- Executes on button press in lp113.
function lp113_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(33)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch121.
function ch121_Callback(hObject, eventdata, handles)
handles.g.chgraphs(34)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp121.
function lp121_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(34)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch122.
function ch122_Callback(hObject, eventdata, handles)

handles.g.chgraphs(35)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp122.
function lp122_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(35)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch123.
function ch123_Callback(hObject, eventdata, handles)
handles.g.chgraphs(36)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp123.
function lp123_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(36)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch131.
function ch131_Callback(hObject, eventdata, handles)
handles.g.chgraphs(37)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp131.
function lp131_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(37)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch132.
function ch132_Callback(hObject, eventdata, handles)
handles.g.chgraphs(38)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp132.
function lp132_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(38)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch133.
function ch133_Callback(hObject, eventdata, handles)
handles.g.chgraphs(39)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

% --- Executes on button press in lp133.
function lp133_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(39)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch141.
function ch141_Callback(hObject, eventdata, handles)
handles.g.chgraphs(40)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp141.
function lp141_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(40)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch142.
function ch142_Callback(hObject, eventdata, handles)
handles.g.chgraphs(41)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

% --- Executes on button press in lp142.
function lp142_Callback(hObject, eventdata, handles)
handles.g.lpgraphs(41)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch143.
function ch143_Callback(hObject, eventdata, handles)
handles.g.chgraphs(42)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp143.
function lp143_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(42)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch151.
function ch151_Callback(hObject, eventdata, handles)
handles.g.chgraphs(43)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

% --- Executes on button press in lp151.
function lp151_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(43)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch152.
function ch152_Callback(hObject, eventdata, handles)
handles.g.chgraphs(44)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp152.
function lp152_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(44)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch153.
function ch153_Callback(hObject, eventdata, handles)
handles.g.chgraphs(45)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp153.
function lp153_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(45)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch161.
function ch161_Callback(hObject, eventdata, handles)
handles.g.chgraphs(46)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp161.
function lp161_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(46)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch162.
function ch162_Callback(hObject, eventdata, handles)
handles.g.chgraphs(47)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp162.
function lp162_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(47)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch163.
function ch163_Callback(hObject, eventdata, handles)
handles.g.chgraphs(48)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp163.
function lp163_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(48)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch171.
function ch171_Callback(hObject, eventdata, handles)
handles.g.chgraphs(49)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp171.
function lp171_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(49)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch172.
function ch172_Callback(hObject, eventdata, handles)
handles.g.chgraphs(50)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp172.
function lp172_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(50)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch173.
function ch173_Callback(hObject, eventdata, handles)
handles.g.chgraphs(51)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp173.
function lp173_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(51)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch181.
function ch181_Callback(hObject, eventdata, handles)
handles.g.chgraphs(52)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp181.
function lp181_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(52)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch182.
function ch182_Callback(hObject, eventdata, handles)
handles.g.chgraphs(53)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp182.
function lp182_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(53)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch183.
function ch183_Callback(hObject, eventdata, handles)
handles.g.chgraphs(54)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp183.
function lp183_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(54)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch191.
function ch191_Callback(hObject, eventdata, handles)

handles.g.chgraphs(55)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp191.
function lp191_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(55)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch192.
function ch192_Callback(hObject, eventdata, handles)

handles.g.chgraphs(56)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp192.
function lp192_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(56)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch193.
function ch193_Callback(hObject, eventdata, handles)
handles.g.chgraphs(57)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp193.
function lp193_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(57)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch201.
function ch201_Callback(hObject, eventdata, handles)
handles.g.chgraphs(58)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp201.
function lp201_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(58)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch202.
function ch202_Callback(hObject, eventdata, handles)
handles.g.chgraphs(59)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp202.
function lp202_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(59)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in ch203.
function ch203_Callback(hObject, eventdata, handles)
handles.g.chgraphs(60)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in lp203.
function lp203_Callback(hObject, eventdata, handles)

handles.g.lpgraphs(60)=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

sen_set=open('sensor_setting.mat');
%disp(sen_set)



% --- Executes on button press in ldar.
function ldar_Callback(hObject, eventdata, handles)

handles.g.ldar=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in linet.
function linet_Callback(hObject, eventdata, handles)
handles.g.linet=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in pbfa.
function pbfa_Callback(hObject, eventdata, handles)
handles.g.pbfa=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

% --- Executes on button press in nldn.
function nldn_Callback(hObject, eventdata, handles)
handles.g.nldn=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in plot_all.
function plot_all_Callback(hObject, eventdata, handles)

g=handles.g;

%Is user forget to enter correct date?
if strcmp(g.YYYY,'yyyy') || strcmp(g.DD,'dd') || strcmp(g.MM,'mm')
    errordlg('Check Date input!','Date Error','modal')
    uiwait
    return
end


% Is user forget to enter correct time?
if isnan(g.hh) || isnan(g.mm) || isnan(g.ss)
    errordlg('Check Time input!','Time Error','modal')
    uiwait
    return
end

if g.t1 > g.t2 || isnan(g.t1) || isnan(g.t2)
    t1=g.t1;
    t2=g.t2;
    errordlg('Check Plotting Time Range!','Time range Error','modal')
    uiwait
    return
end



% Is user forget to click on any graphs?
if (sum(g.chgraphs)+sum(g.lpgraphs)+ g.KSC_FM+g.ldar + g.cglss + g.linet + g.pbfa ...
        + g.nldn + g.dedt + g.dedt_trig + g.nldn2) == 0
    errordlg('No graphs were selected to plot!','Graph selection Error','modal')
    uiwait
    return
end


% Has user accidantly trying to plot large range of ch data?
if nnz(g.chgraphs)~= 0 && (g.t2 - g.t1) > 10
    large_data_msg = 'Trying to plot large range of ch data will freez your computer. Do you really want to continue?';
    btn1 = questdlg(large_data_msg,'Warning!','Yes','No','No');
    if ~strcmp(btn1,'Yes')
        return
    end
end

% Lets setup calibration values before plot
handles=calibration(handles);

% Let's remember last plot graphs
try
   if ~isequal(g , handles.temp.g2)
       handles.temp.g1 = handles.temp.g2;
       handles.temp.g2 = g;
   end
    
catch
    
   handles.temp.g1 = g;
   handles.temp.g2 = g;
   
end

handles = tshift_warning(handles);

% Let's save g data for Undo-Redo opration
% if handles.temp.UR_index == 0
%     handles.temp.UR_data = [handles.temp.UR_data g];
%     handles.temp.UR_index = handles.temp.UR_index + 1;
%     handles.temp.UR_cindex = handles.temp.UR_index;     
% elseif ~isequal(handles.temp.UR_data(handles.temp.UR_index),g)
%     handles.temp.UR_data = [handles.temp.UR_data g];
%     handles.temp.UR_index = handles.temp.UR_index + 1;
%     handles.temp.UR_cindex = handles.temp.UR_index;    
% end  
% 
% if handles.temp.UR_cindex > 1
%     set(handles.undo,'Enable','On')
% end

% Setup auto base folder
if strcmp(handles.sen_set.autoBaseFolder,'on')  
    handles = set_auto_b_folder(handles);
end

plot_all6(handles.g);

guidata(hObject, handles);



% --- Executes on button press in plot_next.
function plot_next_Callback(hObject, eventdata, handles)
g=handles.g;

dts = handles.sen_set.next_prev_dt;

if g.t2 == g.t2max
    t = g.t2;
else
    t = g.t1;
end

% there was extra 300s shift when dts is multiple of 300s for unknown
% reason. I added the following  if statement to prevent that.
if mod(dts,300) == 0
    t=g.t1;
end

try  
   
    hh = floor(t/3600);
    mm = floor(mod((t/60), 60));
    ss = floor(mod(t,60));
    
        
    now = datenum(str2double(g.YYYY{:}),str2double(g.MM{:}),str2double(g.DD{:}),...
        hh,mm,ss);

   
catch
    errordlg('Date/Time input not valid; Check before continue','Plotter Error','modal')
    return
end

delta_t = ([00, 00, dts] * [3600; 60; 1]) / 86400;

new_t = datevec(now + delta_t);

year={num2str(new_t(1))};

if str2double(year) > 2025
    errordlg('Plotter2 only support up to the year 2025; Sumedhe should be done with his phd by then.','Year Error','modal')
    return
end


g.YYYY=year;
l=get(handles.year,'String');
g.YYYYn=find(ismember(l, g.YYYY)==1);


g.MM={sprintf('%2.2i',new_t(2))};
l=get(handles.month,'String');
g.MMn=find(ismember(l, g.MM)==1);

g.DD={sprintf('%2.2i',new_t(3))};
l=get(handles.date,'String');
g.DDn=find(ismember(l, g.DD)==1);

g.hh=new_t(4);
l=get(handles.hour,'String');
g.hhn=find(ismember(l, {sprintf('%2.2i',new_t(4))})==1);


g.mm=floor(new_t(5)/5)*5;
l=get(handles.minute,'String');
g.mmn=find(ismember(l, {sprintf('%2.2i',g.mm)})==1);
g.mmn=g.mmn(1);

g.ss=0;
l=get(handles.sec,'String');
g.ssn=find(ismember(l, {sprintf('%2.2i',g.ss)})==1);

t1minimum=3600*g.hh+60*g.mm+g.ss;
t2maximum=t1minimum+(5-mod(g.mm,5))*60;

g.t1min=t1minimum;
g.t2max=t2maximum;

g.t1 = g.t2;


if g.t1 <= t1minimum 
    g.t1 = t1minimum;
end

g.t2 = g.t1 + dts;

if g.t2 > t2maximum
    g.t2 = t2maximum;
end

handles.g=g;
guidata(hObject, handles);
set_values(handles)
plot_all_Callback(hObject, eventdata, handles)






% --- Executes on button press in plot_prev.
function plot_prev_Callback(hObject, eventdata, handles)

g=handles.g;

dts = handles.sen_set.next_prev_dt;

if g.t1 == g.t1min
    t = g.t1;
else
    t = g.t2;
end

try
   
   
    hh = floor(t/3600);
    mm = floor(mod((t/60), 60));
    ss = floor(mod(t,60));
    
        
    now = datenum(str2double(g.YYYY{:}),str2double(g.MM{:}),str2double(g.DD{:}),...
        hh,mm,ss);

   
catch
    errordlg('Date/Time input not valid; Check before continue','Plotter Error','modal')
    return
end

delta_t = ([00, 00, dts] * [3600; 60; 1]) / 86400;

new_t = datevec(now - delta_t);

year={num2str(new_t(1))};

if str2double(year) > 2025
    errordlg('Plotter2 only support up to the year 2025; Sumedhe should be done with his phd by then.','Year Error','modal')
    return
end


g.YYYY=year;
l=get(handles.year,'String');
g.YYYYn=find(ismember(l, g.YYYY)==1);


g.MM={sprintf('%2.2i',new_t(2))};
l=get(handles.month,'String');
g.MMn=find(ismember(l, g.MM)==1);

g.DD={sprintf('%2.2i',new_t(3))};
l=get(handles.date,'String');
g.DDn=find(ismember(l, g.DD)==1);

g.hh=new_t(4);
l=get(handles.hour,'String');
g.hhn=find(ismember(l, {sprintf('%2.2i',new_t(4))})==1);


g.mm=floor(new_t(5)/5)*5;
l=get(handles.minute,'String');
g.mmn=find(ismember(l, {sprintf('%2.2i',g.mm)})==1);
g.mmn=g.mmn(1);

g.ss=0;
l=get(handles.sec,'String');
g.ssn=find(ismember(l, {sprintf('%2.2i',g.ss)})==1);

t1minimum=3600*g.hh+60*g.mm+g.ss;
t2maximum=t1minimum+(5-mod(g.mm,5))*60;

g.t1min=t1minimum;
g.t2max=t2maximum;

g.t2 = g.t1;
g.t1 = g.t2 - dts;

if g.t2 > t2maximum || g.t2 < t1minimum
    g.t2 = t2maximum;
    g.t1 = g.t2 - dts;
end

if g.t1 < t1minimum;
    g.t1 = t1minimum;
end

handles.g=g;
guidata(hObject, handles);
set_values(handles)
plot_all_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function settings_Callback(hObject, eventdata, handles)
% hObject    handle to settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tools_Callback(hObject, eventdata, handles)
% hObject    handle to tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_37_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
g=handles.g;
try
    about_plotter(g.version,g.modi_date)
catch exception
    errordlg('Not available','About plottor','modal')
end

% --------------------------------------------------------------------
function usage_Callback(hObject, eventdata, handles)
dir=pwd;
fn='\html\plotter_help.html';

try
    open([dir fn])
catch execption
    errordlg('Usage not availble at this time. Please email Sumedhe:skarunar@olemiss.edu for bug fixes/questions','Plotter help','modal')
end

% --------------------------------------------------------------------
function handles = ldar_plots_Callback(hObject, eventdata, handles)

g=handles.g;

sen_set=handles.sen_set;

str={'07. north vs east (CCT)',...
    '01. 3D Color Coded by time',...
    '02. 3D Color Coded by type (CGLSS/DLS)',...
    '03. 3D Color Coded by altitue',...
    '04. 3D Color Coded by distance from a sensor',...
    '05. Altitude vs Time (CCT)',...
    '06. Altitude vs east (CCT)',...    
    '08. north vs altitude (CCT)',...
    '09. Altitude vs north (CCT)',...
    '10. All CCT graphs together',...
    '11. Altitude vs Time CC by type'...
    '12. Nort vs East CC by type'...
    '13. Histogram',...
    '14. Frequency Plot'...
    };


% Map quality
switch sen_set.map_quality 
    case 1; map_q_str = {'Low','Medium','High'};
    case 2; map_q_str = {'Medium','Low','High'};
    case 3; map_q_str = {'High','Low','Medium'};
end

% Find out available RADAR files according to current t1
[radarFiles, radarFullFiles] = getRadarFiles(handles);

[choise, button] = settingsdlg(...
    'Description', 'Select the date/time from the main Plotter window. Default t1 and t2 were also pulled out from the main windiw (CCT - Color Coded by TIme)',...
    'title' , 'LDAR - LINET - PBFA',...
    'Graph Type',str,...
    {'LDAR'; 'ldarOn'},logical(sen_set.ldarOn), ...
    {'CGLSS'; 'cglssOn'},logical(sen_set.cglssOn), ...
    {'LINET'; 'linetOn'},logical(sen_set.linetOn), ...
    {'PBFA'; 'pbfaOn'},logical(sen_set.pbfaOn), ...
    {'PBFA-A'; 'nldnOn'},logical(sen_set.nldnOn), ...
    {'PBFA-O'; 'pbfaOOn'},logical(sen_set.pbfaOOn), ...
    {'NLDN2'; 'nldn2On'},logical(sen_set.nldn2On), ...
    {'Plot Aditional Land Marks'; 'lmOn'},logical(sen_set.lmOn), ...
    {'Background map'; 'mapOn'},logical(sen_set.mapOn), ...
    {'Map quality';'map_quality'},map_q_str,...
    {'RADAR data'; 'radarOn'},logical(sen_set.radarOn),...
    {'RADAR file Name';'radarFn'},radarFiles,...
    {'RADAR elevation angle index';'radEleAngInd'},num2str(sen_set.radEleAngInd),...
    'WindowWidth' , 400,...
    'ControlWidth', 230);

if sum(choise.ldarOn + choise.cglssOn + choise.linetOn ...
        + choise.pbfaOn + choise.lmOn + choise.nldnOn + ...
        choise.nldn2On + choise.radarOn + choise.lmOn + ...
        choise.mapOn + choise.radarOn) == 0
    errordlg('No plots were selected to plot!','Plotter Error','modal')
    return
end

if strcmp(button,'ok')

    sen_set.lmOn=choise.lmOn;
    sen_set.ldarOn=choise.ldarOn;
    sen_set.cglssOn = choise.cglssOn;
    sen_set.linetOn=choise.linetOn;
    sen_set.pbfaOn=choise.pbfaOn;
    sen_set.nldnOn = choise.nldnOn;
    sen_set.pbfaOOn = choise.pbfaOOn;
    sen_set.mapOn = choise.mapOn;
    sen_set.nldn2On = choise.nldn2On;
    sen_set.ldar_graph_type = str2double(choise.GraphType(1:2));

    if isempty(radarFiles)
        sen_set.radarOn = choise.radarOn;
        sen_set.radarFn = '';
        sen_set.radEleAngInd = 1;
    else
        sen_set.radarOn = choise.radarOn;
        nn = find(strcmp(choise.radarFn,radarFiles),1);
        sen_set.radarFn = radarFullFiles{nn};
        sen_set.radEleAngInd = choise.radEleAngInd;
    end
    
    
    switch choise.map_quality
        case 'Medium' ; sen_set.map_quality = 2;
        case 'Low'    ; sen_set.map_quality = 1;
        case 'High'   ; sen_set.map_quality = 3;
    end
    
    save('sensor_setting.mat','-struct', 'sen_set')
    
    handles.g=g;
    handles.sen_set=sen_set;
    guidata(hObject, handles);
    
    ldar_plot_execute(handles)   
    
end


% --------------------------------------------------------------------
function fft_Callback(hObject, eventdata, handles)

sen_set=handles.sen_set;

[chanel, botton] = settingsdlg(...
    'Description', 'Choose Sensor ID and Channel Number and press OK to plot FFT graph. Date and Time Range selection Should be done from the main window',...
    'title' , 'Fast Fourier Transform',...
    'Sensor ID',sen_set.sen_IDs,...
    'Channel No',{'lp1', 'lp2','lp3','ch1','ch2','ch3'});

sn = find(strcmp(chanel.SensorID,sen_set.sen_IDs))*3 + ...
    str2double(chanel.ChannelNo(3))-3;

if strcmp(botton,'ok')
    
    g=handles.g;
    
    if strcmp(chanel.ChannelNo(1:2),'lp')
        
        
        % Generating the file name
        fn=sprintf('%s/%s/%s/%s/%s_%s%s%s_%2.2i%2.2i%2.2i.%s', ...
            sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},chanel.SensorID, ...
            g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,chanel.ChannelNo);
        
        % Check whether the file is exists
        if exist(fn,'file')~=0
            SA_FFT_01(fn,g.t1,g.t2,0)
        else
            str=sprintf('File %s is not found',fn);
            errordlg(str,'FFT plot error','modal')
        end
        
    else
        %Need to do ch FFT
        
        % finding the file name
        fn=sprintf('%s/%s/%s/%s/%s_%s%s%s_%2.2i%2.2i%2.2i.ch%s', ...
            sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},chanel.SensorID, ...
            g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,chanel.ChannelNo(1,3));
        
        hfn=sprintf('%s/%s/%s/%s/%s_%s%s%s_%2.2i%2.2i%2.2i.h%s', ...
            sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},chanel.SensorID, ...
            g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,chanel.ChannelNo(1,3));
        
        % Check whether the files are exists
        if exist(fn,'file') && exist(hfn,'file')
            FA_FFT_01(fn,hfn,g.t1,g.t2,sen_set.t_shift(sn),sen_set,sn)
        else
            str=sprintf('Files not found\n%s\n%s',fn,hfn);
            errordlg(str,'FFT plot error','modal')
        end
        
        
        
        
        
    end
        
end


% --------------------------------------------------------------------
function plot_all_hilbert_Callback(hObject, eventdata, handles)

g=handles.g;

%Is user forget to enter correct date?
if strcmp(g.YYYY,'yyyy') || strcmp(g.DD,'dd') || strcmp(g.MM,'mm')
    errordlg('Check Date input!','Date Error','modal')
    uiwait
    return
end


% Is user forget to enter correct time?
if isnan(g.hh) || isnan(g.mm) || isnan(g.ss)
    errordlg('Check Time input!','Time Error','modal')
    uiwait
    return
end

if g.t1 > g.t2 || isnan(g.t1) || isnan(g.t2)
    errordlg('Check Plotting Time Range!','Time range Error','modal')
    uiwait
    return
end



% Is user forget to click on any graphs?
if nnz(g.chgraphs)==0 && nnz(g.lpgraphs)==0 && g.ldar==0 && g.cglss==0 && g.linet==0 && g.pbfa ==0 && g.nldn==0
    errordlg('No graphs were selected to plot!','Graph selection Error','modal')
    uiwait
    return
end

% Lets setup calibration values before plot
handles=calibration(handles);

plot_all_hilbert2(handles.g)

guidata(hObject, handles);

% --------------------------------------------------------------------
function update_trange_Callback(hObject, eventdata, handles)
% hObject    handle to update_trange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function chpeaks_Callback(hObject, eventdata, handles)
% hObject    handle to chpeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function plot_all_color_coded_dots_Callback(hObject, eventdata, handles)
% hObject    handle to plot_all_color_coded_dots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function calibrate_Callback(hObject, eventdata, handles)
% hObject    handle to calibrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function calch_Callback(hObject, eventdata, handles)
    disp('Select exactly 3 ch graphs from plotter to continue')
    g=handles.g;
    calch(g)


% --------------------------------------------------------------------
function cal_lp_Callback(hObject, eventdata, handles)
        
    calibration4(handles.g)
    


% --------------------------------------------------------------------
function sensor_id_settings_Callback(hObject, eventdata, handles)
ids_positions
%uiwait
handles.sen_set=open('sensor_setting.mat');
guidata(hObject, handles);




% --------------------------------------------------------------------
function base_folder_Callback(hObject, eventdata, handles)

sen_set=handles.sen_set;

if sen_set.base_dir ==0
    sen_set.base_dir='None';
end
msg=sprintf('Choose the bease directory for slow and fast antenna data\n\nCurrent base dir = %s',...
    sen_set.base_dir);
sen_set.base_dir=uigetdir(sen_set.base_dir,msg);
save('sensor_setting.mat','-struct', 'sen_set');

handles.sen_set=sen_set;
guidata(hObject, handles);




% --------------------------------------------------------------------
function radius_Callback(hObject, eventdata, handles)

sen_set=handles.sen_set;

prompt =sprintf('Effective Radius (in meters) for \n    LDAR    LINET\n    PBFA    NLDN');
dlg_title = 'Radius';
num_lines = 1;
def = sen_set.ldar_r;
answer = inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(answer)
    sen_set.ldar_r=answer;
    save('sensor_setting.mat','-struct', 'sen_set')
    handles.sen_set=sen_set;
    guidata(hObject, handles);
end




% --------------------------------------------------------------------
function calibration_para_Callback(hObject, eventdata, handles)

a=handles.sen_set;

str={'01. Use daily average',...
     '02. Use given time range average' ...
     '03. Manual calibration'
     };

%Discription
discrip =sprintf('Select the type of calibration you want\n       Currently Using: "%s"',str{a.cal_type});


[choise, button] = settingsdlg(...
    'Description', discrip ,...
    'title' , 'Calibration Parameters',...
    {'Plot Calibrated values'; 'plot_calibrated'},[logical(a.plot_calibrated) ~logical(a.plot_calibrated) ], ...
    'Calibration Type',str, ...
    {'Cal Range Start (yyyymmdd)', 'cal_start'}, num2str(a.cal_start), ...
    {'Cal Range End (yyyymmdd)', 'cal_end'}, num2str(a.cal_end), ...    
     'WindowWidth' , 400,...
    'ControlWidth', 230);

if strcmp(button,'ok')
    a.cal_type=str2double(choise.CalibrationType(1:2));
    a.cal_start = choise.cal_start;
    a.cal_end = choise.cal_end;
    a.plot_calibrated = choise.plot_calibrated;
    save('sensor_setting.mat','-struct','a')
    handles.sen_set=a;
    handles=calibration(handles);    
    guidata(hObject, handles);
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% --------------------------------------------------------------------
function clear_all_Callback(hObject, eventdata, handles)
g=handles.g;

g.YYYYn=1;
g.YYYY = {'YYYY'};
g.MMn =1;
g.MM = {'MM'};
g.DDn=1;
g.DD = {'DD'};
g.hhn=1;
g.hh = {'hh'};
g.mmn=1;
g.mm = {'mm'};
g.ssn=1;
g.ss ={'ss'};
g.fn=0;
g.pn=0;
g.save_fn=0;
g.save_pn=0;
g.lpgraphs=zeros(1,60);
g.chgraphs=zeros(1,60);
g.linet=0;
g.ldar=0;
g.cglss = 0;
g.pbfa =0 ;
g.nldn =0;
g.t1=0;
g.t2=0;
g.t1min=0;
g.t2max=0;
g.isSaved=1;

handles.g = g;
guidata(hObject, handles);
set_values(handles)




% --------------------------------------------------------------------
function clear_gs_Callback(hObject, eventdata, handles)

g=handles.g;

g.lpgraphs=zeros(1,60);
g.chgraphs=zeros(1,60);
g.linet = 0;
g.ldar =0 ;
g.cglss = 0;
g.pbfa = 0;
g.nldn = 0;

handles.g = g;
guidata(hObject, handles);
set_values(handles)



% --------------------------------------------------------------------
function clear_date_Callback(hObject, eventdata, handles)

g =  handles.g;

g.YYYYn=1;
g.YYYY = {'YYYY'};
g.MMn =1;
g.MM = {'MM'};
g.DDn=1;
g.DD = {'DD'};

handles.g=g;
guidata(hObject, handles);
set_values(handles)


% --------------------------------------------------------------------
function clear_time_Callback(hObject, eventdata, handles)

g = handles.g;

g.hhn=1;
g.hh = {'hh'};
g.mmn=1;
g.mm = {'mm'};
g.ssn=1;
g.ss ={'ss'};


g.t1min = 0;
g.t2max = 0;
g.t1 = 0;
g.t2 = 0;

handles.g = g;
guidata(hObject, handles);
set_values(handles)


% --------------------------------------------------------------------
function new_Callback(hObject, eventdata, handles)
% hObject    handle to new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function open_Callback(hObject, eventdata, handles)

% Default location
def_loc = [handles.sen_set.working_dir '/*.mat'];

[fn,pn] = uigetfile(def_loc);

filename=[pn fn];

if fn==0
    % File did not selected, do nothing!
else
    g=open(filename);
    
    % this will be necessary if user rename the file
    g.save_fn = fn;
    g.save_pn = pn;
    
    g.isSaved = 1;  
    
    % Change the working directory to this directory
    handles.sen_set.working_dir = pn;
        
    % If user open an old version, let's fill the blanks and make it
    % compatible with version 2.0
    if strcmp(g.version, '1.0')
        
        % Remove not obsalete variables for version 2
        g=rmfield(g,{'note','chfft','lpfft','chfilter','lpfilter'});
        
        % Add/modify variables version 2 
        g.lpgraphs=[g.lpgraphs zeros(1,45)];
        g.YYYYn=g.YYYYn+2; % This is necessary as plooter started from 2009 instead of 2007 (my mistake)
        if g.mmn > 14; g.mmn=g.mmn+1; end;
        g.chgraphs=[g.chgraphs zeros(1,45)];
        g.nldn = 0;
        
        g.dedt = 0;
        g.dedt_trig = 0;
        g.dedt_pps = 0;
        g.dedt_opt = 0;
        
        % One of the plotter vertion didn't have pbfa, lets add it
        if ~isfield(g,'pbfa')
            g.pbfa = 0;
        end
  
        % nldn2 started with version 3
        g.nldn = 0;
        
        % Make it version 3 if user want to save 
        g.version= '3.0'; 
        g.modi_date='2013-July-06';       
        
    elseif strcmp(g.version, '2.0')
        % nldn2 started with version number 3
        g.nldn2 = 0;
        
        % Make it version 3 if user want to save 
        g.version= '3.0'; 
        g.modi_date='2013-July-06';  
        
        
    end
    
    % I added dedt later. So lets make sure user can open saved plotter2
    % files before I started these modification
    try 
        if g.de/dt == 0
            % Do nothing
        end
    catch
        g.dedt = 0;
        g.dedt_trig = 0;
        g.dedt_pps = 0;
        g.dedt_opt = 0;
    end
    
    % I added CGLSS later. So lets make sure user can open saved plotter2
    % files before I started these modification
    try 
        if g.cglss == 0
            % Do nothing
        end
    catch
        g.cglss = g.ldar;
    end
    try
        if g.temp.ch_hpf==0
        end
    catch
         g.temp.ch_hpf = 0;
    end

    % I added KSC feild mill plotting later, need to make sure we can plot
    % them
    try
        if g.KSC_FM(1)==0
        end
    catch
         g.KSC_FM = zeros(1,34);
    end
    
    handles.g=g;
    
    guidata(hObject, handles);
    
    set_values(handles)
    
end




% --------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)

g=handles.g;
g.sen_set = handles.sen_set;

if g.save_fn(1)==0
    % suggested filename
    fn=sprintf('%s%s%s_%5.5us_%3.3ums_',...
        g.YYYY{:},g.MM{:},g.DD{:},...
        floor(g.t1),floor((g.t1-floor(g.t1))*1000));
    [g.save_fn,g.save_pn] = uiputfile('*.mat','Save..',fn);
    if g.save_fn(1)~=0
        fn=[g.save_pn g.save_fn];
        save(fn,'-struct', 'g')
        set(handles.plotter2,'Name',['Plotter -' g.save_fn])
        g.isSaved=1;
        handles.sen_set.working_dir = g.save_pn;
    end
else
    fn=[g.save_pn g.save_fn];
    save(fn,'-struct', 'g')
    g.isSaved=1;
end
handles.g=g;
guidata(hObject, handles);


% --------------------------------------------------------------------
function saveas_Callback(hObject, eventdata, handles)

g=handles.g;
g.sen_set = handles.sen_set;

file=sprintf('%s//%s%s%s_%5.5us_%3.3ums_',handles.sen_set.working_dir,...
    g.YYYY{:},g.MM{:},g.DD{:},floor(g.t1),floor((g.t1-floor(g.t1))*1000));
[fn,pn] = uiputfile('*.mat','Save As..',file);
if fn(1)~=0
    g.save_fn=fn;
    g.save_pn=pn;
    fn=[pn fn];
    save(fn,'-struct', 'g')
    set(handles.plotter2,'Name',['Plotter2 -' g.save_fn])
    g.isSaved=1;
    handles.sen_set.working_dir = g.save_pn;
    handles.g=g;
    guidata(hObject, handles);
end

function close_Callback(hObject, eventdata, handles)

    plotter2_CloseRequestFcn(hObject, eventdata, handles)



% --------------------------------------------------------------------
function print_Callback(hObject, eventdata, handles)

printdlg('-crossplatform')


% --------------------------------------------------------------------
function new_plotter_Callback(hObject, eventdata, handles)
%try 
    %close_Callback(hObject,eventdata,handles)
    close(handles.plotter2)
    pause(1)
    plotter2(1)
    
%catch exception
 %   errordlg('Cannot Open a new Plotter!','Plotter error','modal')
%end

% --- Executes when user attempts to close plotter2.
function plotter2_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to plotter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%disp('Working000')
g=handles.g;

if g.isSaved==0
    %disp('Working000')
    answer=questdlg('Save before exit2?','Exit','Yes','No','Yes');
    if strcmp(answer, 'Yes')
        save_Callback(hObject,eventdata,handles)
        %display('Working0')
    elseif strcmp(answer,'No')
        % Save tempary file befor quite
        try
            save('plotter_last_gui_data.mat','-struct', 'g')
        end
        %display('Working1')
        delete(handles.plotter2)
        %delete(hObject)
    else
        % User undisided. Do nothing
    end
else
    % Save tempary file befor quite
    %try
    save('plotter_last_gui_data.mat','-struct', 'g')
    %display('Working2')
    delete(handles.plotter2)
    clear global g
    %delete(hObject)
end

function set_values(handles)

%disp('Working set_values')

a=handles.g;


set(handles.year,'Value',a.YYYYn)
set(handles.month,'Value',a.MMn)
set(handles.date, 'Value',a.DDn)

set(handles.hour,'Value',a.hhn)
set(handles.minute,'Value',a.mmn)
set(handles.sec,'Value',a.ssn)
%set(handles.base_dir,'Value',a.Sidn)

set(handles.t1,'String',num2str(a.t1))
set(handles.t2,'String',num2str(a.t2))

t1min_s=sprintf('t Min = %is',a.t1min);
t2max_s=sprintf('t Max = %is',a.t2max);
set(handles.t1min,'String',t1min_s)
set(handles.t2max,'String',t2max_s)

set(handles.lp011,'Value',a.lpgraphs(1))
set(handles.lp012,'Value',a.lpgraphs(2))
set(handles.lp013,'Value',a.lpgraphs(3))
set(handles.lp021,'Value',a.lpgraphs(4))
set(handles.lp022,'Value',a.lpgraphs(5))
set(handles.lp023,'Value',a.lpgraphs(6))
set(handles.lp031,'Value',a.lpgraphs(7))
set(handles.lp032,'Value',a.lpgraphs(8))
set(handles.lp033,'Value',a.lpgraphs(9))
set(handles.lp041,'Value',a.lpgraphs(10))
set(handles.lp042,'Value',a.lpgraphs(11))
set(handles.lp043,'Value',a.lpgraphs(12))
set(handles.lp051,'Value',a.lpgraphs(13))
set(handles.lp052,'Value',a.lpgraphs(14))
set(handles.lp053,'Value',a.lpgraphs(15))
set(handles.lp061,'Value',a.lpgraphs(16))
set(handles.lp062,'Value',a.lpgraphs(17))
set(handles.lp063,'Value',a.lpgraphs(18))
set(handles.lp071,'Value',a.lpgraphs(19))
set(handles.lp072,'Value',a.lpgraphs(20))
set(handles.lp073,'Value',a.lpgraphs(21))
set(handles.lp081,'Value',a.lpgraphs(22))
set(handles.lp082,'Value',a.lpgraphs(23))
set(handles.lp083,'Value',a.lpgraphs(24))
set(handles.lp091,'Value',a.lpgraphs(25))
set(handles.lp092,'Value',a.lpgraphs(26))
set(handles.lp093,'Value',a.lpgraphs(27))
set(handles.lp101,'Value',a.lpgraphs(28))
set(handles.lp102,'Value',a.lpgraphs(29))
set(handles.lp103,'Value',a.lpgraphs(30))
set(handles.lp111,'Value',a.lpgraphs(31))
set(handles.lp112,'Value',a.lpgraphs(32))
set(handles.lp113,'Value',a.lpgraphs(33))
set(handles.lp121,'Value',a.lpgraphs(34))
set(handles.lp122,'Value',a.lpgraphs(35))
set(handles.lp123,'Value',a.lpgraphs(36))
set(handles.lp131,'Value',a.lpgraphs(37))
set(handles.lp132,'Value',a.lpgraphs(38))
set(handles.lp133,'Value',a.lpgraphs(39))
set(handles.lp141,'Value',a.lpgraphs(40))
set(handles.lp142,'Value',a.lpgraphs(41))
set(handles.lp143,'Value',a.lpgraphs(42))
set(handles.lp151,'Value',a.lpgraphs(43))
set(handles.lp152,'Value',a.lpgraphs(44))
set(handles.lp153,'Value',a.lpgraphs(45))
set(handles.lp161,'Value',a.lpgraphs(46))
set(handles.lp162,'Value',a.lpgraphs(47))
set(handles.lp163,'Value',a.lpgraphs(48))
set(handles.lp171,'Value',a.lpgraphs(49))
set(handles.lp172,'Value',a.lpgraphs(50))
set(handles.lp173,'Value',a.lpgraphs(51))
set(handles.lp181,'Value',a.lpgraphs(52))
set(handles.lp182,'Value',a.lpgraphs(53))
set(handles.lp183,'Value',a.lpgraphs(54))
set(handles.lp191,'Value',a.lpgraphs(55))
set(handles.lp192,'Value',a.lpgraphs(56))
set(handles.lp193,'Value',a.lpgraphs(57))
set(handles.lp201,'Value',a.lpgraphs(58))
set(handles.lp202,'Value',a.lpgraphs(59))
set(handles.lp203,'Value',a.lpgraphs(60))


set(handles.ch011,'Value',a.chgraphs(1))
set(handles.ch012,'Value',a.chgraphs(2))
set(handles.ch013,'Value',a.chgraphs(3))
set(handles.ch021,'Value',a.chgraphs(4))
set(handles.ch022,'Value',a.chgraphs(5))
set(handles.ch023,'Value',a.chgraphs(6))
set(handles.ch031,'Value',a.chgraphs(7))
set(handles.ch032,'Value',a.chgraphs(8))
set(handles.ch033,'Value',a.chgraphs(9))
set(handles.ch041,'Value',a.chgraphs(10))
set(handles.ch042,'Value',a.chgraphs(11))
set(handles.ch043,'Value',a.chgraphs(12))
set(handles.ch051,'Value',a.chgraphs(13))
set(handles.ch052,'Value',a.chgraphs(14))
set(handles.ch053,'Value',a.chgraphs(15))
set(handles.ch061,'Value',a.chgraphs(16))
set(handles.ch062,'Value',a.chgraphs(17))
set(handles.ch063,'Value',a.chgraphs(18))
set(handles.ch071,'Value',a.chgraphs(19))
set(handles.ch072,'Value',a.chgraphs(20))
set(handles.ch073,'Value',a.chgraphs(21))
set(handles.ch081,'Value',a.chgraphs(22))
set(handles.ch082,'Value',a.chgraphs(23))
set(handles.ch083,'Value',a.chgraphs(24))
set(handles.ch091,'Value',a.chgraphs(25))
set(handles.ch092,'Value',a.chgraphs(26))
set(handles.ch093,'Value',a.chgraphs(27))
set(handles.ch101,'Value',a.chgraphs(28))
set(handles.ch102,'Value',a.chgraphs(29))
set(handles.ch103,'Value',a.chgraphs(30))
set(handles.ch111,'Value',a.chgraphs(31))
set(handles.ch112,'Value',a.chgraphs(32))
set(handles.ch113,'Value',a.chgraphs(33))
set(handles.ch121,'Value',a.chgraphs(34))
set(handles.ch122,'Value',a.chgraphs(35))
set(handles.ch123,'Value',a.chgraphs(36))
set(handles.ch131,'Value',a.chgraphs(37))
set(handles.ch132,'Value',a.chgraphs(38))
set(handles.ch133,'Value',a.chgraphs(39))
set(handles.ch141,'Value',a.chgraphs(40))
set(handles.ch142,'Value',a.chgraphs(41))
set(handles.ch143,'Value',a.chgraphs(42))
set(handles.ch151,'Value',a.chgraphs(43))
set(handles.ch152,'Value',a.chgraphs(44))
set(handles.ch153,'Value',a.chgraphs(45))
set(handles.ch161,'Value',a.chgraphs(46))
set(handles.ch162,'Value',a.chgraphs(47))
set(handles.ch163,'Value',a.chgraphs(48))
set(handles.ch171,'Value',a.chgraphs(49))
set(handles.ch172,'Value',a.chgraphs(50))
set(handles.ch173,'Value',a.chgraphs(51))
set(handles.ch181,'Value',a.chgraphs(52))
set(handles.ch182,'Value',a.chgraphs(53))
set(handles.ch183,'Value',a.chgraphs(54))
set(handles.ch191,'Value',a.chgraphs(55))
set(handles.ch192,'Value',a.chgraphs(56))
set(handles.ch193,'Value',a.chgraphs(57))
set(handles.ch201,'Value',a.chgraphs(58))
set(handles.ch202,'Value',a.chgraphs(59))
set(handles.ch203,'Value',a.chgraphs(60))

set(handles.ldar,'Value',a.ldar)
set(handles.cglss,'Value',a.cglss)
set(handles.linet,'Value',a.linet)
set(handles.pbfa,'Value',a.pbfa)
set(handles.nldn,'Value',a.nldn)
set(handles.nldn2,'Value',a.nldn2)

set(handles.dedt,'Value',a.dedt)
set(handles.dedt_trig,'Value',a.dedt_trig)
set(handles.dedt_pps,'Value',a.dedt_pps)
set(handles.dedt_opt,'Value',a.dedt_opt)

if isfield('a','pbfa_old')
    a.pbfa_old = 0;
end
set(handles.pbfa_old,'Value',a.pbfa_old);
    

if handles.sen_set.pbfa_ebar_on
    set(handles.pbfa_error_bar_on,'Checked','on')
else
    set(handles.pbfa_error_bar_on,'Checked','off')
end

if a.save_fn~=0
    set(handles.plotter2,'Name',['Plotter2 -' a.save_fn])
else
    set(handles.plotter2,'Name','Plotter2')
end

set(handles.auto_b_folder,'Checked',handles.sen_set.autoBaseFolder)

function set_sensor_ids(handles)
set(handles.sen01,'String',handles.sen_set.sen_IDs{1})
set(handles.sen02,'String',handles.sen_set.sen_IDs{2})
set(handles.sen03,'String',handles.sen_set.sen_IDs{3})
set(handles.sen04,'String',handles.sen_set.sen_IDs{4})
set(handles.sen05,'String',handles.sen_set.sen_IDs{5})
set(handles.sen06,'String',handles.sen_set.sen_IDs{6})
set(handles.sen07,'String',handles.sen_set.sen_IDs{7})
set(handles.sen08,'String',handles.sen_set.sen_IDs{8})
set(handles.sen09,'String',handles.sen_set.sen_IDs{9})
set(handles.sen10,'String',handles.sen_set.sen_IDs{10})
set(handles.sen11,'String',handles.sen_set.sen_IDs{11})
set(handles.sen12,'String',handles.sen_set.sen_IDs{12})
set(handles.sen13,'String',handles.sen_set.sen_IDs{13})
set(handles.sen14,'String',handles.sen_set.sen_IDs{14})
set(handles.sen15,'String',handles.sen_set.sen_IDs{15})
set(handles.sen16,'String',handles.sen_set.sen_IDs{16})
set(handles.sen17,'String',handles.sen_set.sen_IDs{17})
set(handles.sen18,'String',handles.sen_set.sen_IDs{18})
set(handles.sen19,'String',handles.sen_set.sen_IDs{19})
set(handles.sen20,'String',handles.sen_set.sen_IDs{20})
set(handles.dEdt_sn,'String',handles.sen_set.dEdt.sn)



% --- Executes on button press in close_all_graphs.
function close_all_graphs_Callback(hObject, eventdata, handles)

choise = handles.sen_set.button_fun(1);

    switch choise
        
        case 1; ldar_plot_execute(handles);
        case 2; calibration4(handles.g);
        case 3; calch(handles.g);
        case 4; n_CGLSS_finder3(handles);
        case 5; close(gcf)
        case 6; close all
        case 7; show_cglss_info_Callback(hObject, eventdata, handles);
        case 8; high_z_ldar(handles)
        case 9; distance2ldar
        case 10; find_max_dE(handles);
        case 11; toggle_inputs_Callback(hObject, eventdata, handles);
            
        otherwise
            disp('Error!')
    end

% --- Executes on button press in close_current_graph.
function close_current_graph_Callback(hObject, eventdata, handles)

   choise = handles.sen_set.button_fun(2);

    switch choise
        
        case 1; ldar_plot_execute(handles);
        case 2; calibration4(handles.g);
        case 3; calch(handles.g);
        case 4; n_CGLSS_finder3(handles);
        case 5; close(gcf);
        case 6; close all
        case 7; show_cglss_info_Callback(hObject, eventdata, handles);
        case 8; high_z_ldar(handles)
        case 9; distance2ldar
        case 10; find_max_dE(handles);
        case 11; toggle_inputs_Callback(hObject, eventdata, handles);
            
        otherwise
            disp('Error!')
    end


% --------------------------------------------------------------------
function ldar_tshift_Callback(hObject, eventdata, handles)
sen_set=handles.sen_set;

str2 = {'LDAR' 'PBFA' 'LINET' 'PBFA-A' 'PBFA-O'};

if sen_set.ldar_tshiftOn
    str=sprintf('LDAR LINET PBFA time shift currently sat to %s',sen_set.sen_IDs{sen_set.ldar_tshift_sn});
elseif sen_set.plot_tshiftOn
    str=sprintf('Currently plots are time shifted to %s',str2{sen_set.plot_tshiftType});
else
    str='';
end

[tshifts button] = settingsdlg(...
'description',str, ...
'title' , 'Time shift Settings',...
'separator' , 'Time Shift LDAR to',...
{'Turn ON'; 'ldar_tshift'}, boolean(sen_set.ldar_tshiftOn),...
{'Sensor ID';'sid'}, sen_set.sen_IDs,...
'separator' , 'Time shift plots to',...
{'Turn ON'; 'plot_tshiftOn'}, boolean(sen_set.plot_tshiftOn),...
{'Type';'psid'}, str2);


if strcmp(button,'ok')
    % Let's see user selected both check marks
    if tshifts.ldar_tshift && tshifts.plot_tshiftOn
        % Let user know only one choice can be maid
        estr = sprintf('Only one time shift types allowed! \nLDAR will be time shifted to %s.',tshifts.sid);
        errordlg(estr,'Time shift error!')
        tshifts.plot_tshiftOn = 0;
    end
    
    sen_set.ldar_tshiftOn=tshifts.ldar_tshift;
    sen_set.ldar_tshift_sn=find(ismember(sen_set.sen_IDs, tshifts.sid)==1);
    sen_set.plot_tshiftOn = tshifts.plot_tshiftOn;
    sen_set.plot_tshiftType = find(ismember(str2, tshifts.psid)==1);
    
    save('sensor_setting.mat','-struct','sen_set')
    handles.sen_set=sen_set;
    guidata(hObject, handles);
end




% --------------------------------------------------------------------
function ch_freq_Callback(hObject, eventdata, handles)

a=handles.sen_set;

method={'Summing','Down Sampling'};
freq={'5.00MHz','10.00MHz','5.00MHz','2.50MHz','1.00MHz','0.50MHz','0.25MHz','0.10MHz','0.05MHz','0.01MHz','0.005MHz','0.001MHz'};

if a.chSumMethod ~= 0
    str = sprintf('choose the meothd (Current = %s) ',...
        method{a.chSumMethod});
    yes = [true false];
    
else
    str = sprintf('choose the meothd and the freqency. \n\n        Current method    = None ' );
    yes = [false true];
end

[choise button] = settingsdlg(...
'description',str, ...
'title' , 'Ch Data Freqency',...
{'Include Undersampling'; 'yes'}, yes,...
{'Method'; 'method'}, method,...
{'Frequency';'freq'}, freq,...
{'ch1'; 'ch1'}, false,...
{'ch2'; 'ch2'}, false,...
{'ch3'; 'ch3'}, false);

if strcmp(button,'ok')
    if choise.yes 
        if strcmp(choise.method ,'Summing')
            a.chSumMethod = 1;
        else
            a.chSumMethod = 2;
        end
    else
        a.chSumMethod = 0;
    end
    
    chfreq = str2double(choise.freq(1:end-3))*1e6;
    
    if choise.ch1; for i = 1:20; a.chfreq(i*3-2) = chfreq;    end;    end
    if choise.ch2; for i = 1:20; a.chfreq(i*3-1) = chfreq;    end;    end
    if choise.ch3; for i = 1:20; a.chfreq(i*3)   = chfreq;    end;    end
    
    handles.sen_set=a;
    save('sensor_setting.mat','-struct','a')
    guidata(hObject, handles);       
end


% --- Executes on button press in dedt.
function dedt_Callback(hObject, eventdata, handles)
handles.g.dedt=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in dedt_trig.
function dedt_trig_Callback(hObject, eventdata, handles)
handles.g.dedt_trig=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --- Executes on button press in dedt_opt.
function dedt_opt_Callback(hObject, eventdata, handles)
handles.g.dedt_opt=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);

% --- Executes on button press in dedt_pps.
function dedt_pps_Callback(hObject, eventdata, handles)
handles.g.dedt_pps=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --------------------------------------------------------------------
function dEdt_settings_Callback(hObject, eventdata, handles)
dedt_id_position
uiwait
handles.sen_set=open('sensor_setting.mat');
guidata(hObject, handles);


% --------------------------------------------------------------------
function cal_pbfa_Callback(hObject, eventdata, handles)
    peak_modifier


% --------------------------------------------------------------------
function next_btn_Callback(hObject, eventdata, handles)
   % configure the next button time
   prompt='dt for ''Next'' ''Prev'' buttons in seconds:';
   name='Delta t';
   numlines=1;   
   dt = handles.sen_set.next_prev_dt;   
   defaultanswer={num2str(dt)};
   answer=inputdlg(prompt,name,numlines,defaultanswer);   
   
   if ~isempty(answer)       
       handles.sen_set.next_prev_dt = str2double(answer);
       set(handles.plot_next,'String',['Next ' answer{:} 's'])
       set(handles.plot_prev,'String',['Prev ' answer{:} 's'])
       guidata(hObject, handles);       
   end
   
  

   


% --------------------------------------------------------------------
function show_cglss_info_Callback(hObject, eventdata, handles)
    g = handles.g;
    settings = handles.sen_set;
    
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
    
    % Exit if the file does not exsist
    if ~exist(ldar_fn,'file')
        errordlg(ldar_fn,'File not found')
        return
    end
    
    %% Time shift for LDAR and Linet
    if settings.ldar_tshiftOn==1
        sn=settings.ldar_tshift_sn;
        x0=settings.x(sn)-settings.x0;
        y0=settings.x(sn)-settings.y0;
        z0=settings.x(sn)-settings.z0;
        name = settings.sen_IDs{sn};
        
    else
        x0=0;
        y0=0;
        z0=0;
        name = 'None';
    end
    
    % Load LDAR data
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0,0);
    
    CG(:,11) = int64(CG(:,11));
    
    if ~isnan(CG(1,1))
        fprintf('\nOccured(s) \t\t\tx(m) \t\ty(m) \t\tD(m) \t\tDelay to %s (us) \n',name)
        for i=1:length(CG(:,1))
            fprintf('%6.7f \t%9d\t%9d \t%9d \t\t%.1f \n',...
                CG(i,9),CG(i,6),CG(i,7),CG(i,11),(CG(i,10)-CG(i,9))*1e6)
        end
        
        clipboard('copy', CG(1,11));
    end
    
    
    


% --------------------------------------------------------------------
function find_high_z_ldar_Callback(hObject, eventdata, handles)
%% 
g = handles.g;
settings = handles.sen_set;

try
    max_z = handles.temp.max_z;
catch
    max_z = 15000;
end

% Get user altitude
answer = inputdlg('Enter the Altitude (m)','High z LDAR2',1,{num2str(max_z)});

if isempty(answer)
    return
end

max_z = str2double(answer{1});

% Update this to handles as the next default
handles.temp.max_z = max_z;
guidata(hObject, handles);

high_z_ldar(handles)



% --- Executes on button press in cglss.
function cglss_Callback(hObject, eventdata, handles)
handles.g.cglss=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --------------------------------------------------------------------
function n_stroke_finder_Callback(hObject, eventdata, handles)
%% This function will find n strokes...
% clc
str = 'This tool will find n-strokes flashes.';
str = [str ' Number of strokes will be found withn time threshold T_th'];

% load previous values if there are..
try
    N = handles.temp.N;
    delT = handles.temp.delT;
catch
    N = 1;
    delT =1;
end

[choise button] = settingsdlg(...
    'description',str, ...
    'title' , 'Ch Data Freqency',...
    {'LINET'; 'linet'}, true,...
    {'CGLSS'; 'cglss'}, true,...
    {'Number of strokes'; 'N'},N ,...
    {'Time threshold (s)';'delT'},delT );

% If user didn't press OK, let's exit
if strcmp(button,'ok')
    % lets save values to handle
    handles.temp.N = choise.N;
    handles.temp.delT = choise.delT;
    handles.temp.cglss = choise.cglss;
    handles.temp.linet = choise.linet;
    guidata(hObject, handles);
    n_CGLSS_finder3(handles);
end


% --------------------------------------------------------------------
function customise_LR_buttons_Callback(hObject, eventdata, handles)

str={'01. LDAR plots',...
     '02. LP calibration',...
     '03. Ch Calibration',...
     '04. N-stoke finder',...
     '05. Close current plot',...
     '06. Close All Open plots',...
     '07. Show CGLSS info',...
     '08. Find high altitude LDAR',...
     '09. Same distance CGLSS', ....
     '10. Find max dE',....
     '11. Toggle Inputs', ...
    };

[choise, button] = settingsdlg(...
    'Description', 'Customoise Left Most and Right Most buttons to do special tasks.',...
    'title' , 'Customize Buttons',...
    'Button',{'01. Left Most','02. Right Most'},...
    'Action',str,...
    'WindowWidth' , 400,...
    'ControlWidth', 230);

if ~strcmp(button,'ok')
    return
end


% sen_set.button_fn array holds the info about button customization
% index 1 : Left most button
% index 2 : Right most button
sen_set = handles.sen_set;
ind = str2double(choise.Button(1:2));
sen_set.button_fun(ind) = str2double(choise.Action(1:2));


handles.sen_set = sen_set;
setup_buttons(handles);

save('sensor_setting.mat','-Struct','sen_set')

guidata(hObject, handles);

function setup_buttons(handles)

button_str = {'LDAR','LP cal','Ch cal','N-stroks','X Curr','X All',...
                'CGLSS info','Hi Z LDAR','=R LDAR','Max dE','Toggle'};

set(handles.close_all_graphs,'String',...
    button_str{handles.sen_set.button_fun(1)})

set(handles.close_current_graph,'String',...
    button_str{handles.sen_set.button_fun(2)})

function ldar_plot_execute(handles)

if strcmp(handles.sen_set.autoBaseFolder,'on')
    handles = set_auto_b_folder(handles);
end

g = handles.g;
sen_set = handles.sen_set;

% set the base folder


if g.mm < 30
    ext=0;
else
    ext=30;
end

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);


linet_fn=sprintf('%s/LINET/%s/%s/%s/linet_%s%s%s_%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

% nldn_fn=sprintf('%s/LINET2/%s/%s/%s/linet_%s%s%s_%2.2i%2.2i.txt',...
%     sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
%     g.MM{:},g.DD{:},g.hh,ext);

nldn_fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

%Old PBFA fn
pbfaO_fn=sprintf('%s/PBFA_old/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

nldn2_fn=sprintf('%s/NLDN2/%s/%s/NLDN2_%s%s%s.txt',...
        sen_set.loc_dir,g.YYYY{:},g.MM{:},g.YYYY{:},...
        g.MM{:},g.DD{:});
    
% Graph selection from the user
num=sen_set.ldar_graph_type;

% Do we need ldar time shift to be included?
if sen_set.ldar_tshiftOn==1
    x0=sen_set.x(sen_set.ldar_tshift_sn)-sen_set.x0;
    y0=sen_set.y(sen_set.ldar_tshift_sn)-sen_set.y0;
    z0=sen_set.z(sen_set.ldar_tshift_sn)-sen_set.z0;
else
    x0=0;
    y0=0;
    z0=0;
end

switch num
    case 1
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,sen_set,g)
    case 2
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,2,sen_set,g)
    case 3
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,3,sen_set,g)
    case 4
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,4,sen_set,g)
    case 5
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[1,0,0,0,0],sen_set,g);
    case 6
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[0,1,0,0,0],sen_set);
    case 7
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[0,0,1,0,0],sen_set,g);
    case 8
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[0,0,0,1,0],sen_set,g);
    case 9
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[0,0,0,0,1],sen_set,g);
    case 10
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,0,[0,1,0,0,0],sen_set,g);
    case 11
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,5,sen_set,g)
    case 12
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,6,sen_set,g)
    case 13
        LDAR_Hist1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,sen_set)
    case 14
        LDAR_frequency(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,sen_set)
        
    otherwise
        disp xc
        
end


function high_z_ldar(handles)

g = handles.g;
max_z = handles.temp.max_z;
settings = handles.sen_set;

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

% Exit if the file does not exsist
if ~exist(ldar_fn,'file')
    errordlg(ldar_fn,'File not found')
    return
end

% Time shift for LDAR and Linet
if settings.ldar_tshiftOn==1
    sn=settings.ldar_tshift_sn;
    x0=settings.x(sn)-settings.x0;
    y0=settings.x(sn)-settings.y0;
    z0=settings.x(sn)-settings.z0;
    name = settings.sen_IDs{sn};
    
else
    x0=0;
    y0=0;
    z0=0;
    name = 'None';
end

% Load LDAR data
[CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
    x0,y0,z0,0);

%DLS(:,11) = int64(DLS(:,8));

% Sort according to altitude
DLS = sortrows(DLS,8);

% Choose only heigher altitudes
lol=length(DLS(:,8))- nnz(DLS(:,8) > max_z)+1 ;      % Index of the lower matrix element
ul =length(DLS(:,8))  ;                         % Index of the upper matrix element  

DLS = DLS(lol:ul,:);

% Sort back according to time
DLS = sortrows(DLS,9);

if isempty(DLS) || isnan(DLS(1,1))
    fprintf('\n\nNo points found higher than %.1f m \n',max_z)
    return
end

fprintf('\nPrinting LDAR2 points of Altitudes > %.1fm',max_z) 
fprintf('\nOccured(s) \t\t\tx(m) \t\ty(m) \t\tZ(m) \t\tDelay to %s (us) \n',name)
for i=1:length(DLS(:,1))
    fprintf('%6.7f \t%9d\t%9d \t%9d \t\t%.1f \n',...
        DLS(i,9),DLS(i,6),DLS(i,7),DLS(i,8),(DLS(i,10)-DLS(i,9))*1e6)
end


% --------------------------------------------------------------------
function find_max_dE_Callback(hObject, eventdata, handles)

find_max_dE(handles);


% --------------------------------------------------------------------
function toggle_inputs_Callback(hObject, eventdata, handles)

try
    handles.g       = handles.temp.g1;
    handles.temp.g1 = handles.temp.g2;
    handles.temp.g2 = handles.g;
    guidata(hObject, handles);    
    set_values(handles)
    
catch
    disp('Data not found. Sorry!')
end


function handles = tshift_warning(handles)

if handles.tshift_war == 0
    return
end
sen_set = handles.sen_set;
g = handles.g;

if g.hh > 12
    % Current date
    date = str2double(sprintf('%s%s%s',g.YYYY{:},g.MM{:},g.DD{:}));
    msg = sprintf('There are some time shift conflicts.\n\nSen          Curr            New\n');
else
    % Previous date
    date = sprintf('%s%s%s',g.YYYY{:},g.MM{:},g.DD{:});
    daten = datenum(date,'yyyymmdd');
    date = str2double(datestr(daten-1,'yyyymmdd'));
    msg = sprintf('There are some time shift conflicts. \n\nSince the plotting time is before 12.00h, time shifts information were loaded from the previous date %s.',datestr(daten-1,'yyyy-mmm-dd'));
    msg = sprintf('%s\n\nSen          Curr            New\n',msg);
end


index = find(handles.tshift(:,1)== date);

% replace NaNs with zeros
handles.tshift(isnan(handles.tshift)) = 0;
conflicts = 0;

if ~isempty(index);
    % let's compare with current values
   
    for i=1:20
        if handles.tshift(index,i+1) ~= sen_set.t_shift(i*3)
            msg = sprintf('%s%s          %5.0f          %5.0f\n',...
                msg,sen_set.sen_IDs{i},...
                sen_set.t_shift(i*3),...
                handles.tshift(index,i+1));
            conflicts = 1;
        end
    end
end

if conflicts
    choice = questdlg(msg, ...
	'Time Shift Conflicts..', ...
	'Change','Stop warn me!','Cancel','Change');
else
    return
end

switch choice

    case 'Cancel'
        return
    case 'Stop warn me!'
        handles.tshift_war = 0;
    case 'Change'
        for i=1:20
            sen_set.t_shift(i*3-1) = handles.tshift(index,i+1);
            sen_set.t_shift(i*3-2) = handles.tshift(index,i+1);
            sen_set.t_shift(i*3)   = handles.tshift(index,i+1);
        end
        handles.sen_set = sen_set;
        save('sensor_setting.mat','-struct','sen_set');      
        
end
        

% --------------------------------------------------------------------
function max_altitude_Callback(hObject, eventdata, handles)

sen_set = handles.sen_set;

prompt={'Enter MAX Altitude:',...
        'Enter MIN Altitude:'};
name='Altitudes';
numlines=1;
defaultanswer={num2str(sen_set.max_z),num2str(sen_set.min_z)};
answer=inputdlg(prompt,name,numlines,defaultanswer);

answer = str2double(answer);

if ~isempty(answer);
    sen_set.max_z = answer(1);
    sen_set.min_z = answer(2);
    save('sensor_setting.mat','-struct','sen_set')
    handles.sen_set = sen_set;
    guidata(hObject, handles);  
end
    
    


% --------------------------------------------------------------------
function undo_Callback(hObject, eventdata, handles)

temp = handles.temp;

if temp.UR_cindex == 2
    temp.UR_cindex = temp.UR_cindex - 1;
    handles.g = temp.UR_data(temp.UR_cindex);
    set_values(handles);    
    set(handles.redo,'Enable','on')
    set(handles.undo,'Enable','off')
elseif temp.UR_cindex > 2
    temp.UR_cindex = temp.UR_cindex - 1;
    handles.g = temp.UR_data(temp.UR_cindex);
    set_values(handles);    
    set(handles.redo,'Enable','on')
end

handles.temp = temp;
guidata(hObject, handles);  


% --------------------------------------------------------------------
function redo_Callback(hObject, eventdata, handles)

temp = handles.temp;

if (temp.UR_index - temp.UR_cindex) == 1
    temp.UR_cindex = temp.UR_cindex + 1;
    handles.g = temp.UR_data(temp.UR_cindex);
    set_values(handles);    
    set(handles.redo,'Enable','off')
    set(handles.undo,'Enable','on')
elseif temp.UR_cindex < temp.UR_index
    temp.UR_cindex = temp.UR_cindex + 1;
    handles.g = temp.UR_data(temp.UR_cindex);
    set_values(handles);    
    set(handles.undo,'Enable','on')
end

handles.temp = temp;
guidata(hObject, handles);  


% --------------------------------------------------------------------
function ch_high_pass_filter_Callback(hObject, eventdata, handles)

discr = 'Enter High-Pass filter freqency in Hz. Value "0" will disable the filter';


[values, button] = settingsdlg(...
      'Description', discr,... 
      'title'      , 'Ch High-Pass Filter',...
      {'Frequency';'ch_hpf'}, handles.g.temp.ch_hpf ...
      );
  
if strcmp(button, 'ok')
    handles.g.temp.ch_hpf = values.ch_hpf;
    guidata(hObject, handles)
end


 


% --------------------------------------------------------------------
function long_range_ldar_Callback(hObject, eventdata, handles)

discr = 'Entert t1, t2 values from 0 -- 86400s from midnight and XY plot of LDAR2 will be created with color coded by time';

[ts, button] = settingsdlg(...
      'Description', discr,... 
      'title'      , 'Long Range LDAR2',...
      {'t1';'t1'}, handles.g.t1, ...
      {'t2';'t2'}, handles.g.t2, ...
      {'Histograms','hgrms'},true ...
      );

if ~strcmp(button, 'ok'); return; end

ts.t1 = 60000;
ts.t2 = 86400;

wbh = waitbar(0.1,'Loading LDAR2 data...','name','Long Range LDAR2 plot');

g = handles.g;
settings = handles.sen_set;

t = floor(ts.t1/1800)*1800;

% Time shift for LDAR and Linet
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
DLSt = [];  DLSx = [];  DLSy = [];
CGt  = [];  CGx  = [];  CGy  = [];


while t < ts.t2
   
    % Generating file name for LDAR data
    [hms hh mm] = sec2hhmmss(t);
    if mm < 30
        ext=0;
    else
        ext=30;
    end
    
    dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
        -datenum(str2double(g.YYYY),0,0));
    
    ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,hh,ext);
    
    % Exit if the file does not exsist
    if ~exist(ldar_fn,'file')
        errordlg(ldar_fn,'File not found')
    end
    
    %fprintf('%s\n',ldar_fn)
    
    % Load LDAR2 data
    
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,ts.t1,ts.t2,str2double(settings.ldar_r),...
        x0,y0,z0,0);
   
   DLSt = [DLSt; DLS(:,10)];
   DLSx = [DLSx; DLS(:,6)];
   DLSy = [DLSy; DLS(:,7)];
   
   CGt  = [CGt; CG(:,10)];
   CGx  = [CGx ; CG(:,6)];
   CGy  = [CGy ; CG(:,7)];
    
    t = t+1800;
end

DLSx = DLSx/1000;
DLSy = DLSy/1000;
CGx  = CGx/1000;
CGy  = CGy/1000;


% nothing to plot? exit
if sum(isnan(DLSt)) == length(DLSt) && ...
        sum(isnan(CGt)) == length(CGt)
    disp('No data in given time range')
    delete(wbh)
    return
end

t = [DLSt; CGt];

fg = figure;
set(fg,'visible','off');
hold all


map=colormap;
miv=min(t);
mav=max(t);
clrstep = (mav-miv)/size(map,1);

lg_ldar = 0;
lg_cglss = 0;
lg = {};

% marker size
mz = 4;

% Total number of data points
Lt = length(t);
% Number of points plotted so far
Lp = 0;

if g.ldar
for nc=1:size(map,1)
        iv = find(DLSt>=miv+(nc-1)*clrstep & DLSt<=miv+nc*clrstep);
        Liv = length(iv);
        Lp = Lp + Liv;
        wbval = 0.1+0.9*Lp/Lt;
        wbstr = sprintf('Plotting LDAR2 points ... (%0.2f%%)',wbval*100);
        waitbar(wbval,wbh,wbstr)
        
        hLine = plot(DLSx(iv),DLSy(iv),'*','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
        
        if lg_ldar
            %Exclude from legend
            set(get(get(hLine,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');
        else
            lg = [lg sprintf('LDAR2-%i',length(DLSx))];
            lg_ldar = 1;
        end    
end
end


if g.cglss
for nc=1:size(map,1)
        iv = find(CGt>=miv+(nc-1)*clrstep & CGt<=miv+nc*clrstep);
        Liv = length(iv);
        Lp = Lp + Liv;
        wbval = 0.1+0.9*Lp/Lt;
        wbstr = sprintf('Plotting LDAR2 points ... (%0.2f%%)',wbval*100);
        waitbar(wbval,wbh,wbstr)
        
        hLine = plot(CGx(iv),CGy(iv),'+','color',map(nc,:),'MarkerSize',mz);
        
        if lg_cglss
            %Exclude from legend
            set(get(get(hLine,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');
        else
            lg = [lg sprintf('CGLSS-%i',length(CGx))];
            lg_cglss = 1;
        end    
end
end


% Re-format the colorbar
h=colorbar;


%set(h,'ylim',[1 length(map)]);
yal=linspace(1,length(map),10);
set(h,'ytick',yal);
% Create the yticklabels
ytl=linspace(miv,mav,10);
s=char(10,4);

for i=1:10
    if min(abs(ytl)) >= 0.001
        B=sprintf('%-4.0f',ytl(i));
    else
        B=sprintf('%-3.0E',ytl(i));
    end
    s(i,1:length(B))=B;
end

set(h,'yticklabel',s);

box on
daspect([ 1 1 1])
legend(lg)

if settings.mapOn
    florida_map
end

title(sprintf('%s-%s-%s         %0.1fs -- %0.1fs',g.YYYY{:},g.MM{:},g.DD{:},ts.t1,ts.t2))

delete(wbh)
set(fg,'visible','on')

if ts.hgrms
    figure
    subplot(2,1,1)
    hist(DLSt,miv:clrstep:mav)
    title(sprintf('LDAR2 Histogram %s-%s-%s         %0.1fs -- %0.1fs',g.YYYY{:},g.MM{:},g.DD{:},ts.t1,ts.t2))
    xlabel('Time (s)')
    ylabel('Number of soueces')
    
    
    subplot(2,1,2)
    hist(CGt,miv:clrstep:mav)
    title(sprintf('CGLSS Histogram %s-%s-%s         %0.1fs -- %0.1fs',g.YYYY{:},g.MM{:},g.DD{:},ts.t1,ts.t2))
    xlabel('Time (s)')
    ylabel('Number of soueces')
end


% --------------------------------------------------------------------
function calibrate_cross_lp_Callback(hObject, eventdata, handles)
    ip_cross_calibrator2(handles.g)


% --------------------------------------------------------------------
function pbfa_error_bar_on_Callback(hObject, eventdata, handles)
if strcmp(get(gcbo, 'Checked'),'on')
    set(gcbo, 'Checked', 'off');
    handles.sen_set.pbfa_ebar_on = 0;    
else
    set(gcbo, 'Checked', 'on');
    handles.sen_set.pbfa_ebar_on = 1;
end
sen_set = handles.sen_set;
guidata(hObject, handles)
save('sensor_setting.mat','-struct','sen_set')


% --------------------------------------------------------------------
function tgf_plotter_Callback(hObject, eventdata, handles)
    tgf_plotter


% --- Executes on button press in pbfa_old.
function pbfa_old_Callback(hObject, eventdata, handles)

handles.g.pbfa_old = (get(hObject,'Value'));
handles.g.isSaved = 0;
guidata(hObject, handles);


% --- Executes on button press in nldn2.
function nldn2_Callback(hObject, eventdata, handles)
handles.g.nldn2=(get(hObject,'Value'));
handles.g.isSaved=0;
guidata(hObject, handles);


% --------------------------------------------------------------------
function auto_b_folder_Callback(hObject, eventdata, handles)
sen_set = handles.sen_set;
val = get(handles.auto_b_folder,'Checked');

if strcmp(val,'off')
    set(handles.auto_b_folder,'Checked','on');
    sen_set.autoBaseFolder = 'on';    
else
    set(handles.auto_b_folder,'Checked','off');
    sen_set.autoBaseFolder = 'off';
end

save('sensor_setting.mat','-Struct','sen_set')

handles.sen_set = sen_set;
guidata(hObject, handles);


function sen_set = set_location_folder(sen_set,cname)

switch cname
    case 'hst000068015um.phy.olemiss.edu';    % Dr. Marshall's computer  
        sen_set.loc_dir = '/Users/thomas/Documents/KSC_2011_data/LocData';
    case 'hst000068047um.phy.olemiss.edu';     % Dr. Marshall's new computer
        sen_set.loc_dir = '/Users/thomas/Documents/KSC_2011_data/LocData';        
    otherwise
        sen_set.loc_dir = sen_set.base_dir;
end


function handles = set_auto_b_folder(handles)

sen_set = handles.sen_set;
date = str2double([handles.g.YYYY{:} handles.g.MM{:} handles.g.DD{:}]);
prefix = handles.prefix;
  
if date <= 20110627
    sen_set.base_dir = [prefix '/data/2010-07-01 -- 2010-08-19/data/'];
elseif date <= 20110716
    sen_set.base_dir = [prefix '/data/2011-07-07 -- 2011-07-16/data/'];
elseif date <= 20110804
    sen_set.base_dir = [prefix '/data/2011-07-17 -- 2011-08-04/data/'];
elseif date <= 20110816
    sen_set.base_dir = [prefix '/data/2011-08-05 -- 2011-08-16/data/'];
end

sen_set = set_location_folder(sen_set,handles.cname);
save('sensor_setting.mat','-struct', 'sen_set');
handles.sen_set=sen_set;

function handles = set_prefix(handles)
    switch handles.cname
        case 'sadaq7';        handles.prefix = 'C:';  
        case 'sadaq9';        handles.prefix = '/media'; % Ubuntu server
        case 'sumedhe-HP';    handles.prefix = '\\SADAQ7'; % Sumedhe's laptop       
        case 'hst000068015um.phy.olemiss.edu';    % Dr. Marshall's computer  
            handles.prefix = '/Volumes';
        case 'hst000068047um.phy.olemiss.edu';    % Dr. Marshall's new computer
            handles.prefix = '/Volumes';
        otherwise;
            handles.prefix = '';
    end

% --------------------------------------------------------------------
function PBFA_auto_Callback(hObject, eventdata, handles)

PBFA_auto4(handles)


% --------------------------------------------------------------------
function remove_electronic_v_decay_Callback(hObject, eventdata, handles)

sen_set = handles.sen_set;
val = get(hObject,'Checked');

if strcmp(val,'off')
    set(handles.remove_electronic_v_decay,'Checked','on');
    sen_set.remove_time_decay = 'on';    
else
    set(handles.remove_electronic_v_decay,'Checked','off');
    sen_set.remove_time_decay = 'off';
end

save('sensor_setting.mat','-Struct','sen_set')

handles.sen_set = sen_set;
guidata(hObject, handles);


% --------------------------------------------------------------------
function nearest_pulse_settings_Callback(hObject, eventdata, handles)
a=handles.sen_set;

%Discription
discrip =sprintf('Settings for nearest pulses/flashes finding\nbased on LDAR2 data');

[choise, button] = settingsdlg(...
    'Description', discrip ,...
    'title' , 'Finding nearest pulses',...
    {'Use xyz (if not will use only xy)'; 'use_xyz'},logical(a.nearest.use_xyz), ...
    {'Range in meters', 'R'}, num2str(a.nearest.R), ...
    {'Plot time range (+/-s)', 'dt'}, num2str(a.nearest.dt), ...    
    {'Ingnore previous pulses within (us)', 'idt1'}, num2str(a.nearest.idt1), ...
    {'Ingnore pulses after within (us)', 'idt2'}, num2str(a.nearest.idt2), ...
     'WindowWidth' , 300,...
    'ControlWidth', 100);


if strcmp(button,'ok')
    
    a.nearest.use_xyz = choise.use_xyz;
    a.nearest.R = choise.R;
    a.nearest.dt = choise.dt;
    a.nearest.idt1 = choise.idt1;
    a.nearest.idt2 = choise.idt2;
    
  
    save('sensor_setting.mat','-struct','a')
    handles.sen_set=a;
    handles=calibration(handles);    
    guidata(hObject, handles);
   
end


% --------------------------------------------------------------------
function FM_plotter_Callback(hObject, eventdata, handles)

    field_mill_plotter(handles)

function [radarFiles , radarFullFiles] = getRadarFiles(handles)

g = handles.g;
sen_set = handles.sen_set;

% Radar file folder
rdn = sprintf('%s/netCDF/%s/%s/%s/%s/%2.2i/',...
    sen_set.base_dir,sen_set.radarStations{sen_set.radarStationID},...
    g.YYYY{:},g.MM{:},g.DD{:},g.hh);

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
    end
end

% if within the first 10 mins of the hour, let's load files from the
% previous hour
if g.mm < 10
    dNum1 = datenum(str2double(g.YYYY{:}),str2double(g.MM{:}),str2double(g.DD{:}),g.hh,0,0);
    
    % substract an hour
    dNum2 = dNum1 - 1/24;
    
    dvec = datevec(dNum2);
    
    % Radar file folder
    rdn = sprintf('%s/netCDF/%4.4i/%2.2i/%2.2i/%2.2i/',...
        sen_set.base_dir,dvec(1),dvec(2),dvec(3),dvec(4));
    
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
if g.mm > 50
    dNum1 = datenum(str2double(g.YYYY{:}),str2double(g.MM{:}),str2double(g.DD{:}),g.hh,0,0);
    
    % add an hour
    dNum2 = dNum1 + 1/24;
    
    dvec = datevec(dNum2);
    
    % Radar file folder
    rdn = sprintf('%s/netCDF/%4.4i/%2.2i/%2.2i/%2.2i/',...
        sen_set.base_dir,dvec(1),dvec(2),dvec(3),dvec(4));
    
    files = dir(rdn);
    
    L = length(files);
    
    % We just need the first file of this folder
    if L > 2
        fn =  files(3).name;
        times =[times, str2double(fn(14:15))*3600+str2double(fn(16:17))*60+str2double(fn(18:19))];
        radarFiles = [radarFiles; fn];
        radarFullFiles = [radarFullFiles; [rdn fn]];
    end
end

if ~isempty(times)
    dts = abs(times - g.t1);
    [dts I] = sort(dts);
    radarFiles = radarFiles(I);
    radarFullFiles = radarFullFiles(I);
else
    radarFiles = '';
    radarFullFiles = '';
end


% --------------------------------------------------------------------
function location_base_folder_Callback(hObject, eventdata, handles)

sen_set=handles.sen_set;

if sen_set.loc_dir ==0
    sen_set.loc_dir='None';
end

msg=sprintf('Choose the bease directory for slow and fast antenna data\n\nCurrent base dir = %s',...
    sen_set.loc_dir);

sen_set.loc_dir=uigetdir(sen_set.loc_dir,msg);
save('sensor_setting.mat','-struct', 'sen_set');

handles.sen_set=sen_set;
guidata(hObject, handles);


% --------------------------------------------------------------------
function radar_station_Callback(hObject, eventdata, handles)

sen_set = handles.sen_set;

%Discription
discrip =sprintf('Set RADAR station. Current station is sat to %s.',...
    sen_set.radarStations{sen_set.radarStationID});

[choise, button] = settingsdlg(...
    'Description', discrip ,...
    'title' , 'RADAR Station setup',...
    {'Station ID', 'station_ID'}, sen_set.radarStations, ...
     'WindowWidth' , 200,...
    'ControlWidth', 100);


if strcmp(button,'ok')
    sen_set.radarStationID = find(strcmp(sen_set.radarStations,choise.station_ID));
    
    % Show Current radar station on settings
    set(handles.radar_station,'Label',['Radara station [' ...
        sen_set.radarStations{sen_set.radarStationID} ']']);
    
    save('sensor_setting.mat','-struct', 'sen_set');
    
    handles.sen_set=sen_set;
    guidata(hObject, handles);
end
    


% --------------------------------------------------------------------
function ldar_box_Callback(hObject, eventdata, handles)

sen_set = handles.sen_set;
ldar_box = sen_set.ldar_box;

%Discription
discrip =sprintf('Set LDAR Box coners. Setting this will only loads LDARs/CGLSS/PBFA,etc in box.');

[choise, button] = settingsdlg(...
    'Description', discrip ,...
    'title' , 'LDAR box setup',...
    {'Turn ON'; 'on'}, boolean(ldar_box(1)),...
    {'Left-Bottom x (km)', 'lbx'}, ldar_box(2), ...
    {'Left-Bottom y (km)', 'lby'}, ldar_box(3), ...
    {'Horizontal lengh (km)', 'hl'}, ldar_box(4), ...
    {'Vertical height (km)', 'vh'}, ldar_box(5), ...
     'WindowWidth' , 250,...
     'ControlWidth', 75);

if strcmp(button,'ok')
    sen_set.ldar_box(1) = choise.on;
    sen_set.ldar_box(2) = choise.lbx;
    sen_set.ldar_box(3) = choise.lby;
    sen_set.ldar_box(4) = choise.hl;
    sen_set.ldar_box(5) = choise.vh;
    
    save('sensor_setting.mat','-struct', 'sen_set')
    
    handles.sen_set = sen_set;
    guidata(hObject, handles);    
end
    


% --- Executes on button press in checkbox282.
function checkbox282_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox282 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox282


% --- Executes on button press in checkbox283.
function checkbox283_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox283 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox283


% --- Executes on button press in checkbox284.
function checkbox284_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox284 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox284


% --- Executes on button press in checkbox285.
function checkbox285_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox285 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox285


% --- Executes on selection change in popupmenu17.
function popupmenu17_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu17 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu17


% --- Executes during object creation, after setting all properties.
function popupmenu17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function xy_marker_size_Callback(hObject, eventdata, handles)
% hObject    handle to xy_marker_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sen_set=handles.sen_set;

prompt =sprintf('x-y plot marker size (an integer, default 4) for \n    LDAR    LINET\n    PBFA    NLDN');
dlg_title = 'Marker Size';
num_lines = 1;
def = sen_set.xy_marker_size;
answer = inputdlg(prompt,dlg_title,num_lines,def);

if ~isempty(answer)
    sen_set.xy_marker_size = answer;
    save('sensor_setting.mat','-struct', 'sen_set')
    handles.sen_set=sen_set;
    guidata(hObject, handles);
end


% --------------------------------------------------------------------
function xt_marker_size_Callback(hObject, eventdata, handles)
% hObject    handle to xt_marker_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


sen_set=handles.sen_set;

prompt =sprintf('x-t plot marker size (an integer, default 4) for \n    LDAR    LINET\n    PBFA    NLDN');
dlg_title = 'Marker Size';
num_lines = 1;
def = sen_set.xt_marker_size;
answer = inputdlg(prompt,dlg_title,num_lines,def);

if ~isempty(answer)
    sen_set.xt_marker_size = answer;
    save('sensor_setting.mat','-struct', 'sen_set')
    handles.sen_set=sen_set;
    guidata(hObject, handles);
end
