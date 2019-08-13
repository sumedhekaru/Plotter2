function varargout = pulse_modeler(varargin)
% PULSE_MODELER MATLAB code for pulse_modeler.fig
%      PULSE_MODELER, by itself, creates a new PULSE_MODELER or raises the existing
%      singleton*.
%
%      H = PULSE_MODELER returns the handle to a new PULSE_MODELER or the handle to
%      the existing singleton*.
%
%      PULSE_MODELER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PULSE_MODELER.M with the given input arguments.
%
%      PULSE_MODELER('Property','Value',...) creates a new PULSE_MODELER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pulse_modeler_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pulse_modeler_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATAc, GUIHANDLES

% Edit the above text to modify the response to help pulse_modeler

% Last Modified by GUIDE v2.5 20-Aug-2015 08:12:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pulse_modeler_OpeningFcn, ...
                   'gui_OutputFcn',  @pulse_modeler_OutputFcn, ...
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


% --- Executes just before pulse_modeler is made visible.
function pulse_modeler_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pulse_modeler (see VARARGIN)

% Choose default command line output for pulse_modeler
handles.output = hObject;



try
    % Try to open last gui data
    d = open('pulse_modeler_last_gui_data.mat');
    handles.d = d;
    update_GUI_vals(handles);
catch
    % Reset to defaults
    handles = push_Reset_Callback(hObject, eventdata, handles);
end

% open matlabpool if it is not enable

% % This code seems to not work with R2015a
% if matlabpool('size') == 0
%     % close any exisiting fist
%     matlabpool close force
%     try; matlabpool open; end;
% end
   
% %Used the follwoing code instead of above few line for R2015a
% if isempty(gcp('nocreate'))
%     % close any exisiting fist
%     %parpool close force
%     try; parpool open; end;
% end


handles.best_set.d = [];
handles.best_set.diff = NaN;
handles.best_set.norm_diff = NaN;

handles.diffs.CurrDiff1= NaN;
handles.diffs.CurrNormDiff1= NaN;
handles.diffs.CurrDiff2= NaN;
handles.diffs.CurrNormDiff2= NaN;
handles.diffs.PrevDiff = NaN;
handles.diffs.PrevNormDiff = NaN;
handles.diffs.BestDiff = NaN;
handles.diffs.BestNormDiff = NaN;
handles.diffs.d1 = [];
handles.diffs.d2 = [];

% Settings
handles.sen_set = open('sensor_setting.mat');

% Turn off variable initiation parfor warning
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pulse_modeler wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function update_GUI_vals(handles)

h = handles;
d = h.d;
set(h.TxtMsg,'String',d.msg)
set(h.dir,'String',d.dir)
set(h.sns_str,'String',d.sns_str)
set(h.edit_dt,'String',d.dt)
set(h.edit_dh,'String',d.dh)
set(h.edit_T1,'String',d.T1)
set(h.edit_T2,'String',d.T2)
set(h.edit_x0,'String',d.x0)
set(h.edit_y0,'String',d.y0)
set(h.edit_H1,'String',d.H1)
set(h.edit_Hm,'String',d.Hm)
set(h.edit_H2,'String',d.H2)
set(h.edit_t1,'String',d.t1)
set(h.edit_t2,'String',d.t2)
set(h.edit_v,'String',d.v)
set(h.edit_lamda1,'String',d.lamda1)
set(h.edit_lamda2,'String',d.lamda2)
set(h.edit_AlpInd,'String',d.AlpInd)

set(h.checkDynamic,'Value',d.dynamicMod)
set(h.popModalType,'Value',d.modal_num)
set(h.popCurrType,'Value',d.currType_num)

set(h.edit_dx0,'String',d.dx0)
set(h.edit_dy0,'String',d.dy0)
set(h.edit_dH1,'String',d.dH1)
set(h.edit_dHm,'String',d.dHm)
set(h.edit_dH2,'String',d.dH2)
set(h.edit_dt1,'String',d.dt1)
set(h.edit_dt2,'String',d.dt2)
set(h.edit_dv,'String',d.dv)
set(h.edit_dLamda1,'String',d.dLamda1)
set(h.edit_dLamda2,'String',d.dLamda2)
set(h.edit_dAlpInd,'String',d.dAlpInd)


try
    set(h.edit_tt2,'String',d.tt2)
    set(h.edit_dtt2,'String',d.dtt2)
    set(h.edit_t3,'String',d.t3)
    set(h.edit_dt3,'String',d.dt3)
end
    

enable_items(handles)




% --- Outputs from this function are returned to the command line.
function varargout = pulse_modeler_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
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
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function create_data_Callback(hObject, eventdata, handles)
% hObject    handle to create_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function full_screen_Callback(hObject, eventdata, handles)
% hObject    handle to full_screen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function dynamic_modeling_Callback(hObject, eventdata, handles)
% hObject    handle to dynamic_modeling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function data_folder_Callback(hObject, eventdata, handles)
% hObject    handle to data_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function open_Callback(hObject, eventdata, handles)

d = handles.d;

[filename,pathname] = uigetfile(sprintf('%s%s-.mat',d.dir,d.modal));

if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   % remember current values
   d0 = handles.d;
   handles.d = open(fullfile(pathname,filename));
   try
        update_GUI_vals(handles);
   catch
       % restore defaults
       handles.d = d0;
       update_GUI_vals(handles);
       msgbox('Contents of the file you choosed can not be verified.','Open failier')
   end
end
guidata(hObject, handles)

% --------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)

disp('Use save as instead...')


% --------------------------------------------------------------------
function save_as_Callback(hObject, eventdata, handles)

d = handles.d;

fn = sprintf('%s%s-nDiff-%0.2e.mat',d.dir,d.modal,handles.diffs.CurrNormDiff1);

[filename,pathname] = uiputfile(fn);

if isequal(filename,0) || isequal(pathname,0)
   % disp('User selected Cancel')
else
   % make sure to have save some modeling parameters too.
   
   save(fullfile(pathname,filename),'-Struct','d')
end




% --------------------------------------------------------------------
function save_figure_as_Callback(hObject, eventdata, handles)
% hObject    handle to save_figure_as (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_H1_Callback(hObject, eventdata, handles)
handles.d.H1 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_H1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_H1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_dH1_Callback(hObject, eventdata, handles)
handles.d.dH1 = str2double(get(hObject,'String'));

if handles.d.dynamicMod
    d1 = handles.d;
    d2 = handles.d;
    d1.H1 = d1.H1 - d1.dH1;
    d2.H1 = d2.H1 + d2.dH1;
    handles.d1 = d1;
    handles.d2 = d2;
    handles = modal_fineTune(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dH1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dH1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_H1.
function fineTune_H1_Callback(hObject, eventdata, handles)
d1 = handles.d;
d2 = handles.d;
d1.H1 = d1.H1 - d1.dH1;
d2.H1 = d2.H1 + d2.dH1;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);


function edit_dHm_Callback(hObject, eventdata, handles)
handles.d.dHm = str2double(get(hObject,'String'));

if handles.d.dynamicMod
    d1 = handles.d;
    d2 = handles.d;
    d1.Hm = d1.Hm - d1.dHm;
    d2.Hm = d2.Hm + d2.dHm;
    handles.d1 = d1;
    handles.d2 = d2;
    handles = modal_fineTune(handles);
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_dHm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dHm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_Hm.
function fineTune_Hm_Callback(hObject, eventdata, handles)
d1 = handles.d;
d2 = handles.d;
d1.Hm = d1.Hm - d1.dHm;
d2.Hm = d2.Hm + d2.dHm;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);


function edit_Hm_Callback(hObject, eventdata, handles)
handles.d.Hm = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_Hm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Hm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dH2_Callback(hObject, eventdata, handles)
handles.d.dH2 = str2double(get(hObject,'String'));

if handles.d.dynamicMod
    d1 = handles.d;
    d2 = handles.d;
    d1.H2 = d1.H2 - d1.dH2;
    d2.H2 = d2.H2 + d2.dH2;
    handles.d1 = d1;
    handles.d2 = d2;
    handles = modal_fineTune(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dH2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dH2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_H2.
function fineTune_H2_Callback(hObject, eventdata, handles)
d1 = handles.d;
d2 = handles.d;
d1.H2 = d1.H2 - d1.dH2;
d2.H2 = d2.H2 + d2.dH2;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);



function edit_H2_Callback(hObject, eventdata, handles)
handles.d.H2 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_H2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_H2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dt1_Callback(hObject, eventdata, handles)

handles.d.dt1 = str2double(get(hObject,'String'));

if handles.d.dynamicMod
    d1 = handles.d;
    d2 = handles.d;
    d1.t1 = d1.t1 - d1.dt1;
    d2.t1 = d2.t1 + d2.dt1;
    handles.d1 = d1;
    handles.d2 = d2;
    handles = modal_fineTune(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dt1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_t1.
function fineTune_t1_Callback(hObject, eventdata, handles)
d1 = handles.d;
d2 = handles.d;
d1.t1 = d1.t1 - d1.dt1;
d2.t1 = d2.t1 + d2.dt1;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);


function edit_t1_Callback(hObject, eventdata, handles)
handles.d.t1 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_t1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_t1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_dt_Callback(hObject, eventdata, handles)
handles.d.dt = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_dh_Callback(hObject, eventdata, handles)
handles.d.dh = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_dh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dx0_Callback(hObject, eventdata, handles)
handles.d.dx0 = str2double(get(hObject,'String'));

if handles.d.dynamicMod
    d1 = handles.d;
    d2 = handles.d;
    d1.x0 = d1.x0 - d1.dx0;
    d2.x0 = d2.x0 + d2.dx0;
    handles.d1 = d1;
    handles.d2 = d2;
    handles = modal_fineTune(handles);
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_dx0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dx0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_x0.
function fineTune_x0_Callback(hObject, eventdata, handles)

d1 = handles.d;
d2 = handles.d;
d1.x0 = d1.x0 - d1.dx0;
d2.x0 = d2.x0 + d2.dx0;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);



function edit_x0_Callback(hObject, eventdata, handles)
handles.d.x0 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dy0_Callback(hObject, eventdata, handles)
handles.d.dy0 = str2double(get(hObject,'String'));

if handles.d.dynamicMod
    d1 = handles.d;
    d2 = handles.d;
    d1.y0 = d1.y0 - d1.dy0;
    d2.y0 = d2.y0 + d2.dy0;
    handles.d1 = d1;
    handles.d2 = d2;
    handles = modal_fineTune(handles);
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_dy0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dy0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_y0.
function fineTune_y0_Callback(hObject, eventdata, handles)
d1 = handles.d;
d2 = handles.d;
d1.y0 = d1.y0 - d1.dy0;
d2.y0 = d2.y0 + d2.dy0;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);



function edit_y0_Callback(hObject, eventdata, handles)
handles.d.y0 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_y0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dt2_Callback(hObject, eventdata, handles)
handles.d.dt2 = str2double(get(hObject,'String'));

if handles.d.dynamicMod
    d1 = handles.d;
    d2 = handles.d;
    d1.t2 = d1.t2 - d1.dt2;
    d2.t2 = d2.t2 + d2.dt2;
    handles.d1 = d1;
    handles.d2 = d2;
    handles = modal_fineTune(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_t2.
function fineTune_t2_Callback(hObject, eventdata, handles)
d1 = handles.d;
d2 = handles.d;
d1.t2 = d1.t2 - d1.dt2;
d2.t2 = d2.t2 + d2.dt2;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);



function edit_t2_Callback(hObject, eventdata, handles)
handles.d.t2 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_t2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_t2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dLamda1_Callback(hObject, eventdata, handles)
handles.d.dLamda1 = str2double(get(hObject,'String'));

if handles.d.dynamicMod
    d1 = handles.d;
    d2 = handles.d;
    d1.lamda1 = d1.lamda1 - d1.dLamda1;
    d2.lamda1 = d2.lamda1 + d2.dLamda1;
    handles.d1 = d1;
    handles.d2 = d2;
    handles = modal_fineTune(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dLamda1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dLamda1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_lamda1.
function fineTune_lamda1_Callback(hObject, eventdata, handles)
d1 = handles.d;
d2 = handles.d;
d1.lamda1 = d1.lamda1 - d1.dLamda1;
d2.lamda1 = d2.lamda1 + d2.dLamda1;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);



function edit_lamda1_Callback(hObject, eventdata, handles)
handles.d.lamda1 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_lamda1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lamda1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dLamda2_Callback(hObject, eventdata, handles)
handles.d.dLamda2 = str2double(get(hObject,'String'));

if handles.d.dynamicMod
    d1 = handles.d;
    d2 = handles.d;
    d1.lamda2 = d1.lamda2 - d1.dLamda2;
    d2.lamda2 = d2.lamda2 + d2.dLamda2;
    handles.d1 = d1;
    handles.d2 = d2;
    handles = modal_fineTune(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dLamda2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dLamda2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_lamda2.
function fineTune_lamda2_Callback(hObject, eventdata, handles)
d1 = handles.d;
d2 = handles.d;
d1.lamda2 = d1.lamda2 - d1.dLamda2;
d2.lamda2 = d2.lamda2 + d2.dLamda2;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);



function edit_lamda2_Callback(hObject, eventdata, handles)
handles.d.lamda2 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_lamda2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lamda2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dAlpInd_Callback(hObject, eventdata, handles)
handles.d.dAlpInd = str2double(get(hObject,'String'));

if handles.d.dynamicMod
    d1 = handles.d;
    d2 = handles.d;
    d1.AlpInd = d1.AlpInd - d1.dAlpInd;
    d2.AlpInd = d2.AlpInd + d2.dAlpInd;
    handles.d1 = d1;
    handles.d2 = d2;
    handles = modal_fineTune(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dAlpInd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dAlpInd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_AppInd.
function fineTune_AppInd_Callback(hObject, eventdata, handles)
d1 = handles.d;
d2 = handles.d;
d1.AlpInd = d1.AlpInd - d1.dAlpInd;
d2.AlpInd = d2.AlpInd + d2.dAlpInd;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);


function edit_AlpInd_Callback(hObject, eventdata, handles)
handles.d.AlpInd = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_AlpInd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AlpInd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in push_modal.
function push_modal_Callback(hObject, eventdata, handles)

handles = model(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkDynamic.
function checkDynamic_Callback(hObject, eventdata, handles)

handles.d.dynamicMod = get(hObject,'Value');
guidata(hObject, handles);


function dir_Callback(hObject, eventdata, handles)
% hObject    handle to dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dir as text
%        str2double(get(hObject,'String')) returns contents of dir as a double


% --- Executes during object creation, after setting all properties.
function dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in getDir.
function getDir_Callback(hObject, eventdata, handles)

dirN = uigetdir(handles.d.dir);

if ischar(dirN)
    handles.d.dir = [dirN '/'];
    set(handles.dir,'String',handles.d.dir)
    guidata(hObject, handles);
end





% --- Executes on selection change in popModalType.
function popModalType_Callback(hObject, eventdata, handles)

contents = cellstr(get(hObject,'String'));
handles.d.modal_num = get(hObject,'Value');
handles.d.modal = contents{get(hObject,'Value')};

% set some values disabled according to modal type
enable_items(handles)

if handles.d.dynamicMod
    handles = model(handles);
end

guidata(hObject, handles);

function enable_items(handles)

switch handles.d.modal
    case 'MTLL'
        % we don't need Hm, lamda1, Lmda2
        
        set(handles.text_H1,'enable', 'on')
        set(handles.edit_H1,'enable', 'on')
        set(handles.pm_H1,'enable', 'on')
        set(handles.edit_dH1,'enable', 'on')
        set(handles.fineTune_H1,'enable', 'on')
        
        set(handles.text_Hm,'enable', 'off')
        set(handles.edit_Hm,'enable', 'off')
        set(handles.pm_Hm,'enable', 'off')
        set(handles.edit_dHm,'enable', 'off')
        set(handles.fineTune_Hm,'enable', 'off')
        
        set(handles.text_lamda1,'String', 'Lamda 1 (m)')
        set(handles.text_lamda1,'enable', 'off')
        set(handles.edit_lamda1,'enable', 'off')
        set(handles.pm_lamda1,'enable', 'off')
        set(handles.edit_dLamda1,'enable', 'off')
        set(handles.fineTune_lamda1,'enable', 'off')
        
        set(handles.text_lamda2,'String', 'Lamda 2 (m)')
        set(handles.text_lamda2,'enable', 'off')
        set(handles.edit_lamda2,'enable', 'off')
        set(handles.pm_lamda2,'enable', 'off')
        set(handles.edit_dLamda2,'enable', 'off')
        set(handles.fineTune_lamda2,'enable', 'off')
        
        set(handles.edit_tt2,'enable', 'off')
        set(handles.edit_dtt2,'enable', 'off')
        set(handles.pm_tt2,'enable', 'off')
        set(handles.fineTune_tt2,'enable', 'off')
        
        
        set(handles.edit_t3,'enable', 'off')
        set(handles.edit_dt3,'enable', 'off')
        set(handles.pm_t3,'enable', 'off')
        set(handles.fineTune_t3,'enable', 'off')
        
    case 'MTLE'
        set(handles.text_H1,'enable', 'on')
        set(handles.edit_H1,'enable', 'on')
        set(handles.pm_H1,'enable', 'on')
        set(handles.edit_dH1,'enable', 'on')
        set(handles.fineTune_H1,'enable', 'on')
        
        set(handles.text_Hm,'enable', 'off')
        set(handles.edit_Hm,'enable', 'off')
        set(handles.pm_Hm,'enable', 'off')
        set(handles.edit_dHm,'enable', 'off')
        set(handles.fineTune_Hm,'enable', 'off')
        
        set(handles.text_lamda1,'String', 'Lamda 1 (m)')
        set(handles.text_lamda1,'enable', 'on')
        set(handles.edit_lamda1,'enable', 'on')
        set(handles.pm_lamda1,'enable', 'on')
        set(handles.edit_dLamda1,'enable', 'on')
        set(handles.fineTune_lamda1,'enable', 'on')
        
        set(handles.text_lamda2,'String', 'Lamda 2 (m)')
        set(handles.text_lamda2,'enable', 'off')
        set(handles.edit_lamda2,'enable', 'off')
        set(handles.pm_lamda2,'enable', 'off')
        set(handles.edit_dLamda2,'enable', 'off')
        set(handles.fineTune_lamda2,'enable', 'off')
        
        set(handles.edit_tt2,'enable', 'off')
        set(handles.edit_dtt2,'enable', 'off')
        set(handles.pm_tt2,'enable', 'off')
        set(handles.fineTune_tt2,'enable', 'off')
        
        
        set(handles.edit_t3,'enable', 'off')
        set(handles.edit_dt3,'enable', 'off')
        set(handles.pm_t3,'enable', 'off')
        set(handles.fineTune_t3,'enable', 'off')
        
    case 'MTLEI'
        set(handles.text_H1,'enable', 'on')
        set(handles.edit_H1,'enable', 'on')
        set(handles.pm_H1,'enable', 'on')
        set(handles.edit_dH1,'enable', 'on')
        set(handles.fineTune_H1,'enable', 'on')
        
        set(handles.text_Hm,'enable', 'off')
        set(handles.edit_Hm,'enable', 'off')
        set(handles.pm_Hm,'enable', 'off')
        set(handles.edit_dHm,'enable', 'off')
        set(handles.fineTune_Hm,'enable', 'off')
        
        set(handles.text_lamda1,'String', 'Lamda 1 (m)')
        set(handles.text_lamda1,'enable', 'on')
        set(handles.edit_lamda1,'enable', 'on')
        set(handles.pm_lamda1,'enable', 'on')
        set(handles.edit_dLamda1,'enable', 'on')
        set(handles.fineTune_lamda1,'enable', 'on')
        
        set(handles.text_lamda2,'String', 'Lamda 2 (m)')
        set(handles.text_lamda2,'enable', 'off')
        set(handles.edit_lamda2,'enable', 'off')
        set(handles.pm_lamda2,'enable', 'off')
        set(handles.edit_dLamda2,'enable', 'off')
        set(handles.fineTune_lamda2,'enable', 'off')
        
        set(handles.edit_tt2,'enable', 'off')
        set(handles.edit_dtt2,'enable', 'off')
        set(handles.pm_tt2,'enable', 'off')
        set(handles.fineTune_tt2,'enable', 'off')
        
        
        set(handles.edit_t3,'enable', 'off')
        set(handles.edit_dt3,'enable', 'off')
        set(handles.pm_t3,'enable', 'off')
        set(handles.fineTune_t3,'enable', 'off')
        
    case 'MTLK'
        
        set(handles.text_H1,'enable', 'on')
        set(handles.edit_H1,'enable', 'on')
        set(handles.pm_H1,'enable', 'on')
        set(handles.edit_dH1,'enable', 'on')
        set(handles.fineTune_H1,'enable', 'on')
        
        set(handles.text_Hm,'enable', 'off')
        set(handles.edit_Hm,'enable', 'off')
        set(handles.pm_Hm,'enable', 'off')
        set(handles.edit_dHm,'enable', 'off')
        set(handles.fineTune_Hm,'enable', 'off')
        
        set(handles.text_lamda1,'String', 'a')
        set(handles.text_lamda1,'enable', 'on')
        set(handles.edit_lamda1,'enable', 'on')
        set(handles.pm_lamda1,'enable', 'on')
        set(handles.edit_dLamda1,'enable', 'on')
        set(handles.fineTune_lamda1,'enable', 'on')
        
        set(handles.text_lamda2,'String', 'Beeta')
        set(handles.text_lamda2,'enable', 'on')
        set(handles.edit_lamda2,'enable', 'on')
        set(handles.pm_lamda2,'enable', 'on')
        set(handles.edit_dLamda2,'enable', 'on')
        set(handles.fineTune_lamda2,'enable', 'on')
        
        set(handles.edit_tt2,'enable', 'off')
        set(handles.edit_dtt2,'enable', 'off')
        set(handles.pm_tt2,'enable', 'off')
        set(handles.fineTune_tt2,'enable', 'off')
        
        
        set(handles.edit_t3,'enable', 'off')
        set(handles.edit_dt3,'enable', 'off')
        set(handles.pm_t3,'enable', 'off')
        set(handles.fineTune_t3,'enable', 'off')
        
    case 'MTLEL'
        set(handles.text_H1,'enable', 'on')
        set(handles.edit_H1,'enable', 'on')
        set(handles.pm_H1,'enable', 'on')
        set(handles.edit_dH1,'enable', 'on')
        set(handles.fineTune_H1,'enable', 'on')
        
        set(handles.text_Hm,'enable', 'off')
        set(handles.edit_Hm,'enable', 'off')
        set(handles.pm_Hm,'enable', 'off')
        set(handles.edit_dHm,'enable', 'off')
        set(handles.fineTune_Hm,'enable', 'off')
        
        set(handles.text_lamda1,'String', 'Lamda 1 (m)')
        set(handles.text_lamda1,'enable', 'on')
        set(handles.edit_lamda1,'enable', 'on')
        set(handles.pm_lamda1,'enable', 'on')
        set(handles.edit_dLamda1,'enable', 'on')
        set(handles.fineTune_lamda1,'enable', 'on')
        
        set(handles.text_lamda2,'String', 'Lamda 2 (m)')
        set(handles.text_lamda2,'enable', 'off')
        set(handles.edit_lamda2,'enable', 'off')
        set(handles.pm_lamda2,'enable', 'off')
        set(handles.edit_dLamda2,'enable', 'off')
        set(handles.fineTune_lamda2,'enable', 'off')
        
        set(handles.edit_tt2,'enable', 'on')
        set(handles.edit_dtt2,'enable', 'on')
        set(handles.pm_tt2,'enable', 'on')
        set(handles.fineTune_tt2,'enable', 'on')
        
        
        set(handles.edit_t3,'enable', 'on')
        set(handles.edit_dt3,'enable', 'on')
        set(handles.pm_t3,'enable', 'on')
        set(handles.fineTune_t3,'enable', 'on')
        
        
    otherwise
        % Do nothing        
        
end


% --- Executes during object creation, after setting all properties.
function popModalType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popModalType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popCurrType.
function popCurrType_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
handles.d.currType_num = get(hObject,'Value');
handles.d.currType = contents{get(hObject,'Value')}; 

if handles.d.dynamicMod
    handles = model(handles);
end

guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function popCurrType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popCurrType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sns_str_Callback(hObject, eventdata, handles)

sns_str = get(hObject,'String');

if isempty(str2num(sns_str))
    errordlg('Please use comma seperated numbers','Error: Sensor Numbers')
    set(handles.sns_str,'String',handles.d.sns_str)
    return
end

handles.d.sns_str = sns_str;
guidata(hObject, handles);


% hObject    handle to sns_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sns_str as text
%        str2double(get(hObject,'String')) returns contents of sns_str as a double


% --- Executes during object creation, after setting all properties.
function sns_str_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sns_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pm_x0_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.pm_x0,'Value',0.5)
d = handles.d;

if val > 0.5
    d.x0 = d.x0 + d.dx0;
else
    d.x0 = d.x0 - d.dx0;
end

set(handles.edit_x0,'String',d.x0)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pm_x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pm_y0_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.pm_y0,'Value',0.5)
d = handles.d;

if val > 0.5
    d.y0 = d.y0 + d.dy0;
else
    d.y0 = d.y0 - d.dy0;
end

set(handles.edit_y0,'String',d.y0)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pm_y0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pm_H1_Callback(hObject, eventdata, handles)

val = get(hObject,'Value');
set(handles.pm_H1,'Value',0.5)
d = handles.d;

if val > 0.5
    d.H1 = d.H1 + d.dH1;
else
    d.H1 = d.H1 - d.dH1;
end

set(handles.edit_H1,'String',d.H1)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);

% hObject    handle to pm_H1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function pm_H1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_H1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pm_Hm_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.pm_Hm,'Value',0.5)
d = handles.d;

if val > 0.5
    d.Hm = d.Hm + d.dHm;
else
    d.Hm = d.Hm - d.dHm;
end

set(handles.edit_Hm,'String',d.Hm)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function pm_Hm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_Hm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pm_H2_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.pm_H2,'Value',0.5)
d = handles.d;

if val > 0.5
    d.H2 = d.H2 + d.dH2;
else
    d.H2 = d.H2 - d.dH2;
end

set(handles.edit_H2,'String',d.H2)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pm_H2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_H2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pm_t1_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.pm_t1,'Value',0.5)
d = handles.d;

if val > 0.5
    d.t1 = d.t1 + d.dt1;
else
    d.t1 = d.t1 - d.dt1;
end

set(handles.edit_t1,'String',d.t1)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pm_t1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_t1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pm_t2_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.pm_t2,'Value',0.5)
d = handles.d;

if val > 0.5
    d.t2 = d.t2 + d.dt2;
else
    d.t2 = d.t2 - d.dt2;
end

set(handles.edit_t2,'String',d.t2)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pm_t2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_t2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pm_lamda1_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.pm_lamda1,'Value',0.5)
d = handles.d;

if val > 0.5
    d.lamda1 = d.lamda1 + d.dLamda1;
else
    d.lamda1 = d.lamda1 - d.dLamda1;
end

set(handles.edit_lamda1,'String',d.lamda1)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pm_lamda1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_lamda1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pm_lamda2_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.pm_lamda2,'Value',0.5)
d = handles.d;

if val > 0.5
    d.lamda2 = d.lamda2 + d.dLamda2;
else
    d.lamda2 = d.lamda2 - d.dLamda2;
end

set(handles.edit_lamda2,'String',d.lamda2)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pm_lamda2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_lamda2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function pm_AlpInd_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.pm_AlpInd,'Value',0.5)
d = handles.d;

if val > 0.5
    d.AlpInd = d.AlpInd + d.dAlpInd;
else
    d.AlpInd = d.AlpInd - d.dAlpInd;
end

set(handles.edit_AlpInd,'String',d.AlpInd)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pm_AlpInd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_AlpInd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in RecallBest.
function RecallBest_Callback(hObject, eventdata, handles)

if ~isnan(handles.best_set.diff)
    handles.d = handles.best_set.d;
    update_GUI_vals(handles)
    set(handles.TxtMsg,'String','Best set was restored.')
else
    set(handles.TxtMsg,'String','There is nothing to recall.')
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in push_Reset.
function handles = push_Reset_Callback(hObject, eventdata, handles)

% Defult values
d.version = 1;
d.msg = 'Ready...';     d.dir = '';         d.sns_str = '';
d.dt = 0.5;             d.dh = 50;          d.T1 = 50;
d.T2 = 100;             d.modal = 'MTLE';   d.modal_num = 1;
d.curType = 'Positive'; d.currType_num = 1; d.x0 = 0;
d.y0 = 0;               d.H1 = 5000;        d.Hm = 5500;
d.H2 = 6000;            d.t1 = 10;          d.t2 = 25;
d.v  = 1.0;             d.lamda1 = 200;     d.lamda2 = 200;
d.AlpInd = 10;          d.dynamicMod = 1;

d.dx0 = 50;             d.dy0 = 50;         d.dH1 = 50;
d.dH2 = 50;             d.dHm = 50;         d.dt1 = 0.5;
d.dt2 = 0.5;            d.dv  = 1;          d.dLamda1 = 20;
d.dLamda2 = 20;         d.dAlpInd = 0.5;

d.tt2 = 20;             d.dtt2 = 0.5;
d.t3 = 1000;            d.dt3 = 1;

handles.d = d;

% Update GUI values
update_GUI_vals(handles)

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in push_done.
function push_done_Callback(hObject, eventdata, handles)
% hObject    handle to push_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_T1_Callback(hObject, eventdata, handles)
handles.d.T1 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_T1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_T1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_T2_Callback(hObject, eventdata, handles)
handles.d.T2 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_T2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_T2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dv_Callback(hObject, eventdata, handles)
handles.d.dv = str2double(get(hObject,'String'));

if handles.d.dynamicMod
    d1 = handles.d;
    d2 = handles.d;
    d1.v = d1.v - d1.dv;
    d2.v = d2.v + d2.dv;
    handles.d1 = d1;
    handles.d2 = d2;
    handles = modal_fineTune(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_v.
function fineTune_v_Callback(hObject, eventdata, handles)
d1 = handles.d;
d2 = handles.d;
d1.v = d1.v - d1.dv;
d2.v = d2.v + d2.dv;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);


function edit_v_Callback(hObject, eventdata, handles)
handles.d.v = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function pm_v_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.pm_v,'Value',0.5)
d = handles.d;

if val > 0.5
    d.v = d.v + d.dv;
else
    d.v = d.v - d.dv;
end

set(handles.edit_v,'String',d.v)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pm_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function handles = prev_diffs(handles)

if isnan(handles.diffs.CurrDiff2) || ...
   handles.diffs.CurrDiff1 < handles.diffs.CurrDiff2 || ...
   handles.diffs.CurrNormDiff1 < handles.diffs.CurrNormDiff2

    handles.diffs.PrevDiff = handles.diffs.CurrDiff1;
    handles.diffs.PrevNormDiff = handles.diffs.CurrNormDiff1;
    handles.diffs.Prevd = handles.diffs.d1;
    
else
    handles.diffs.PrevDiff = handles.diffs.CurrDiff2;
    handles.diffs.PrevNormDiff = handles.diffs.CurrNormDiff2;
    handles.diffs.Prevd = handles.diffs.d2;
end
    

function handles = model(handles)
tic
set(handles.TxtMsg,'String','Working on...')
try
    drawnow
end

% save previous differences
handles = prev_diffs(handles);

d = handles.d;

% Initial parameters

if d.currType_num == 1
    maxA=1000;
else
    maxA = -1000;
end


sns = str2num(d.sns_str);

%% Starting program
% sensor locations
sns_x =  1.0e+04 *[ -0.7524, 0.3372, -0.3149, -0.4266, -1.1658, -0.2701, -2.7640, ...
                 -6.0091, 0.1825, -5.7394, -2.0637];
             
sns_y =    1.0e+04 *[1.6555, 0.4446, -0.6838, 0.9545, -1.7020, 0.2631, 4.9254, ...
                    -3.3983, -5.3008, 1.1923, 0.1569];

sns_IDs = {'K02','K14','K24','WSB','BCC','K17','EDW','STC','FLT','OVD'};


sns_x = sns_x(sns);
sns_y = sns_y(sns);
sns_IDs = sns_IDs(sns);

% Load real data
L0 = length(sns);

tmin = inf;
tleng = inf;

for i = 1:L0
    rd(i) = open([d.dir sns_IDs{i} '.mat']);
    [tmp1 tmpind] = min(rd(i).y);
    tmp1 = rd(i).t(tmpind);
    tmp2 = range(rd(i).t);
    
    if tmp1 < tmin;      tmin = tmp1;    end    
    if tmp2 < tleng;     tleng = tmp2;  end
end

% Round tmin to lowest milisecond
tleng = tleng*1.2;
tmin = round(tmin*1e6)/1e6;
tleng = round(tleng*2e5)/2e5;

handles.tmin = tmin;
handles.tleng = tleng;

try temp = d.tt2;
catch
    d.tt2 = 20;
    d.dtt2 = .1;
    d.t3 = 1000;
    d.dt3 = 1;
end
    
    
   


% Generate data
cdd = generate_data(maxA,d.t1*1e-6,d.t2*1e-6,d.tt2*1e-6,d.t3*1e-6,d.v*1e7,d.H1,d.Hm,d.H2,d.AlpInd,d.x0,d.y0,...
    d.dt*1e-6,d.lamda1,d.lamda2,sns_x,sns_y,d.dh,rd,L0,d.T1*1e-6,d.T2*1e-6,d.modal);


% Get fit values
[diff, norm_diff, scale] = compare_data(rd, L0,cdd);

toffset = sec2hhmmss(tmin);

switch d.modal
    case 'MTLE'
        tit_str = sprintf('Manual Para search (%s)\n---------------------------\nT-offset = %s UT\nT-duration = %0.0f \\mus \nt-step = % 6.1f \\mus \ndh     = % 6i m\n\nmaxA   = % 6.1f kA \nCharge Mom   = % 6.2f C m \nH1     = % 6.0f m\nH2     = % 6.0f m\nx0     = % 6.0f m\ny0     = % 6.0f m\nt1     = % 6.1f \\mus \nt2     = % 6.1f \\mus \nv      = % 6.1e m/s \nlamda = % 6i m\nalpInd = % 6.1f\n\nFit values...\n\\rho = % 6.2e\n\\rho_{norm} = % 6.2e\n',...
                d.modal,toffset,tleng*1e6,d.dt,d.dh,maxA*scale/1000,cdd(1).P*scale,d.H1,d.H2,d.x0,d.y0,d.t1,d.t2,d.v*1e7,d.lamda1,d.AlpInd,diff,norm_diff);
    case 'MTLEI'
      if d.H1 < d.H2
            PeakCurrNew = maxA*scale/1000*exp((d.H2-d.H1)/d.lamda1);
       else
            PeakCurrNew = maxA*scale/1000*exp((d.H1-d.H2)/d.lamda1);
       end
       
        tit_str = sprintf('Manual Para search (%s) \n---------------------------\nT-offset = %s UT\nT-duration = %0.0f \\mus\nt-step = % 6.1f \\mus \ndh     = % 6i m \n\nA   = % 4.1e kA \nPeakCurr   = % 6.1f kA \nCharge Mom   = % 6.2f C m \nH1     = % 6.0f m \nH2     = % 6.0f m \nx0     = % 6.0f m \ny0     = % 6.0f m \nt1     = % 6.1f \\mus \nt2     = % 6.1f \\mus \nv      = % 6.1e m/s \nlamda = % 6i m \nalpInd = % 6.1f\n\nFit values...\n\\rho = % 6.2e\n\\rho_{norm} = % 6.2e\n',...
                d.modal,toffset,tleng*1e6,d.dt,d.dh,maxA*scale/1000,PeakCurrNew,cdd(1).P*scale,d.H1,d.H2,d.x0,d.y0,d.t1,d.t2,d.v*1e7,d.lamda1,d.AlpInd,diff,norm_diff);
    case 'MTLL'
        tit_str = sprintf('Manual Para search (%s) \n---------------------------\nT-offset = %s UT\nT-duration = %0.0f \\mus\nt-step = % 6.1f \\mus \ndh     = % 6i m \n\nmaxA   = % 6.1f kA \nCharge Mom   = % 6.2f C m\nH1     = % 6.0f m \nH2     = % 6.0f m \nx0     = % 6.0f m \ny0     = % 6.0f m \nt1     = % 6.1f \\mus \nt2     = % 6.1f \\mus \nv      = % 6.1e m/s \nalpInd = % 6.1f\n\nFit values...\n\\rho = % 6.2e\n\\rho_{norm} = % 6.2e\n',...
                d.modal,toffset,tleng*1e6,d.dt,d.dh,maxA*scale/1000,cdd(1).P*scale,d.H1,d.H2,d.x0,d.y0,d.t1,d.t2,d.v*1e7,d.AlpInd,diff,norm_diff);
    case 'MTLEID'
        tit_str = sprintf('Manual Para search (%s) \n---------------------------\nT-offset = %s UT\nT-duration = %0.0f \\mus\nt-step = % 6.1f \\mus \ndh     = % 6i m \nmaxA   = % 6.1f kA \nH1     = % 6.0f\nHm     = % 6.0f\nH2     = % 6.0f\nx0     = % 6.0f\ny0     = % 6.0f\nt1     = % 6.1f\nt2     = % 6.1f\nv      = % 6.1e\nlamda1 = % 6i\nlamda2 = % 6i\nalpInd = % 6.1f\n\nFit values...\n\\rho = % 6.2e\nNdiff = % 6.2e\n',...
                d.modal,toffset,tleng*1e6,d.dt,d.dh,maxA*scale/1000,d.H1,d.Hm,d.H2,d.x0,d.y0,d.t1,d.t2,d.v*1e7,d.lamda1,d.lamda2,d.AlpInd,diff,norm_diff);
    case 'MTLK'
        set(handles.edit_Hm,'String',sprintf('%.0f',cdd.Hm))
        d.Hm = cdd.Hm;
        tit_str = sprintf('Manual Para search (%s) \n---------------------------\nT-offset = %s UT\nT-duration = %0.0f \\mus\nt-step = % 6.1f \\mus \ndh     = % 6i m \nmaxA   = % 6.1f kA \nCharge Mom   = % 6.2f C m\nH1     = % 6.0f m \nHm     = % 6.0f m \nH2     = % 6.0f m \nx0     = % 6.0f m \ny0     = % 6.0f m \nt1     = % 6.1f \\mus \nt2     = % 6.1f \\mus \nv      = % 6.1e m/s \nalpInd = % 6.1f\na    = % 6.1f\nb    = % 6.1f\nFit values...\n\\rho = % 6.2e\n\\rho_{norm} = % 6.2e\n',...
                d.modal,toffset,tleng*1e6,d.dt,d.dh,maxA*scale/1000,cdd(1).P*scale,d.H1,d.Hm,d.H2,d.x0,d.y0,d.t1,d.t2,d.v*1e7,d.AlpInd,d.lamda1,d.lamda2,diff,norm_diff);
    case 'MTLEL'
            tit_str = sprintf('Manual Para search (%s) \n---------------------------\nT-offset = %s UT\nT-duration = %0.0f \\mus\nt-step = % 6.1f \\mus \ndh     = % 6i m \n\nA   = % 4.1e kA \nCharge Mom   = % 6.2f C m \nH_1     = % 6.0f m \nH_2     = % 6.0f m \nx_0     = % 6.0f m \ny_0     = % 6.0f m \nt_1     = % 6.1f \\mus \nt_2     = % 6.1f \\mus \nt^''_2     = % 6.1f \\mus \nt_3     = % 6.1f \\mus \nv      = % 6.1e m/s \n\\lambda = % 6i m \n\\alpha_{Ind} = % 6.1f\n\nFit values...\n\\rho = % 6.2e\n\\rho_{norm} = % 6.2e\n',...
                d.modal,toffset,tleng*1e6,d.dt,d.dh,maxA*scale/1000,cdd(1).P*scale,d.H1,d.H2,d.x0,d.y0,d.t1,d.t2,d.tt2,d.t3,d.v*1e7,d.lamda1,d.AlpInd,diff,norm_diff);
end

% Let's save current results in to d
d.diff      = diff;
d.norm_diff = norm_diff;
d.scalse    = scale;
d.maxA      = maxA*scale;
d.P         = cdd(1).P*scale;

handles.d = d;

% Plot figures
handles = plot_figures(rd,cdd,sns_IDs,tit_str,handles);

handles.diffs.CurrDiff1 = diff;
handles.diffs.CurrNormDiff1 = norm_diff;
handles.diffs.d1 = d;
handles.diffs.d2 = [];

% Display difference values and find best values
handles = disp_diff_values(handles);
toc




       
       
%% Functions       
function data = generate_data(maxA,t1,t2,tt2,t3,v,H1,Hm,H2,alpInd,x0,y0,...
           t_step,lamda,lamda2,sns_x,sns_y,dh,rd,L0,T1,T2,method)
       
   clc
   scale = zeros(1,L0);
   
   % if it is MTLK we need to estimate Hm first
   if strcmp(method,'MTLK')       
       x = 0:0.0001:1;
       I = x.^(lamda-1).*((1-x.^lamda).^(lamda2-1));
       I = I/max(I);
       [~, ind] = max(I);
       
       if H2 < H1
           Hm = H1 + (H2-H1)/(length(x)-1)*ind;
       else
           Hm = H2 - (H1-H2)/(length(x)-1)*ind;
       end
       
       fprintf('Hm = %0.f\n',Hm)
       
       %H1 = H2 - (Hm-H2)/(x(ind) - 1);
       data.x = x;
       data.I = I;
       data.z = H1 + (H2-H1)*x;
       data.Hm = Hm;
       
       % % Another way to find H1
       % H1 = (Hm - H2*((lamda - 1)/(lamda*lamda2 - 1))^(1/lamda))/...
       %        (1  -  ((lamda - 1)/(lamda*lamda2 - 1))^(1/lamda))
   end
  
   for j=1:L0
       t = [];
       E_stat = [];
       E_ind = [];
       E_rad = [];
       P = [];
       
       % Changing parameters
       alpha = alpInd/t2;
       
       % Pulse arrival time to the sensor location
       t0 = sqrt((x0-sns_x(j))^2+(y0-sns_y(j))^2+H1^2)/3e8;
       rdt = rd(j).t;
       rdy = rd(j).y;
       
       switch method
           case 'MTLE'
               [t,E_stat,E_ind,E_rad,P] = IBP_mod_MTLE(maxA,t1,t2,v,H1,H2,x0,y0,...
                   t0 - T1, t0 + T2 ,t_step,lamda,sns_x(j),sns_y(j),dh,alpha);
           case 'MTLL'
               [t,E_stat,E_ind,E_rad,P] = IBP_mod_MTLL(maxA,t1,t2,v,H1,H2,x0,y0,...
                   t0 - T1, t0 + T2 ,t_step,lamda,sns_x(j),sns_y(j),dh,alpha);

            case 'MTLEI'
               [t,E_stat,E_ind,E_rad,P] = IBP_mod_MTLEI(maxA,t1,t2,v,H1,H2,x0,y0,...
                   t0 - T1, t0 + T2 ,t_step,lamda,sns_x(j),sns_y(j),dh,alpha);
           case 'MTLEID'
               [t,E_stat,E_ind,E_rad] = IBP_mod_MTLEID(maxA,t1,t2,v,H1,Hm,H2,x0,y0,...
                   t0 - T1, t0 + T2 ,t_step,lamda,lamda2,sns_x(j),sns_y(j),dh,alpha); 
           case 'MTLK'
               [t,E_stat,E_ind,E_rad,P] = IBP_mod_MTLK(maxA,t1,t2,v,H1,Hm,H2,x0,y0,...
                   t0 - T1, t0 + T2 ,t_step,lamda,lamda2,sns_x(j),sns_y(j),dh,alpha);
           case 'MTLEL'
               [t,E_stat,E_ind,E_rad,P,Is,didts] = IBP_mod_MTLEL(maxA,t1,t2,tt2,t3,v,H1,H2,x0,y0,...
                   t0 - T1, t0 + T2 ,t_step,lamda,sns_x(j),sns_y(j),dh,alpha);
                data(j).Is = Is;
                data(j).didts = didts;
           otherwise
               disp('This model is not working yet')                  
       end
       
       E_tot = E_stat+E_ind+E_rad;
       
%        figure
%        hold all
%        plot(rdy)
%        plot(E_tot)
       
       [hpk1, dt1] = max(E_tot);
       [hpk2, dt2] = max(rdy);
       
       [lpk1, dt3] = min(E_tot);
       [lpk2, dt4] = min(rdy);
   
       %scale(j) = mean([hpk2/hpk1,lpk2/lpk1]);
       scale(j) = (hpk2-lpk2)/(hpk1-lpk1);
      
       
       % time shift
       if hpk2 > -lpk2
            data(j).t = t + rdt(dt2) - t(dt1);
       else
            data(j).t = t + rdt(dt4) - t(dt3);
       end        
       
       % save data
       data(j).y = E_tot;  
       data(j).P = P;

   end
   
   data(1).scale = mean(scale);
   
   
function handles = plot_figures(rd,cdd,sns_IDs,tit_str,handles)

L = length(sns_IDs);

%add_static = [0.1728 -0.8462 -4.0494 -0.7778 0 0 0 0 0];
add_static = zeros(1:10);

tmin  = handles.tmin;
tleng = handles.tleng;

d = handles.d;


if strcmp(d.modal,'MTLK') || strcmp(d.modal,'MTLEL')
    L = L + 1;
end

if     L == 1; m = 1; n = 2; 
elseif L == 2; m = 1; n = 3;
elseif L == 3; m = 2; n = 2;
elseif L <= 5; m = 2; n = 3;
elseif L <= 7; m = 2; n = 4;
elseif L <= 8; m = 3; n = 3;
else m = 3; n = 4;
end


% check figure exist?
fg = findobj('type','figure','name','Pulse Modeler');

if isempty(fg)
    fg = figure(1431);
    sz = get(0, 'ScreenSize');
    set(fg, 'Position', [0 0 sz(3)*0.9 sz(4)*0.9 ],'name','Pulse Modeler')
    tools2fig
else
    clf
end

% Plotting I vs z
if strcmp(d.modal,'MTLK') 
    L = L - 1;
    subplot(m,n,L+1); hold all; box on;
    plot(cdd(1).z/1000,cdd(1).I)

    xlabel('Altitude (km)')
    ylabel('Normalized peak Current')
    %[Hm, ind] = max(I);
    %Hm = z(ind);
    %handles.d.Hm = Hm;
    %update_GUI_vals(handles);   
    
end


% Plotting I vs z
if strcmp(d.modal,'MTLEL') 
    L = L - 1;
    
    [mm, mind] = min(cdd(1).Is);
    tshift_ii = cdd(1).t(mind) - tmin;
    
    subplot(m,n,L+1); hold all; box on;
    plot((cdd(1).t-tmin-tshift_ii)*1e6,cdd(1).Is/1000)

    xlabel('Time (us)')
    ylabel('Current (A)')
    %[Hm, ind] = max(I);
    %Hm = z(ind);
    %handles.d.Hm = Hm;
    %update_GUI_vals(handles);   
    
end

sn_nums = str2num(d.sns_str);

for i = 1:L
    subplot(m,n,i); hold all; box on;
    plot((rd(i).t-tmin)*1e6,rd(i).y);
    plot((cdd(i).t-tmin)*1e6,cdd(i).y*cdd(1).scale + add_static(i));
    [mmtemp mmind] = min(rd(i).y);    
    %xlim(([min(rd(i).t) max(rd(i).t)]-tmin)*1e6)
    %xlim(([0 tleng]+min(rd(i).t)-tmin)*1e6)
    %xlim(([-tleng/2 tleng/2]+rd(i).t(mmind)-tmin)*1e6)
    
    % Xlimit
    if ~d.auto_ylim
        l_ylim = str2num(d.lower_ylim);
        u_ylim = str2num(d.upper_ylim);
        
        if length(d.lower_ylim) == 1
            yl1 = l_ylim;
            yl2 = u_ylim;            
        else
            yl1 = l_ylim(i);
            yl2 = u_ylim(i);
        end
        
        ylim([yl1 yl2])
            
    end
   
    
    xlim(([cdd(i).t(1) cdd(i).t(end)]-tmin)*1e6)
    
    D = round(sqrt((d.x0-handles.sen_set.x(sn_nums(i))).^2 ...
                   +(d.y0-handles.sen_set.y(sn_nums(i))).^2)/1000);
    
    title([sns_IDs{i} '     D = ' num2str(D) ' km'])
    %legend('Real','Calc','Location','northwest')
    legend('Real','Calc')
    plot(xlim,[0 0],'r-.')
end


subplot(m,n,m*n); box off; axis off;
text(0.1,0.5,tit_str) 
set(gca,'XTick',[])
set(gca,'YTick',[])

%title([rd(1).rd '    ' sec2hhmmss(rd(1).t(1))])
set(fg,'visible','on')

function [diff, norm_diff, scale] = compare_data(rd, L0,data)
  
   diffs = nan(1,L0);
   norm_diffs = nan(1,L0);
   scale = data(1).scale;
   
   for j = 1:L0
        
       % Get calculated data for real data
       caly = interp1(data(j).t,data(j).y,rd(j).t)*scale;
       
       % differnce
       diffs(j) =  sqrt(nanmean((caly - rd(j).y).^2));
       
       factor = 1/range(rd(j).y);
       
       norm_diffs(j) = sqrt(nanmean((caly - rd(j).y).^2))*factor;

   end
   
   diff = mean(diffs);
   norm_diff = mean(norm_diffs);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
    d = handles.d;
    save('pulse_modeler_last_gui_data.mat','-Struct','d')
    
delete(hObject);

function handles = modal_fineTune(handles)
tic
set(handles.TxtMsg,'String','Working on...')
drawnow
d = handles.d1;

% save previous differences
handles = prev_diffs(handles);

% Initial parameters
maxA=1000;
sns = str2num(d.sns_str);

%% Starting program
% sensor locations
sns_x =  1.0e+04 *[ -0.7524, 0.3372, -0.3149, -0.4266, -1.1658, -0.2701, -2.7640, ...
                 -6.0091, 0.1825, -5.7394, -2.0637];
             
sns_y =    1.0e+04 *[1.6555, 0.4446, -0.6838, 0.9545, -1.7020, 0.2631, 4.9254, ...
                    -3.3983, -5.3008, 1.1923, 0.1569];

sns_IDs = {'K02','K14','K24','WSB','BCC','K17','EDW','STC','FLT','OVD'};


sns_x = sns_x(sns);
sns_y = sns_y(sns);
sns_IDs = sns_IDs(sns);

% Load real data
L0 = length(sns);

for i = 1:L0
    rd(i) = open([d.dir sns_IDs{i} '.mat']);
end

% Generate data
cdd1 = generate_data(maxA,d.t1*1e-6,d.t2*1e-6,d.v*1e7,d.H1,d.Hm,d.H2,d.AlpInd,d.x0,d.y0,...
    d.dt*1e-6,d.lamda1,d.lamda2,sns_x,sns_y,d.dh,rd,L0,d.T1*1e-6,d.T2*1e-6,d.modal);

% Get fit values
[diff1, norm_diff1, scale1] = compare_data(rd, L0,cdd1);

d = handles.d2;

% Generate data
cdd2 = generate_data(maxA,d.t1*1e-6,d.t2*1e-6,d.v*1e7,d.H1,d.Hm,d.H2,d.AlpInd,d.x0,d.y0,...
    d.dt*1e-6,d.lamda1,d.lamda2,sns_x,sns_y,d.dh,rd,L0,d.T1*1e-6,d.T2*1e-6,d.modal);

% Get fit values
[diff2, norm_diff2, scale2] = compare_data(rd, L0,cdd2);


handles.diffs.CurrDiff1 = diff1;
handles.diffs.CurrNormDiff1 = norm_diff1;
handles.diffs.CurrDiff2 = diff2;
handles.diffs.CurrNormDiff2 = norm_diff2;
handles.diffs.d1 = handles.d1;
handles.diffs.d2 = handles.d2;

handles = disp_diff_values(handles);


% witch one is better?
if norm_diff1 < handles.diffs.PrevNormDiff && norm_diff1 < norm_diff2
    diff = diff1;
    norm_diff = norm_diff1;
    scale = scale1;
    cdd = cdd1;
    d = handles.d1;
	handles.d = handles.d1;
	update_GUI_vals(handles);
	
elseif norm_diff2 < handles.diffs.PrevNormDiff && norm_diff2 < norm_diff1
    diff = diff2;
    norm_diff = norm_diff2;
    scale = scale2;
    cdd = cdd2;
	handles.d = handles.d2;
    d = handles.d2;
	update_GUI_vals(handles);
else
    toc
    return
end
    
switch d.modal
    case 'MTLE'
        tit_str = sprintf('Manual Para search (%s) \nmaxA   = % 6.1f kA \nCharge Mom   = % 6.2f C m \nH1     = % 6.0f\nH2     = % 6.0f\nx0     = % 6.0f\ny0     = % 6.0f\nt-step = % 6.1f\ndh     = % 6i\nt1     = % 6.1f\nt2     = % 6.1f\nv      = % 6.1e\nlamda = % 6i\nalpInd = % 6.1f\n\nFit values...\ndiff = % 6.2e\nNdiff = % 6.2e\n',...
                d.modal,maxA*scale/1000,cdd(1).P*scale,d.H1,d.H2,d.x0,d.y0,d.dt,d.dh,d.t1,d.t2,d.v*1e7,d.lamda1,d.AlpInd,diff,norm_diff);
    case 'MTLEI'
        if d.H1 > d.H2
            PeakCurrNew = maxA*scale/1000*exp((d.H2-d.H1)/d.lamda1);
        else
            PeakCurrNew = maxA*scale/1000*exp((d.H1-d.H2)/d.lamda1);
        end
        PeakCurrNew
        tit_str = sprintf('Manual Para search (%s) \nA   = % 6.1f kA \nPeakCurr   = % 6.1f kA \nCharge Mom   = % 6.2f C m \nH1     = % 6.0f\nH2     = % 6.0f\nx0     = % 6.0f\ny0     = % 6.0f\nt-step = % 6.1f\ndh     = % 6i\nt1     = % 6.1f\nt2     = % 6.1f\nv      = % 6.1e\nlamda = % 6i\nalpInd = % 6.1f\n\nFit values...\ndiff = % 6.2e\nNdiff = % 6.2e\n',...
                d.modal,maxA*scale/1000,PeakCurrNew,cdd(1).P*scale,d.H1,d.H2,d.x0,d.y0,d.dt,d.dh,d.t1,d.t2,d.v*1e7,d.lamda1,d.AlpInd,diff,norm_diff);
    case 'MTLL'
        tit_str = sprintf('Manual Para search (%s) \nmaxA   = % 6.1f kA \nCharge Mom   = % 6.2f C m\nH1     = % 6.0f\nH2     = % 6.0f\nx0     = % 6.0f\ny0     = % 6.0f\nt-step = % 6.1f\ndh     = % 6i\nt1     = % 6.1f\nt2     = % 6.1f\nv      = % 6.1e\nalpInd = % 6.1f\n\nFit values...\ndiff = % 6.2e\nNdiff = % 6.2e\n',...
                d.modal,maxA*scale/1000,cdd(1).P*scale,d.H1,d.H2,d.x0,d.y0,d.dt,d.dh,d.t1,d.t2,d.v*1e7,d.AlpInd,diff,norm_diff);
    case 'MTLEID'
        tit_str = sprintf('Manual Para search (%s) \nmaxA   = % 6.1f kA\nH1     = % 6.0f\nHm     = % 6.0f\nH2     = % 6.0f\nx0     = % 6.0f\ny0     = % 6.0f\nt-step = % 6.1f\ndh     = % 6i\nt1     = % 6.1f\nt2     = % 6.1f\nv      = % 6.1e\nlamda1 = % 6i\nlamda2 = % 6i\nalpInd = % 6.1f\n\nFit values...\ndiff = % 6.2e\nNdiff = % 6.2e\n',...
                d.modal,maxA*scale/1000,d.H1,d.Hm,d.H2,d.x0,d.y0,d.dt,d.dh,d.t1,d.t2,d.v*1e7,d.lamda1,d.lamda2,d.AlpInd,diff,norm_diff);
    case 'MTLK'
        d.H1 = cdd(1).z(1);
        tit_str = sprintf('Manual Para search (%s) \nmaxA   = % 6.1f kA \nCharge Mom   = % 6.2f C m \nH1     = % 6.0f\nHm     = % 6.0f\nH2     = % 6.0f\nx0     = % 6.0f\ny0     = % 6.0f\nt-step = % 6.1f\ndh     = % 6i\nt1     = % 6.1f\nt2     = % 6.1f\nv      = % 6.1e\nalpInd = % 6.1f\na    = % 6.1f\nb    = % 6.1f\n\nFit values...\ndiff = % 6.2e\nNdiff = % 6.2e\n',...
                d.modal,maxA*scale/1000,cdd(1).P*scale,d.H1,d.Hm,d.H2,d.x0,d.y0,d.dt,d.dh,d.t1,d.t2,d.v*1e7,d.AlpInd,d.lamda1,d.lamda2,diff,norm_diff);
            
end


% Let's save current results in to d
d.diff      = diff;
d.norm_diff = norm_diff;
d.scalse    = scale;
d.maxA      = maxA*scale;
d.P         = cdd(1).P*scale;

handles.d = d;
% Plot figures
handles = plot_figures(rd,cdd,sns_IDs,tit_str,handles);
toc




function handles = disp_diff_values(handles)

% display differences
set(handles.CurrDiff1,'String',handles.diffs.CurrDiff1)
set(handles.CurrNormDiff1,'String',handles.diffs.CurrNormDiff1)
set(handles.CurrDiff2,'String',handles.diffs.CurrDiff2)
set(handles.CurrNormDiff2,'String',handles.diffs.CurrNormDiff2)
set(handles.PrevDiff,'String',handles.diffs.PrevDiff)
set(handles.PrevNormDiff,'String',handles.diffs.PrevNormDiff)
set(handles.BestDiff,'String',handles.diffs.BestDiff)
set(handles.BestNormDiff,'String',handles.diffs.BestNormDiff)

if isnan(handles.diffs.CurrDiff2)
    diff = handles.diffs.CurrDiff1;
    NormDiff = handles.diffs.CurrNormDiff1;
    str = 'The change ';
    d = handles.diffs.d1;
elseif  handles.diffs.CurrNormDiff1 < handles.diffs.CurrNormDiff2
    diff = handles.diffs.CurrDiff1;
    NormDiff = handles.diffs.CurrNormDiff1;
    str = 'The reduction ';
    d = handles.diffs.d1;
else
    diff = handles.diffs.CurrDiff2;
    NormDiff = handles.diffs.CurrNormDiff2;
    str = 'The increase ';
    d = handles.diffs.d2;
end
    
    
% Is the current set better?
if isnan(handles.diffs.PrevDiff)
    set(handles.TxtMsg,'String','Done!')
elseif NormDiff < handles.diffs.PrevNormDiff
    set(handles.TxtMsg,'String',['Done! ' str 'made it better.'])
elseif NormDiff > handles.diffs.PrevNormDiff
   set(handles.TxtMsg,'String',['Done! ' str 'made it worst.'])    
else
    set(handles.TxtMsg,'String',['Done! ' str 'made no difference.'])
end
    
% Find the best set
if isnan(handles.best_set.diff) || isnan(handles.best_set.norm_diff)
    handles.best_set.d = d;
    handles.best_set.diff = handles.diffs.CurrDiff1;
    handles.best_set.norm_diff = handles.diffs.CurrNormDiff1;
    handles.diffs.BestDiff = handles.diffs.CurrDiff1;
    handles.diffs.BestNormDiff = handles.diffs.CurrNormDiff1;
    set(handles.BestDiff,'String',handles.diffs.BestDiff)
    set(handles.BestNormDiff,'String',handles.diffs.BestNormDiff)
    
elseif handles.best_set.norm_diff > NormDiff
    
    set(handles.TxtMsg,'String',['Done! '  str ' made it better than the best.'])
    handles.best_set.d = d;
    handles.best_set.diff = diff;
    handles.best_set.norm_diff = NormDiff;
    
    handles.diffs.BestDiff = diff;
    handles.diffs.BestNormDiff = NormDiff;
    set(handles.BestDiff,'String',handles.diffs.BestDiff)
    set(handles.BestNormDiff,'String',handles.diffs.BestNormDiff)
    
end


function PrevDiff_Callback(hObject, eventdata, handles)
% hObject    handle to PrevDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PrevDiff as text
%        str2double(get(hObject,'String')) returns contents of PrevDiff as a double


% --- Executes during object creation, after setting all properties.
function PrevDiff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PrevDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BestDiff_Callback(hObject, eventdata, handles)
% hObject    handle to BestDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BestDiff as text
%        str2double(get(hObject,'String')) returns contents of BestDiff as a double


% --- Executes during object creation, after setting all properties.
function BestDiff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BestDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PrevNormDiff_Callback(hObject, eventdata, handles)
% hObject    handle to PrevNormDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PrevNormDiff as text
%        str2double(get(hObject,'String')) returns contents of PrevNormDiff as a double


% --- Executes during object creation, after setting all properties.
function PrevNormDiff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PrevNormDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BestNormDiff_Callback(hObject, eventdata, handles)
% hObject    handle to BestNormDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BestNormDiff as text
%        str2double(get(hObject,'String')) returns contents of BestNormDiff as a double


% --- Executes during object creation, after setting all properties.
function BestNormDiff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BestNormDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CurrDiff1_Callback(hObject, eventdata, handles)
% hObject    handle to CurrDiff1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrDiff1 as text
%        str2double(get(hObject,'String')) returns contents of CurrDiff1 as a double


% --- Executes during object creation, after setting all properties.
function CurrDiff1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrDiff1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CurrDiff2_Callback(hObject, eventdata, handles)
% hObject    handle to CurrDiff2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrDiff2 as text
%        str2double(get(hObject,'String')) returns contents of CurrDiff2 as a double


% --- Executes during object creation, after setting all properties.
function CurrDiff2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrDiff2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CurrNormDiff1_Callback(hObject, eventdata, handles)
% hObject    handle to CurrNormDiff1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrNormDiff1 as text
%        str2double(get(hObject,'String')) returns contents of CurrNormDiff1 as a double


% --- Executes during object creation, after setting all properties.
function CurrNormDiff1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrNormDiff1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CurrNormDiff2_Callback(hObject, eventdata, handles)
% hObject    handle to CurrNormDiff2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrNormDiff2 as text
%        str2double(get(hObject,'String')) returns contents of CurrNormDiff2 as a double


% --- Executes during object creation, after setting all properties.
function CurrNormDiff2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrNormDiff2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function plot_3d_current_Callback(hObject, eventdata, handles)

d = handles.d;
d.maxA = 1000;

% Convert numbers to real numbers
d.t1 = d.t1*1e-6;
d.t2 = d.t2*1e-6;
d.v  = d.v*1e7;
d.dt = d.dt*1e-6;


Z = d.H1:1:d.H2;
T = 0:0.1e-6:d.t2;
iii = nan(length(T),length(Z));

i = 0;

wbh = waitbar(0,'Plase wait...','name','3D Current generator');

L = length(T);

for t=T
    
    i=i+1;   
  
    ii = current(Z,t-(Z-d.H1)/d.v,d);
    
    iii(i,:) = ii;
    
    waitbar(i/L,wbh)
end

delete(wbh)

iii = iii/max(iii(:));
Z = Z/1000;
T = T*1e6;
figure
H = surf(Z,T,iii);
set(H, 'linestyle', 'none');
shading interp;
xlim([min(Z) max(Z)])
ylim([min(T) max(T)])
xlabel('Altitude (km)')
ylabel('Time (\mus)')
zlabel('Normalized current')
tools2fig



% --------------------------------------------------------------------
function linear_charge_density_Callback(hObject, eventdata, handles)


% Choosing files to load
need_files = 1;

d = handles.d;

ds = [];
clc


while need_files
   
    
    [fn,  dn] = uigetfile([handles.d.dir '\*.mat'],'Choose a file');
        
    if isequal(fn,0); 
        need_files = 0;
        break
    end;
    
    d = open([dn fn]);
    
    ds = [d ds];
end


% dn = 'C:\Users\Sumedhe\Desktop\IBP_modeling_2013\CG-20120814_77636_333\';
% fn = {'MTLK--final-nDiff-3.34e-02.mat'
%     'MTLEI-nDiff-5.22e-02.mat'
%       'MTLE-final-nDiff-3.14e-02.mat'
%       'MTLL-final-nDiff-3.28e-02.mat'};
%   
% for i = 1:4
%     d = open([dn fn{i}]);
%     ds = [d ds];
% end
    

zLowest = 0;

for i=1:length(ds)
    if strcmp(ds(i).modal,'MTLL');
        zLowest = ds(i).H1;
    end
end


% setup the figure
fg = figure;
set(fg,'visible','off')
tools2fig
hold all
box on
%ylabel('Altitude (m)')
%xlabel('\rho_L (\muC m^{-1})')

lg = {};

for n = 1:length(ds)
    
    d = ds(n);
    try
        maxA_temp = d.maxA;
    catch
        % Let's use 1kA current but lets put astrict in the legend to say
        % the peak current info not available
        
        d.maxA = 1000;
    end
    
    % Convert numbers to real numbers
    d.t1 = d.t1*1e-6;
    d.t2 = d.t2*1e-6;
    d.v  = d.v*1e7;
    d.dt = d.dt*1e-6;
    
    
    d.dh = d.v*d.dt;
    
        
    Z = d.H1:d.dh/10:d.H2;
    T = 0:d.dt/10:d.t2*2;

    L1 = length(T);
    L2 = length(Z);
    
    iii = nan(L1,L2);
    
    i = 0;
    
    wbh = waitbar(0,'Plase wait...','name','\rho_L finder');
    
    
    for t=T
        
        i=i+1;
        
        ii = current(Z,t-(Z-d.H1)/d.v,d);
        
        iii(i,:) = ii;
        try
            waitbar(i/L1,wbh);
        catch
            disp('User stopped processing')
            return
        end
    end
    
    Q_passed = zeros(1,L2);
    
    for i=1:L2
        Q_passed(i) = trapz(T,iii(:,i));
    end
    
    delete(wbh)

    rho = -diff(Q_passed)/(d.dh/10);
    
    % charge diposited
    Q_dipo = - -diff(Q_passed);
    
    if strcmp(d.modal, 'MTLEI')        
        %Z = Z + zLowest - d.H1;
        %Z = Z + 1100;
        %lg_tmp = sprintf('%s',d.modal,1100); %zLowest - d.H1);
        
        fprintf('Maximum coordinates of MTLEI will go to\n\t%0.1f\t%0.1f\n',...
            -rho(end)*1e3,Z(end))
        
   end
    
    lg_tmp = d.modal;
    
    if strcmp(d.modal, 'MTLK')        
        % Find the minimum coordinates
        fprintf('Minimum coordinates of MTLK will go to\n\t%0.1f\t%0.1f\n',...
            -rho(1)*1e3,Z(1))
        
        % Calculate linier charge density from Hm to H2
        [mm ind_m] = min(abs(Z - d.Hm));
        totC2 = sum(Q_dipo(ind_m:end));
        lcd2 = totC2/(d.H2 - d.Hm);
        
        fprintf('Total & Linear charge density of upper part of MTLK\n\t%0.2e\t\t%0.2e\n',totC2,lcd2)
    end
    
   if strcmp(d.modal, 'MTLE')        
        % Find the minimum coordinates
        fprintf('Minimum coordinates of MTLE will go to\n\t%0.1f\t%0.1f\n',...
            -rho(1)*1e3,Z(1))
    end
    
    
    if ~exist('maxA_temp','var')
        lg_tmp = [lg_tmp '*'];
    end
    
    lg = [lg lg_tmp];
    
    plot(rho*1e3,Z(1:end-1))
    %plot(Q_dipo,Z(1:end-1))
    
    % Total charge
    charges(n).modal = d.modal;
    charges(n).Q_dipo = Q_dipo;
    charges(n).H1 = d.H1;
    charges(n).H2 = d.H2;
    
    
%     % Calculate static E-changes due to charge distribution
%     sen_set = open('sensor_setting.mat');
%     E = zeros(1,10);
%     for iI=1:10
%         % Considering only charge
%         EE1 = sum(1/2/pi/8.85e-12*Q_dipo.*Z(1:end-1)./(Z(1:end-1).^2+(d.x0 - sen_set.x(iI))^2+(d.y0 - sen_set.y(iI))^2).^1.5);        
%         
%         % Considering equal and opasite charge at H1
%         EE2 = sum(1/2/pi/8.85e-12*Q_dipo.*Z(1:end-1)./(Z(1:end-1).^2+(d.x0 - sen_set.x(iI))^2+(d.y0 - sen_set.y(iI))^2).^1.5 + ...
%             - 1/2/pi/8.85e-12*Q_dipo.*d.H1./(d.H1.^2+(d.x0 - sen_set.x(iI))^2+(d.y0 - sen_set.y(iI))^2).^1.5);        
%         
%         
%         fprintf('%s\t%0.2f\t%0.2f\n',sen_set.sen_IDs{iI},EE1,EE2)
%     end
%     
    
end

% Print total charge deposited
fprintf('Modal\tQ-total\trho-mean\n')
for n = 1:length(charges)
    fprintf('%s\t%0.1e\t\t%0.1e\n',...
    charges(n).modal,sum((charges(n).Q_dipo)),...
        sum((charges(n).Q_dipo))/(charges(n).H2-charges(n).H1));
    temp1(n) = sum(charges(n).Q_dipo);
    temp2(n) = sum((charges(n).Q_dipo))/(charges(n).H2-charges(n).H1);    
end
fprintf('Mean\t%0.1e\t\t%0.1e\n\n',mean(temp1),mean(temp2));

% Print total charge deposited
fprintf('Modal\tQ-total\trho-mean\n')
for n = 1:length(charges)
    fprintf('%s\t%0.1e\t\t%0.1e\n',...
    charges(n).modal,sum(abs(charges(n).Q_dipo)),...
        sum(abs(charges(n).Q_dipo))/(charges(n).H2-charges(n).H1));
    temp1(n) = sum(abs(charges(n).Q_dipo));
    temp2(n) = sum(abs(charges(n).Q_dipo))/(charges(n).H2-charges(n).H1);
end
fprintf('Mean\t%0.1e\t\t%0.1e\n\n',mean(temp1),mean(temp2));

legend(lg)


% for IBP-01 xlim([-1 1]) ylim([5200 6000])
% For IBP-02 
xlim([-5 5]); 
ylim([5000 7000]);
% For IBP-03 xlim([-5 7]) ylim([5000 6700])
% For IBP-04 xlim([-5 5]) ylim([5000 6500])
% For IBP-05 xlim([-1 7]) ylim([5000 6500])
% For IBP-06 xlim([-6 8.5]) ylim([5750 7300])


% GET TICKS
X=get(gca,'Xtick');
Y=get(gca,'Ytick');

% GET LABELS
XL=get(gca,'XtickLabel');
YL=get(gca,'YtickLabel');

% GET OFFSETS
Xoff=diff(get(gca,'XLim'))./40;
Yoff=diff(get(gca,'YLim'))./40;

ylimv = ylim;
xlimv = xlim;
% DRAW AXIS LINEs
plot(get(gca,'XLim'),[ylimv(1) ylimv(1)],'k');
plot([0 0],get(gca,'YLim'),'k');

% Plot new ticks  
for i=1:length(X)
    plot([X(i) X(i)],[ylimv(1) ylimv(1)+Yoff],'-k');
end;
for i=1:length(Y)
   plot([-Xoff, Xoff]/2,[Y(i) Y(i)],'-k');
end;

% ADD TICK LABELS
text(X,zeros(size(X))-2.*Yoff+ylimv(1),XL);
text(zeros(size(Y))-3.5.*Xoff,Y,YL);

% ADD LABELS
text(mean(xlimv),-4*Yoff+ylimv(1),'\rho_L (mC m^{-1})')
text(-5*Xoff,mean(ylimv),'Altitude (m)','Rotation',90)

box off;
% axis square;
axis off;
set(fg,'color','w');
set(fg,'visible','on')




function [ii, didt]=current(z,t,a)


ii      = zeros(size(z));
didt    = zeros(size(z));
a.alpha = a.AlpInd/a.t2;
a.k     = (a.t2-a.t1)/a.t1;

for j=1:length(t)
    
    if t(j) <= a.t1
        
        switch a.modal
            case 'MTLEI'
                ii(j)=a.maxA*exp((z(j)-a.H1)/a.lamda1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
                didt(j)=exp((z(j)-a.H1)/a.lamda1)*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
            case 'TL'
                ii(j)=a.maxA*exp(-a.alpha^2*(t(j)-a.t1)^2);
                didt(j)=(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
            case 'MTLL'
                ii(j)=a.maxA*(1 - (z(j)-a.H1)/(a.H2-a.H1))*exp(-a.alpha^2*(t(j)-a.t1)^2);
                didt(j)=(1 - (z(j)-a.H1)/(a.H2-a.H1))*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
            case 'MTLE'
                ii(j)  = a.maxA*exp(-(z(j)-a.H1)/a.lamda1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
                didt(j)= exp(-(z(j)-a.H1)/a.lamda1)*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
                
            case 'MTLK'
                ii(j)=a.maxA*((z(j)-a.H1)/(a.H2-a.H1)).^(a.lamda1-1).*((1-((z(j)-a.H1)/(a.H2-a.H1)).^a.lamda1).^(a.lamda2-1))*exp(-a.alpha^2*(t(j)-a.t1)^2);
                didt(j)=((z(j)-a.H1)/(a.H2-a.H1)).^(a.lamda1-1).*((1-((z(j)-a.H1)/(a.H2-a.H1)).^a.lamda1).^(a.lamda2-1))*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
                
               % ii(j)=a.maxA*get_power(((z(j)-a.H1)/(a.H2-a.H1)),(a.lamda1-1)).*get_power(get_power(1-((z(j)-a.H1)/(a.H2-a.H1)),a.lamda1),(a.lamda2-1))*exp(-a.alpha^2*(t(j)-a.t1)^2);
               % didt(j)=get_power(((z(j)-a.H1)/(a.H2-a.H1)),(a.lamda1-1)).*get_power(get_power(1-((z(j)-a.H1)/(a.H2-a.H1)),a.lamda1),(a.lamda2-1))*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);

        end
        
    else
        switch a.modal
            case 'MTLEI'
                ii(j)=a.maxA*exp((z(j)-a.H1)/a.lamda1)*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                didt(j)=exp((z(j)-a.H1)/a.lamda1)*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            case 'TL'
                ii(j)=a.maxA*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                didt(j)=(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            case 'MTLL'
                ii(j)=a.maxA*(1 - (z(j)-a.H1)/(a.H2-a.H1))*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                didt(j)=(1 - (z(j)-a.H1)/(a.H2-a.H1))*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            case 'MTLE'
                ii(j)=a.maxA*exp(-(z(j)-a.H1)/a.lamda1)*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                didt(j)=exp(-(z(j)-a.H1)/a.lamda1)*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            case 'MTLK'
                ii(j)=a.maxA*((z(j)-a.H1)/(a.H2-a.H1)).^(a.lamda1-1).*((1-((z(j)-a.H1)/(a.H2-a.H1)).^a.lamda1).^(a.lamda2-1))*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                didt(j)=((z(j)-a.H1)/(a.H2-a.H1)).^(a.lamda1-1).*((1-((z(j)-a.H1)/(a.H2-a.H1)).^a.lamda1).^(a.lamda2-1))*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        
                %ii(j)=a.maxA*get_power(((z(j)-a.H1)/(a.H2-a.H1)),(a.lamda1-1)).*get_power(get_power(1-((z(j)-a.H1)/(a.H2-a.H1)),a.lamda1),(a.lamda2-1))*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                %didt(j)=get_power(((z(j)-a.H1)/(a.H2-a.H1)),(a.lamda1-1)).*get_power(get_power(1-((z(j)-a.H1)/(a.H2-a.H1)),a.lamda1),(a.lamda2-1))*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        end
    end    
end

function res = get_power(x,n)

res = exp(n*log(x));



function edit_dtt2_Callback(hObject, eventdata, handles)
handles.d.dtt2 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dtt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dtt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_tt2.
function fineTune_tt2_Callback(hObject, eventdata, handles)
d1 = handles.d;
d2 = handles.d;
d1.tt2 = d1.tt2 - d1.dtt2;
d2.tt2 = d2.tt2 + d2.dtt2;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);



function edit_tt2_Callback(hObject, eventdata, handles)
handles.d.tt2 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_tt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function pm_tt2_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.pm_tt2,'Value',0.5)
d = handles.d;

if val > 0.5
    d.tt2 = d.tt2 + d.dtt2;
else
    d.tt2 = d.tt2 - d.dtt2;
end

set(handles.edit_x0,'String',d.x0)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pm_tt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_tt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_dt3_Callback(hObject, eventdata, handles)
handles.d.dt3 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dt3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dt3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fineTune_t3.
function fineTune_t3_Callback(hObject, eventdata, handles)
d1 = handles.d;
d2 = handles.d;
d1.t3 = d1.t3 - d1.dt3;
d2.t3 = d2.t3 + d2.dt3;
handles.d1 = d1;
handles.d2 = d2;
handles = modal_fineTune(handles);
guidata(hObject, handles);



function edit_t3_Callback(hObject, eventdata, handles)
handles.d.t3 = str2double(get(hObject,'String'));


if handles.d.dynamicMod
    handles = model(handles);
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_t3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_t3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function pm_t3_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.pm_t3,'Value',0.5)
d = handles.d;

if val > 0.5
    d.t3 = d.t3 + d.dt3;
else
    d.t3 = d.t3 - d.dt3;
end

set(handles.edit_x0,'String',d.x0)
handles.d = d;


if handles.d.dynamicMod
    handles = model(handles);
end    

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pm_t3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_t3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function figure_y_limits_Callback(hObject, eventdata, handles)

try
    auto_ylim = handles.d.auto_ylim;
catch
   handles.d.auto_ylim = true;
   handles.d.upper_ylim = 100;
   handles.d.lower_ylim = 100;
   guidata(hObject, handles);
end
  
auto_ylim = boolean(handles.d.auto_ylim);

guidata(hObject, handles);
[a, button] = settingsdlg(...
      'Description', 'If you use a single value, all plots will  use that y limits. To use different values for each plot, sperate ylimits by spaces.',... 
      'title'      , 'Set y limits',...
      {'Auto y limits'; 'auto_ylim'}, [auto_ylim, auto_ylim],...
      {'Upper limit';'upper_ylim'}, handles.d.upper_ylim,...
      {'Lower limit';'lower_ylim'}, handles.d.lower_ylim);
  

if strcmp(button,'ok')
    
   handles.d.auto_ylim = a.auto_ylim;
   handles.d.upper_ylim = a.upper_ylim;
   handles.d.lower_ylim = a.lower_ylim;
   guidata(hObject, handles);
end
    
 
