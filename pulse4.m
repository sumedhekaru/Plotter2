function varargout = pulse4(varargin)
% PULSE4 M-file for pulse4.fig
%      PULSE4, by itself, creates a new PULSE4 or raises the existing
%      singleton*.
%
%      H = PULSE4 returns the handle to a new PULSE4 or the handle to
%      the existing singleton*.
%
%      PULSE4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PULSE4.M with the given input arguments.
%
%      PULSE4('Property','Value',...) creates a new PULSE4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pulse4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pulse4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pulse4

% Last Modified by GUIDE v2.5 25-Oct-2010 08:57:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pulse4_OpeningFcn, ...
                   'gui_OutputFcn',  @pulse4_OutputFcn, ...
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


% --- Executes just before pulse4 is made visible.
function pulse4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pulse4 (see VARARGIN)

% Choose default command line output for pulse4
handles.output = hObject;

% Loading default values (or last values)

try
    fn=sprintf('%spulse_para.mat',tempdir);
    data=open(fn);
    handles.data=data;
catch
    handles.data=struct( ...
        'maxA',         {''}     ,...
        't1',           {''}     ,...
        't2',           {''}     ,...
        'v',            {''}     ,...
        'H1',           {''}     ,...
        'H2',           {''}     ,...
        'x',            {''}     ,...
        'y',            {''}     ,...
        'T1',           {''}     ,...
        'T2',           {''}     ,...
        'lamda',        {''}     ,...
        't_step',       {''}     ,...
        'N',            {''}     ,...
        't_reflect',    {''}     ,...
        'coeff_reflect',{''}     ,...
        't_shift',      {''}     ,...
        'is_saved',     {''}     ,...
        'save_fn',      {''}    ,...
        'file_name',    {''}   ,...
        'plots_on',     {[1 1 1 1]}, ...
        'addi_plots',   {[1 1 1 1]} ...
        );
end
set(handles.max_current,'String',handles.data.maxA)
set(handles.increase_time,'String',handles.data.t1)
set(handles.total_t,'String',handles.data.t2)
set(handles.velocity,'String',handles.data.v)
set(handles.initial_height,'String',handles.data.H1)
set(handles.final_height,'String',handles.data.H2)
set(handles.x,'String',num2str(handles.data.x))
set(handles.y,'String',num2str(handles.data.y))
set(handles.plot_begining_time,'String',handles.data.T1)
set(handles.plot_end_time,'String',handles.data.T2)
set(handles.lamda,'String',handles.data.lamda)
set(handles.time_step,'String',handles.data.t_step)
set(handles.number_of_reflections,'String',handles.data.N)
set(handles.reflection_time,'String',handles.data.t_reflect)
set(handles.reflection_coeff,'String',handles.data.coeff_reflect)
set(handles.manual_time_shift,'String',handles.data.t_shift)
set(handles.radiation_on,'Value',handles.data.plots_on(1))
set(handles.static_on,'Value',handles.data.plots_on(2))
set(handles.Induction_on,'Value',handles.data.plots_on(3))
set(handles.total,'Value',handles.data.plots_on(4))
set(handles.components_on,'Value',handles.data.addi_plots(1))
set(handles.bcurrents_on,'Value',handles.data.addi_plots(2))
set(handles.current_h_on,'Value',handles.data.addi_plots(3))
set(handles.Info_sheet_on,'Value',handles.data.addi_plots(4))

% Use my own close function
set(handles.figure1,'CloseRequestFcn',@closeGUI);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pulse4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pulse4_OutputFcn(hObject, eventdata, handles) 
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
%handles.data

% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,pn] = uigetfile('*.mat');

filename=[pn fn];

if fn==0
    % File did not selected, do nothing!
else
    data=open(filename);
    
    try
        
        handles.data=data;
        
        set(handles.max_current,'String',handles.data.maxA)
        set(handles.increase_time,'String',handles.data.t1)
        set(handles.total_t,'String',handles.data.t2)
        set(handles.velocity,'String',handles.data.v)
        set(handles.initial_height,'String',handles.data.H1)
        set(handles.final_height,'String',handles.data.H2)
        set(handles.x,'String',num2str(handles.data.x))
        set(handles.y,'String',num2str(handles.data.y))
        set(handles.plot_begining_time,'String',handles.data.T1)
        set(handles.plot_end_time,'String',handles.data.T2)
        set(handles.lamda,'String',handles.data.lamda)
        set(handles.time_step,'String',handles.data.t_step)
        set(handles.number_of_reflections,'String',handles.data.N)
        set(handles.reflection_time,'String',handles.data.t_reflect)
        set(handles.reflection_coeff,'String',handles.data.coeff_reflect)
        set(handles.manual_time_shift,'String',handles.data.t_shift)
        set(handles.radiation_on,'Value',handles.data.plots_on(1))
        set(handles.static_on,'Value',handles.data.plots_on(2))
        set(handles.Induction_on,'Value',handles.data.plots_on(3))
        set(handles.total,'Value',handles.data.plots_on(4))
        set(handles.components_on,'Value',handles.data.addi_plots(1))
        set(handles.bcurrents_on,'Value',handles.data.addi_plots(2))
        set(handles.current_h_on,'Value',handles.data.addi_plots(3))
        set(handles.Info_sheet_on,'Value',handles.data.addi_plots(4))
        
        % If user have renamed file (use the new name)
        handles.data.file_name=filename;
        
        set(handles.figure1,'Name',['Pulse4 -' fn(1:end-4)])
        handles.data.is_saved=1;
        
    catch
        errordlg('Content could not recognise!','File Open Error','modal')
    end
    
end
guidata(hObject, handles)


% --------------------------------------------------------------------
function Save_me_Callback(hObject, eventdata, handles)
% hObject    handle to Save_me (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.data.file_name,'')==1
    [save_fn,save_pn] = uiputfile('*.mat','Save');
    if save_fn(1)~=0
        fn=[save_pn save_fn];
        data=handles.data;
        save(fn,'-struct', 'data')
        handles.data.file_name=fn;
        set(handles.figure1,'Name',['Pulse4 -' save_fn(1:end-4)])
        handles.data.is_saved=1;
    end
    
else
    data=handles.data;
    save(handles.data.file_name,'-struct', 'data')
    handles.data.is_saved=1;
end

guidata(hObject, handles);


% --------------------------------------------------------------------
function save_as_Callback(hObject, eventdata, handles)
% hObject    handle to save_as (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    [save_fn,save_pn] = uiputfile('*.mat','Save As..');
    if save_fn(1)~=0
        fn=[save_pn save_fn];
        data=handles.data;
        save(fn,'-struct', 'data')
        handles.data.file_name=fn;
        set(handles.figure1,'Name',['Pulse4 -' save_fn(1:end-4)])
        g.isSaved=1;
    end



% --------------------------------------------------------------------
function quit_me_Callback(hObject, eventdata, handles)
% hObject    handle to quit_me (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
closeGUI(hObject,eventdata)


function max_current_Callback(hObject, eventdata, handles)
% hObject    handle to max_current (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.maxA=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function max_current_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_current (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%set(hObject,'String',handles.data.maxA)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function increase_time_Callback(hObject, eventdata, handles)
% hObject    handle to increase_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.t1=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function increase_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to increase_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function velocity_Callback(hObject, eventdata, handles)
% hObject    handle to velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.v=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function velocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initial_height_Callback(hObject, eventdata, handles)
% hObject    handle to initial_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.H1=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);        


% --- Executes during object creation, after setting all properties.
function initial_height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initial_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function final_height_Callback(hObject, eventdata, handles)
% hObject    handle to final_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.H2=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function final_height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to final_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_Callback(hObject, eventdata, handles)
% hObject    handle to x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str=regexp(get(hObject,'String'), ' ', 'split');
handles.data.x=str2double(str);
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_Callback(hObject, eventdata, handles)
% hObject    handle to y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str=regexp(get(hObject,'String'), ' ', 'split');
handles.data.y=str2double(str);
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_begining_time_Callback(hObject, eventdata, handles)
% hObject    handle to plot_begining_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.T1=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function plot_begining_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_begining_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plot_end_time_Callback(hObject, eventdata, handles)
% hObject    handle to plot_end_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.T2=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function plot_end_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_end_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lamda_Callback(hObject, eventdata, handles)
% hObject    handle to lamda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.lamda=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function lamda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lamda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function time_step_Callback(hObject, eventdata, handles)
% hObject    handle to time_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.t_step=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function time_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_of_reflections_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_reflections (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.N=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function number_of_reflections_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_reflections (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function reflection_time_Callback(hObject, eventdata, handles)
% hObject    handle to reflection_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.t_reflect=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function reflection_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reflection_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function reflection_coeff_Callback(hObject, eventdata, handles)
% hObject    handle to reflection_coeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.coeff_reflect=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function reflection_coeff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reflection_coeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function manual_time_shift_Callback(hObject, eventdata, handles)
% hObject    handle to manual_time_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.t_shift=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function manual_time_shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to manual_time_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in components_on.
function components_on_Callback(hObject, eventdata, handles)
% hObject    handle to components_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.addi_plots(1)=get(hObject,'Value');
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in bcurrents_on.
function bcurrents_on_Callback(hObject, eventdata, handles)
% hObject    handle to bcurrents_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.addi_plots(2)=get(hObject,'Value');
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in current_h_on.
function current_h_on_Callback(hObject, eventdata, handles)
% hObject    handle to current_h_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.addi_plots(3)=get(hObject,'Value');
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Info_sheet_on.
function Info_sheet_on_Callback(hObject, eventdata, handles)
% hObject    handle to Info_sheet_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.addi_plots(4)=get(hObject,'Value');
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in radiation_on.
function radiation_on_Callback(hObject, eventdata, handles)
% hObject    handle to radiation_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.plots_on(1)=get(hObject,'Value');
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in static_on.
function static_on_Callback(hObject, eventdata, handles)
% hObject    handle to static_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.plots_on(2)=get(hObject,'Value');
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Induction_on.
function Induction_on_Callback(hObject, eventdata, handles)
% hObject    handle to Induction_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.plots_on(3)=get(hObject,'Value');
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in total.
function total_Callback(hObject, eventdata, handles)
% hObject    handle to total (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.plots_on(4)=get(hObject,'Value');
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in simulate.
function simulate_Callback(hObject, eventdata, handles)
% hObject    handle to simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.data
% global fg
% figure(fg);
clc
% x0 = -16310;
% y0 = -12810;
data = handles.data;
handles.s_num =  1;
data = pulse_simulator(handles.data);




% sen_set = open('sensor_setting.mat');
% x = sen_set.x; y = sen_set.y; z = sen_set.z;

% x = [x(1) x(2) x(4)];
% y = [y(1) y(2) y(4)];

% tt = [69424.3224608 69424.3224616 69424.3224608]; times for upward
% illuminations

% tt = [68941.64519002 68941.64519002 68941.64519002];

% for i =1:3
%     subplot(2,2,i)
%     hold all    
%     handles.data.x = x(i) - x0;
%     handles.data.y = y(i) - y0;
%     handles.data.t_shift = tt(i);
%     [t E_tot] = pulse_simulator(handles.data);
%     
%     plot(t,E_tot)
% end


% --- Executes on button press in close_pulse.
function close_pulse_Callback(hObject, eventdata, handles)
% hObject    handle to close_pulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
closeGUI(hObject,eventdata)

% --------------------------------------------------------------------
function help_me_Callback(hObject, eventdata, handles)
str=['Real help is not availble at this time. However you can go to ' ...
    'File --> Open and then open the file called "pulse_parameter_test.mat".' ...
    'This example contains the values that leads to a quick result' ...
    'Better result may be obtained tweeking the values'];

msgbox(str)



function total_t_Callback(hObject, eventdata, handles)
% hObject    handle to total_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.t2=str2double(get(hObject,'String'));
handles.data.is_saved=0;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function total_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to total_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function closeGUI(hObject,eventdata)
    handles= guidata(hObject);
    % Just don't ask questions but dump current values in a default folder
    % just before exit.
    data=handles.data;
    data.is_saved = 0;
    fn=sprintf('%spulse_para.mat',tempdir);
    save(fn, '-struct','data')
    delete(handles.figure1)
    

        
