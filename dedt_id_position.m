function varargout = dedt_id_position(varargin)
% DEDT_ID_POSITION MATLAB code for dedt_id_position.fig
%      DEDT_ID_POSITION, by itself, creates a new DEDT_ID_POSITION or raises the existing
%      singleton*.
%
%      H = DEDT_ID_POSITION returns the handle to a new DEDT_ID_POSITION or the handle to
%      the existing singleton*.
%
%      DEDT_ID_POSITION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEDT_ID_POSITION.M with the given input arguments.
%
%      DEDT_ID_POSITION('Property','Value',...) creates a new DEDT_ID_POSITION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dedt_id_position_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dedt_id_position_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dedt_id_position

% Last Modified by GUIDE v2.5 08-Sep-2011 19:00:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dedt_id_position_OpeningFcn, ...
                   'gui_OutputFcn',  @dedt_id_position_OutputFcn, ...
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


% --- Executes just before dedt_id_position is made visible.
function dedt_id_position_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dedt_id_position (see VARARGIN)

% Sumedhes data
handles.sen_set=open('sensor_setting.mat');
handles.g.sn=NaN;

%%% Show old values to user'
dEdt=handles.sen_set.dEdt;

set(handles.sid,'String',dEdt.sn)
set(handles.tshift,'String',dEdt.tshift(1))

set(handles.vshift1,'String',dEdt.vshift(1))
set(handles.vshift2,'String',dEdt.vshift(2))
set(handles.vshift3,'String',dEdt.vshift(3))
set(handles.vshift4,'String',dEdt.vshift(4))

set(handles.vgain1,'String',dEdt.gain(1))
set(handles.vgain2,'String',dEdt.gain(2))
set(handles.vgain3,'String',dEdt.gain(3))
set(handles.vgain4,'String',dEdt.gain(4))

set(handles.x,'String',dEdt.x)
set(handles.y,'String',dEdt.y)
set(handles.z,'String',dEdt.z)






% Choose default command line output for dedt_id_position
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dedt_id_position wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dedt_id_position_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function tshift_Callback(hObject, eventdata, handles)
% hObject    handle to tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tshift as text
%        str2double(get(hObject,'String')) returns contents of tshift as a double


% --- Executes during object creation, after setting all properties.
function tshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vshift1_Callback(hObject, eventdata, handles)
% hObject    handle to vshift1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vshift1 as text
%        str2double(get(hObject,'String')) returns contents of vshift1 as a double


% --- Executes during object creation, after setting all properties.
function vshift1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vshift1 (see GCBO)
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

% Hints: get(hObject,'String') returns contents of x as text
%        str2double(get(hObject,'String')) returns contents of x as a double


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



function sid_Callback(hObject, eventdata, handles)
% hObject    handle to sid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sid as text
%        str2double(get(hObject,'String')) returns contents of sid as a double


% --- Executes during object creation, after setting all properties.
function sid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sid (see GCBO)
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

% Hints: get(hObject,'String') returns contents of y as text
%        str2double(get(hObject,'String')) returns contents of y as a double


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



function z_Callback(hObject, eventdata, handles)
% hObject    handle to z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z as text
%        str2double(get(hObject,'String')) returns contents of z as a double


% --- Executes during object creation, after setting all properties.
function z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on selection change in sn.
% function sn_Callback(hObject, eventdata, handles)
% % hObject    handle to sn (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = cellstr(get(hObject,'String')) returns sn contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from sn
% l=get(hObject,'String');
% sn=str2double(l(get(hObject,'Value')));
% handles.g.sn=sn;
% if ~isnan(handles.g.sn)
%     set(handles.sid,'String',handles.sen_set.sen_IDs{sn})
%     set(handles.tshift,'String',handles.sen_set.t_shift(sn*3))
%     
%     set(handles.vshift1,'String',handles.sen_set.vshift(sn*3-2))
%     set(handles.vshift2,'String',handles.sen_set.vshift(sn*3-1))
%     set(handles.vshift3,'String',handles.sen_set.vshift(sn*3))
%     
%     set(handles.vgain1,'String',handles.sen_set.gain(sn*3-1))
%     set(handles.vgain2,'String',handles.sen_set.gain(sn*3-2))
%     set(handles.vgain3,'String',handles.sen_set.gain(sn*3))
%     
%     set(handles.x,'String',handles.sen_set.x(sn))
%     set(handles.y,'String',handles.sen_set.y(sn))
%     set(handles.z,'String',handles.sen_set.z(sn))
% end
% guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function sn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
    
sen_set=handles.sen_set;


sen_set.dEdt.sn=get(handles.sid,'String');

ts = str2double(get(handles.tshift,'String'));
sen_set.dEdt.tshift=[ts ts ts ts];

sen_set.dEdt.vshift(1)=str2double(get(handles.vshift1,'String'));
sen_set.dEdt.vshift(2)=str2double(get(handles.vshift2,'String'));
sen_set.dEdt.vshift(3)=str2double(get(handles.vshift3,'String'));
sen_set.dEdt.vshift(4)=str2double(get(handles.vshift4,'String'));

sen_set.dEdt.gain(1)=str2double(get(handles.vgain1,'String'));
sen_set.dEdt.gain(2)=str2double(get(handles.vgain2,'String'));
sen_set.dEdt.gain(3)=str2double(get(handles.vgain3,'String'));
sen_set.dEdt.gain(4)=str2double(get(handles.vgain4,'String'));

sen_set.dEdt.x=str2double(get(handles.x,'String'));
sen_set.dEdt.y=str2double(get(handles.y,'String'));
sen_set.dEdt.z=str2double(get(handles.z,'String'));

set(findall(0,'Tag','dEdt_sn'),'String',sen_set.dEdt.sn)

save('sensor_setting.mat','-struct', 'sen_set')

delete(gcf)

       
    
   
    



% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(gcf)


function vgain1_Callback(hObject, eventdata, handles)
% hObject    handle to vgain1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vgain1 as text
%        str2double(get(hObject,'String')) returns contents of vgain1 as a double


% --- Executes during object creation, after setting all properties.
function vgain1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vgain1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vshift2_Callback(hObject, eventdata, handles)
% hObject    handle to vshift2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vshift2 as text
%        str2double(get(hObject,'String')) returns contents of vshift2 as a double


% --- Executes during object creation, after setting all properties.
function vshift2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vshift2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vshift3_Callback(hObject, eventdata, handles)
% hObject    handle to vshift3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vshift3 as text
%        str2double(get(hObject,'String')) returns contents of vshift3 as a double


% --- Executes during object creation, after setting all properties.
function vshift3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vshift3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vgain2_Callback(hObject, eventdata, handles)
% hObject    handle to vgain2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vgain2 as text
%        str2double(get(hObject,'String')) returns contents of vgain2 as a double


% --- Executes during object creation, after setting all properties.
function vgain2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vgain2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vgain3_Callback(hObject, eventdata, handles)
% hObject    handle to vgain3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vgain3 as text
%        str2double(get(hObject,'String')) returns contents of vgain3 as a double


% --- Executes during object creation, after setting all properties.
function vgain3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vgain3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vshift4_Callback(hObject, eventdata, handles)
% hObject    handle to vshift4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vshift4 as text
%        str2double(get(hObject,'String')) returns contents of vshift4 as a double


% --- Executes during object creation, after setting all properties.
function vshift4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vshift4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vgain4_Callback(hObject, eventdata, handles)
% hObject    handle to vgain4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vgain4 as text
%        str2double(get(hObject,'String')) returns contents of vgain4 as a double


% --- Executes during object creation, after setting all properties.
function vgain4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vgain4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
