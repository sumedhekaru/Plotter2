function varargout = ids_positions(varargin)
% IDS_POSITIONS MATLAB code for ids_positions.fig
%      IDS_POSITIONS, by itself, creates a new IDS_POSITIONS or raises the existing
%      singleton*.
%
%      H = IDS_POSITIONS returns the handle to a new IDS_POSITIONS or the handle to
%      the existing singleton*.
%
%      IDS_POSITIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IDS_POSITIONS.M with the given input arguments.
%
%      IDS_POSITIONS('Property','Value',...) creates a new IDS_POSITIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ids_positions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ids_positions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ids_positions

% Last Modified by GUIDE v2.5 07-Oct-2015 13:37:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ids_positions_OpeningFcn, ...
                   'gui_OutputFcn',  @ids_positions_OutputFcn, ...
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


% --- Executes just before ids_positions is made visible.
function ids_positions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ids_positions (see VARARGIN)

% Sumedhes data
handles.sen_set=open('sensor_setting.mat');
handles.g.sn=NaN;


% Choose default command line output for ids_positions
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ids_positions wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ids_positions_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on selection change in sn.
function sn_Callback(hObject, eventdata, handles)
% hObject    handle to sn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sn contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sn

l=get(hObject,'String');
sn=str2double(l(get(hObject,'Value')));
handles.g.sn=sn;
if ~isnan(handles.g.sn)
    set(handles.sid,'String',handles.sen_set.sen_IDs{sn})
    set(handles.tshift,'String',handles.sen_set.t_shift(sn*3))
    
    set(handles.vshift1,'String',handles.sen_set.vshift(sn*3-2))
    set(handles.vshift2,'String',handles.sen_set.vshift(sn*3-1))
    set(handles.vshift3,'String',handles.sen_set.vshift(sn*3))
    
    set(handles.vgain1,'String',handles.sen_set.gain(sn*3-2))
    set(handles.vgain2,'String',handles.sen_set.gain(sn*3-1))
    set(handles.vgain3,'String',handles.sen_set.gain(sn*3))
    
    set(handles.ch1_f,'String',handles.sen_set.chfreq(sn*3-2)/1e6)
    set(handles.ch2_f,'String',handles.sen_set.chfreq(sn*3-1)/1e6)
    set(handles.ch3_f,'String',handles.sen_set.chfreq(sn*3)/1e6)
    
    set(handles.de_trend_ch1,'Value',boolean(handles.sen_set.deTrend(sn*3-2)));
    set(handles.de_trend_ch2,'Value',boolean(handles.sen_set.deTrend(sn*3-1)));
    set(handles.de_trend_ch3,'Value',boolean(handles.sen_set.deTrend(sn*3)));
    
    % Notch filter values
    set(handles.notch_ch1,'String',sprintf('(%0.0f,%0.0f,%0.0f)',...
        handles.sen_set.notch_center(sn*3-2),...
        handles.sen_set.notch_bw(sn*3-2),...
        handles.sen_set.notch_at(sn*3-2)));
    
    set(handles.notch_ch2,'String',sprintf('(%0.0f,%0.0f,%0.0f)',...
        handles.sen_set.notch_center(sn*3-1),...
        handles.sen_set.notch_bw(sn*3-1),...
        handles.sen_set.notch_at(sn*3-1)));
    
    set(handles.notch_ch3,'String',sprintf('(%0.0f,%0.0f,%0.0f)',...
        handles.sen_set.notch_center(sn*3),...
        handles.sen_set.notch_bw(sn*3),...
        handles.sen_set.notch_at(sn*3)));
    
    set(handles.hp_filt_ch1,'String',handles.sen_set.hp_filt(sn*3-2));
    set(handles.hp_filt_ch2,'String',handles.sen_set.hp_filt(sn*3-1));
    set(handles.hp_filt_ch3,'String',handles.sen_set.hp_filt(sn*3));
              
    set(handles.intigrate_ch1,'Value',boolean(handles.sen_set.intigrate(sn*3-2)));
    set(handles.intigrate_ch2,'Value',boolean(handles.sen_set.intigrate(sn*3-1)));
    set(handles.intigrate_ch3,'Value',boolean(handles.sen_set.intigrate(sn*3)));
    
    set(handles.integration_gain_ch1,'String',handles.sen_set.integration_gain(sn*3-2));
    set(handles.integration_gain_ch2,'String',handles.sen_set.integration_gain(sn*3-1));
    set(handles.integration_gain_ch3,'String',handles.sen_set.integration_gain(sn*3));
    
    set(handles.plot_type_ch1,'Value',handles.sen_set.plot_type(sn*3-2));
    set(handles.plot_type_ch2,'Value',handles.sen_set.plot_type(sn*3-1));
    set(handles.plot_type_ch3,'Value',handles.sen_set.plot_type(sn*3));
    
    set(handles.x,'String',handles.sen_set.x(sn))
    set(handles.y,'String',handles.sen_set.y(sn))
    set(handles.z,'String',handles.sen_set.z(sn))
end
guidata(gcbo, handles);


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


% --- Executes on button press in add.
function add_Callback(hObject, eventdata, handles)
    
    sn = handles.g.sn;
    
    sen_set=handles.sen_set;
    
    if ~isnan(sn)
        sen_set.sen_IDs{sn}=get(handles.sid,'String');
        
        ts = str2double(get(handles.tshift,'String'));
        sen_set.t_shift(sn*3-2:sn*3)=[ts ts ts];
        
        sen_set.vshift(sn*3-2)=str2double(get(handles.vshift1,'String'));
        sen_set.vshift(sn*3-1)=str2double(get(handles.vshift2,'String'));
        sen_set.vshift(sn*3)=str2double(get(handles.vshift3,'String'));
        
        sen_set.gain(sn*3-2)=str2double(get(handles.vgain1,'String'));
        sen_set.gain(sn*3-1)=str2double(get(handles.vgain2,'String'));
        sen_set.gain(sn*3)=str2double(get(handles.vgain3,'String'));
                     
        sen_set.chfreq(sn*3-2)=str2double(get(handles.ch1_f,'String'))*1e6;
        sen_set.chfreq(sn*3-1)=str2double(get(handles.ch2_f,'String'))*1e6;
        sen_set.chfreq(sn*3)=str2double(get(handles.ch3_f,'String'))*1e6;
                
        sen_set.deTrend(sn*3-2)= get(handles.de_trend_ch1,'Value');
        sen_set.deTrend(sn*3-1)= get(handles.de_trend_ch2,'Value');
        sen_set.deTrend(sn*3)= get(handles.de_trend_ch3,'Value');
        
        sen_set = set_notch_filter_values(handles,sen_set,sn);
               
        sen_set.hp_filt(sn*3-2)=str2double(get(handles.hp_filt_ch1,'String'));
        sen_set.hp_filt(sn*3-1)=str2double(get(handles.hp_filt_ch2,'String'));
        sen_set.hp_filt(sn*3)=str2double(get(handles.hp_filt_ch3,'String'));
        
        sen_set.intigrate(sn*3-2)= get(handles.intigrate_ch1,'Value');
        sen_set.intigrate(sn*3-1)= get(handles.intigrate_ch2,'Value');
        sen_set.intigrate(sn*3)= get(handles.intigrate_ch3,'Value');
        
        sen_set.integration_gain(sn*3-2)= str2double(get(handles.integration_gain_ch1,'String'));
        sen_set.integration_gain(sn*3-1)= str2double(get(handles.integration_gain_ch2,'String'));
        sen_set.integration_gain(sn*3)= str2double(get(handles.integration_gain_ch3,'String'));
        
        
        sen_set.plot_type(sn*3-2)= get(handles.plot_type_ch1,'Value');
        sen_set.plot_type(sn*3-1)= get(handles.plot_type_ch2,'Value');
        sen_set.plot_type(sn*3)= get(handles.plot_type_ch3,'Value');
               
        sen_set.x(sn)=str2double(get(handles.x,'String'));
        sen_set.y(sn)=str2double(get(handles.y,'String'));
        sen_set.z(sn)=str2double(get(handles.z,'String'));
        
       
        set(findall(0,'Tag',sprintf('sen%2.2i',sn)),'String',sen_set.sen_IDs{sn})
        
              
        msg = sprintf('Sensor # %i: entry was succesfully added/modified',sn);
        
        
        set(handles.msg,'String',msg,'FontWeight','bold','ForegroundColor',[0 0 1])
        handles.sen_set=sen_set;
        guidata(gcbo, handles);
        
        % Update plotter handle structure
        hp = guidata(findall(0,'Tag','plotter2'));
        hp.sen_set = sen_set;
        guidata(findall(0,'Tag','plotter2'),hp)
        
        % Save the settings
        save('sensor_setting.mat','-Struct','sen_set')
        
        pause(2)
        try
            set(handles.msg,'String','Change values and then press "Add" to make it permanant',...
                'FontWeight','normal','ForegroundColor',[0 0 0])
        end
            
    else
        
        set(handles.msg,'String','Select sensor # first','FontWeight','bold','ForegroundColor',[1 0 0])
       
        pause(2)
        try
            set(handles.msg,'String','Change values and then press "Add" to make it permanant',...
                'FontWeight','normal','ForegroundColor',[0 0 0])
        end
        
       
    end
    



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



function ch1_f_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch1_f as text
%        str2double(get(hObject,'String')) returns contents of ch1_f as a double


% --- Executes during object creation, after setting all properties.
function ch1_f_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch2_f_Callback(hObject, eventdata, handles)
% hObject    handle to ch2_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch2_f as text
%        str2double(get(hObject,'String')) returns contents of ch2_f as a double


% --- Executes during object creation, after setting all properties.
function ch2_f_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch2_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch3_f_Callback(hObject, eventdata, handles)
% hObject    handle to ch3_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch3_f as text
%        str2double(get(hObject,'String')) returns contents of ch3_f as a double


% --- Executes during object creation, after setting all properties.
function ch3_f_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch3_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in de_trend_ch1.
function de_trend_ch1_Callback(hObject, eventdata, handles)
% hObject    handle to de_trend_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of de_trend_ch1


% --- Executes on button press in de_trend_ch2.
function de_trend_ch2_Callback(hObject, eventdata, handles)
% hObject    handle to de_trend_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of de_trend_ch2


% --- Executes on button press in de_trend_ch3.
function de_trend_ch3_Callback(hObject, eventdata, handles)
% hObject    handle to de_trend_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of de_trend_ch3



function hp_filt_ch1_Callback(hObject, eventdata, handles)
% hObject    handle to hp_filt_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hp_filt_ch1 as text
%        str2double(get(hObject,'String')) returns contents of hp_filt_ch1 as a double


% --- Executes during object creation, after setting all properties.
function hp_filt_ch1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hp_filt_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hp_filt_ch2_Callback(hObject, eventdata, handles)
% hObject    handle to hp_filt_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hp_filt_ch2 as text
%        str2double(get(hObject,'String')) returns contents of hp_filt_ch2 as a double


% --- Executes during object creation, after setting all properties.
function hp_filt_ch2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hp_filt_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hp_filt_ch3_Callback(hObject, eventdata, handles)
% hObject    handle to hp_filt_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hp_filt_ch3 as text
%        str2double(get(hObject,'String')) returns contents of hp_filt_ch3 as a double


% --- Executes during object creation, after setting all properties.
function hp_filt_ch3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hp_filt_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in intigrate_ch1.
function intigrate_ch1_Callback(hObject, eventdata, handles)
% hObject    handle to intigrate_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of intigrate_ch1
if get(hObject,'Value')
    set(handles.plot_type_ch1,'Value',2)
else
    set(handles.plot_type_ch1,'Value',1)
end


% --- Executes on button press in intigrate_ch2.
function intigrate_ch2_Callback(hObject, eventdata, handles)
% hObject    handle to intigrate_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of intigrate_ch2

if get(hObject,'Value')
    set(handles.plot_type_ch2,'Value',2)
else
    set(handles.plot_type_ch2,'Value',1)
end

% --- Executes on button press in intigrate_ch3.
function intigrate_ch3_Callback(hObject, eventdata, handles)
% hObject    handle to intigrate_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of intigrate_ch3
if get(hObject,'Value')
    set(handles.plot_type_ch3,'Value',2)
else
    set(handles.plot_type_ch3,'Value',1)
end


% --- Executes on selection change in plot_type_ch1.
function plot_type_ch1_Callback(hObject, eventdata, handles)
% hObject    handle to plot_type_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_type_ch1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_type_ch1


% --- Executes during object creation, after setting all properties.
function plot_type_ch1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_type_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plot_type_ch2.
function plot_type_ch2_Callback(hObject, eventdata, handles)
% hObject    handle to plot_type_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_type_ch2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_type_ch2


% --- Executes during object creation, after setting all properties.
function plot_type_ch2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_type_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plot_type_ch3.
function plot_type_ch3_Callback(hObject, eventdata, handles)
% hObject    handle to plot_type_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_type_ch3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_type_ch3


% --- Executes during object creation, after setting all properties.
function plot_type_ch3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_type_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function integration_gain_ch1_Callback(hObject, eventdata, handles)
% hObject    handle to integration_gain_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of integration_gain_ch1 as text
%        str2double(get(hObject,'String')) returns contents of integration_gain_ch1 as a double


% --- Executes during object creation, after setting all properties.
function integration_gain_ch1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to integration_gain_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function integration_gain_ch2_Callback(hObject, eventdata, handles)
% hObject    handle to integration_gain_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of integration_gain_ch2 as text
%        str2double(get(hObject,'String')) returns contents of integration_gain_ch2 as a double


% --- Executes during object creation, after setting all properties.
function integration_gain_ch2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to integration_gain_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function integration_gain_ch3_Callback(hObject, eventdata, handles)
% hObject    handle to integration_gain_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of integration_gain_ch3 as text
%        str2double(get(hObject,'String')) returns contents of integration_gain_ch3 as a double


% --- Executes during object creation, after setting all properties.
function integration_gain_ch3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to integration_gain_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function notch_ch1_Callback(hObject, eventdata, handles)
% hObject    handle to notch_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of notch_ch1 as text
%        str2double(get(hObject,'String')) returns contents of notch_ch1 as a double


% --- Executes during object creation, after setting all properties.
function notch_ch1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to notch_ch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function notch_ch2_Callback(hObject, eventdata, handles)
% hObject    handle to notch_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of notch_ch2 as text
%        str2double(get(hObject,'String')) returns contents of notch_ch2 as a double


% --- Executes during object creation, after setting all properties.
function notch_ch2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to notch_ch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function notch_ch3_Callback(hObject, eventdata, handles)
% hObject    handle to notch_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of notch_ch3 as text
%        str2double(get(hObject,'String')) returns contents of notch_ch3 as a double

% --- Executes during object creation, after setting all properties.
function notch_ch3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to notch_ch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sen_set = set_notch_filter_values(handles,sen_set,sn)

try
    fil_str = get(handles.notch_ch3,'String');

    inds = strfind(fil_str,',');
    cf = str2double(fil_str(2:inds(1)-1));
    bw = str2double(fil_str(inds(1)+1:inds(2)-1));
    at = str2double(fil_str(inds(2)+1:end-1));

    sen_set.notch_center(sn*3)= cf;
    sen_set.notch_bw(sn*3) = bw;
    sen_set.notch_at(sn*3) = at;

    
catch
    set(handles.notch_ch3,'String','(0,100,-25)');
    sen_set.notch_center(sn*3)= 0;
    sen_set.notch_bw(sn*3) = 100;
    sen_set.notch_at(sn*3) = -25;    
end


try
    fil_str = get(handles.notch_ch2,'String');
    
    inds = strfind(fil_str,',');
    cf = str2double(fil_str(2:inds(1)-1));
    bw = str2double(fil_str(inds(1)+1:inds(2)-1));
    at = str2double(fil_str(inds(2)+1:end-1));
    
    sen_set.notch_center(sn*3-1)= cf;
    sen_set.notch_bw(sn*3-1) = bw;
    sen_set.notch_at(sn*3-1) = at;
    
catch
    set(handles.notch_ch2,'String','(0,100,-25)');
    sen_set.notch_center(sn*3-1)= 0;
    sen_set.notch_bw(sn*3-1) = 100;
    sen_set.notch_at(sn*3-1) = -25;
end


try
    fil_str = get(handles.notch_ch1,'String');
    
    inds = strfind(fil_str,',');
    cf = str2double(fil_str(2:inds(1)-1));
    bw = str2double(fil_str(inds(1)+1:inds(2)-1));
    at = str2double(fil_str(inds(2)+1:end-1));
    
    sen_set.notch_center(sn*3-2)= cf;
    sen_set.notch_bw(sn*3-2) = bw;
    sen_set.notch_at(sn*3-2) = at;
    
catch
    set(handles.notch_ch1,'String','(0,100,-25)');
    sen_set.notch_center(sn*3-2)= 0;
    sen_set.notch_bw(sn*3-2) = 100;
    sen_set.notch_at(sn*3-2) = -25;
end
    
