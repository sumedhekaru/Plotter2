function varargout = tgf_plotter(varargin)
% TGF_PLOTTER MATLAB code for tgf_plotter.fig
%      TGF_PLOTTER, by itself, creates a new TGF_PLOTTER or raises the existing
%      singleton*.
%
%      H = TGF_PLOTTER returns the handle to a new TGF_PLOTTER or the handle to
%      the existing singleton*.
%
%      TGF_PLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TGF_PLOTTER.M with the given input arguments.
%
%      TGF_PLOTTER('Property','Value',...) creates a new TGF_PLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tgf_plotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tgf_plotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tgf_plotter

% Last Modified by GUIDE v2.5 08-Oct-2012 22:58:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tgf_plotter_OpeningFcn, ...
                   'gui_OutputFcn',  @tgf_plotter_OutputFcn, ...
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


% --- Executes just before tgf_plotter is made visible.
function tgf_plotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tgf_plotter (see VARARGIN)

try
    b = open([tempdir 'tgf_plotter_last_gui_data.mat']);
    handles.b = b;
catch
    handles.b.fn = '';
    handles.b.dn = '';
    handles.b.radialOn = 1;
    handles.b.azimuOn = 1;
    handles.b.B_r_tshift = 0;
    handles.b.B_r_bshift = 0;
    handles.b.B_r_gain = 1;
    handles.b.B_phi_tshift = 0;
    handles.b.B_phi_bshift = 0;
    handles.b.B_phi_gain = 1;
    handles.b.new_figure = 0;
    handles.b.save_fn = '';
    handles.b.save_dn = '';
end

handles.temp.B_r_plot = NaN;
handles.temp.B_phi_plot = NaN;

set_values(handles)

% Choose default command line output for tgf_plotter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tgf_plotter wait for user response (see UIRESUME)
% uiwait(handles.tgf_plotter);


% --- Outputs from this function are returned to the command line.
function varargout = tgf_plotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function file_name_Callback(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of file_name as text
%        str2double(get(hObject,'String')) returns contents of file_name as a double


% --- Executes during object creation, after setting all properties.
function file_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in choose_file.
function choose_file_Callback(hObject, eventdata, handles)
b = handles.b;

[FileName,PathName] = uigetfile([b.dn '\*.mat']);

if ~isequal(FileName,0)
    b.dn = PathName;
    b.fn = FileName;
    set(handles.file_name,'String',FileName)
    handles.b = b;
    guidata(hObject, handles);
end


function B_r_gain_Callback(hObject, eventdata, handles)

handles.b.B_r_gain = str2double(get(hObject,'String'));
if ~isnan(handles.temp.B_r_plot)
    handles = update_B_r(handles);
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function B_r_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_r_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B_r_bshift_Callback(hObject, eventdata, handles)
handles.b.B_r_bshift = str2double(get(hObject,'String'));
if ~isnan(handles.temp.B_r_plot)
    handles = update_B_r(handles);
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function B_r_bshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_r_bshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B_r_tshift_Callback(hObject, eventdata, handles)
handles.b.B_r_tshift = str2double(get(hObject,'String'));
if ~isnan(handles.temp.B_r_plot)
    handles = update_B_r(handles);
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function B_r_tshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_r_tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B_phi_tshift_Callback(hObject, eventdata, handles)
handles.b.B_phi_tshift = str2double(get(hObject,'String'));
if ~isnan(handles.temp.B_phi_plot)
    handles = update_B_phi(handles);
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function B_phi_tshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_phi_tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B_phi_bshift_Callback(hObject, eventdata, handles)
handles.b.B_phi_bshift = str2double(get(hObject,'String'));

if ~isnan(handles.temp.B_phi_plot)
   handles = update_B_phi(handles);
end

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function B_phi_bshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_phi_bshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B_phi_gain_Callback(hObject, eventdata, handles)
handles.b.B_phi_gain = str2double(get(hObject,'String'));
if ~isnan(handles.temp.B_phi_plot)
    handles = update_B_phi(handles);
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function B_phi_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_phi_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radial.
function radial_Callback(hObject, eventdata, handles)
handles.b.radialOn = get(hObject,'Value');
guidata(hObject, handles);



% --- Executes on button press in azimuthal.
function azimuthal_Callback(hObject, eventdata, handles)

handles.b.azimuOn = get(hObject,'Value');
guidata(hObject, handles);


% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function open_Callback(hObject, eventdata, handles)
b = handles.b;
[fn, pn] = uigetfile([b.save_dn '\*.tgf']);

b = load([pn fn],'-mat');
handles.b = b;
set_values(handles)
guidata(hObject, handles);

% --------------------------------------------------------------------
function save_as_Callback(hObject, eventdata, handles)
b = handles.b;

if ~isempty(b.save_dn)
    working_dir = b.save_dn;
elseif ~isempty(b.dn)
    working_dir = b.dn;
else
    working_dir = '';
end
    
try
    h=guidata(findall(0,'Tag','plotter2'));
    g = h.g;
    file=sprintf('%s//%s%s%s_%5.5us_%3.3ums_.tgf',working_dir,...
        g.YYYY{:},g.MM{:},g.DD{:},floor(g.t1),floor((g.t1-floor(g.t1))*1000));
catch
    file = working_dir;
end

[fn,pn] = uiputfile('*.tgf','Save As..',file);

if fn(1)~=0
    b.save_fn=fn;
    b.save_dn=pn;
    fn=[pn fn];
    save(fn,'-struct', 'b')
    handles.b=b;
    guidata(hObject, handles);
end

function set_values(handles)

b = handles.b;

set(handles.radial,'Value',b.radialOn)
set(handles.azimuthal,'Value',b.azimuOn)
set(handles.file_name,'String',b.fn)
set(handles.B_r_tshift,'String',b.B_r_tshift)
set(handles.B_r_bshift,'String',b.B_r_bshift)
set(handles.B_r_gain,'String',b.B_r_gain)
set(handles.B_phi_tshift,'String',b.B_phi_tshift)
set(handles.B_phi_bshift,'String',b.B_phi_bshift)
set(handles.B_phi_gain,'String',b.B_phi_gain)
set(handles.new_figure,'Value',b.new_figure)


% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)

b = handles.b;
% Get data
TGFd = TGF_extractor([b.dn b.fn]);
handles.TGFd = TGFd;

if b.new_figure
    figure
    lg = {};
    box on
    xlabel('Time (s)')
    ylabel('Magnetic Field (nT)')
    tools2fig
else
    figure(gcf)
    lg = get(legend(gca),'String');
end
hold all

if b.radialOn
    x = TGFd.t + b.B_r_tshift;
    y = TGFd.B_r*b.B_r_gain + b.B_r_bshift;
    handles.temp.B_r_plot = plot(x,y);
    
    if b.B_r_gain == 1 && b.B_r_bshift == 0
        lgn = ['B_r ' TGFd.station '-' TGFd.sensor];
    elseif b.B_r_gain == 1
        lgn =  ['B_r ' TGFd.station '-' TGFd.sensor sprintf(' (+%0.2f)', b.B_r_bshift)];
    elseif b.B_r_bshift == 0
        lgn = ['B_r ' TGFd.station '-' TGFd.sensor sprintf(' (x %0.2f)',b.B_r_gain)];
    else
        lgn = ['B_r ' TGFd.station '-' TGFd.sensor sprintf(' (x %0.2f+%0.2f)',b.B_r_gain,b.B_r_bshift)];
    end
else
    lgn = [];
end


if b.azimuOn
    x = TGFd.t + b.B_phi_tshift;
    y = TGFd.B_phi*b.B_phi_gain + b.B_phi_bshift;
    handles.temp.B_phi_plot = plot(x,y);
    
    if b.B_phi_gain == 1 && b.B_phi_bshift == 0
        lgn2 = ['B_\phi ' TGFd.station '-' TGFd.sensor];
    elseif b.B_phi_gain == 1
        lgn2 = ['B_\phi ' TGFd.station '-' TGFd.sensor sprintf(' (+%0.2f)', b.B_phi_bshift)];
    elseif b.B_r_bshift == 0
        lgn2 = ['B_\phi ' TGFd.station '-' TGFd.sensor sprintf(' (x %0.2f)',b.B_phi_gain)];
    else
        lgn2 = ['B_\phi ' TGFd.station '-' TGFd.sensor sprintf(' (x %0.2f+%0.2f)',b.B_phi_gain,b.B_phi_bshift)];
    end
else
    lgn2 = [];
end


n1 = find(strcmp(lg, 'LINET'));
n2 = find(strcmp(lg, 'PBFA'));
n3 = find(strcmp(lg,'LDAR2'));
n4 = find(strcmp(lg,'CGLSS'));

n = min([n1 n2 n3 n4]);

lgt1 = lg(1:n-1);
lgt2 = lg(n:end);

handles.temp.B_phi_lg = lgn2;
handles.temp.B_r_lg = lgn;


legend([lgt1 lgn lgn2 lgt2])
guidata(hObject, handles);


% --- Executes on button press in new_figure.
function new_figure_Callback(hObject, eventdata, handles)
handles.b.new_figure = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes when user attempts to close tgf_plotter.
function tgf_plotter_CloseRequestFcn(hObject, eventdata, handles)

b = handles.b;

try
   save([tempdir 'tgf_plotter_last_gui_data.mat'],'-struct', 'b')
end

delete(hObject)

function handles = update_B_phi(handles)

    lgnd = get(legend(gca),'String');
    N = find(strcmp(lgnd,handles.temp.B_phi_lg));
    
    b = handles.b;
    TGFd = handles.TGFd;
    color = get(handles.temp.B_phi_plot,'Color');
    
    if b.B_phi_gain == 1 && b.B_phi_bshift == 0
        lg = ['B_\phi ' TGFd.station '-' TGFd.sensor];
    elseif b.B_phi_gain == 1
        lg = ['B_\phi ' TGFd.station '-' TGFd.sensor sprintf(' (+%0.2f)', b.B_phi_bshift)];
    elseif b.B_r_bshift == 0
        lg  = ['B_\phi ' TGFd.station '-' TGFd.sensor sprintf(' (x %0.2f)',b.B_phi_gain)];
    else
        lg = ['B_\phi ' TGFd.station '-' TGFd.sensor sprintf(' (x %0.2f+%0.2f)',b.B_phi_gain,b.B_phi_bshift)];
    end
    
    %legend(handles.temp.B_phi_plot,lg)
      lgnd{N} = lg;
    legend(lgnd)
    
    
    x = TGFd.t + b.B_phi_tshift;
    y = TGFd.B_phi*b.B_phi_gain + b.B_phi_bshift;
    
    h = plot(x,y,'Color',color,'DisplayName',lg);
    delete(handles.temp.B_phi_plot)
    handles.temp.B_phi_plot = h;
    handles.temp.B_phi_lg = lg;
    
    
function handles = update_B_r(handles)

     
    lgnd = get(legend(gca),'String');
    N = find(strcmp(lgnd,handles.temp.B_r_lg));
   
 
    b = handles.b;
    TGFd = handles.TGFd;
    color = get(handles.temp.B_r_plot,'Color');
    
    if b.B_r_gain == 1 && b.B_r_bshift == 0
        lg = ['B_r ' TGFd.station '-' TGFd.sensor];
    elseif b.B_phi_gain == 1
        lg = ['B_r ' TGFd.station '-' TGFd.sensor sprintf(' (+%0.2f)', b.B_r_bshift)];
    elseif b.B_r_bshift == 0
        lg  = ['B_r ' TGFd.station '-' TGFd.sensor sprintf(' (x %0.2f)',b.B_r_gain)];
    else
        lg = ['B_r ' TGFd.station '-' TGFd.sensor sprintf(' (x %0.2f+%0.2f)',b.B_r_gain,b.B_r_bshift)];
    end
   
    %legend(handles.temp.B_phi_plot,lg)
    lgnd{N} = lg;
    legend(lgnd)
    
    x = TGFd.t + b.B_r_tshift;
    y = TGFd.B_r*b.B_r_gain + b.B_r_bshift;
    
    h = plot(x,y,'Color',color);
    delete(handles.temp.B_r_plot)
    handles.temp.B_r_plot = h;
    handles.temp.B_r_lg = lg;
