function varargout = field_mill_plotter(varargin)
% FIELD_MILL_PLOTTER MATLAB code for field_mill_plotter.fig
%      FIELD_MILL_PLOTTER, by itself, creates a new FIELD_MILL_PLOTTER or raises the existing
%      singleton*.
%
%      H = FIELD_MILL_PLOTTER returns the handle to a new FIELD_MILL_PLOTTER or the handle to
%      the existing singleton*.
%
%      FIELD_MILL_PLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIELD_MILL_PLOTTER.M with the given input arguments.
%
%      FIELD_MILL_PLOTTER('Property','Value',...) creates a new FIELD_MILL_PLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before field_mill_plotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to field_mill_plotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help field_mill_plotter

% Last Modified by GUIDE v2.5 25-Nov-2013 11:21:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @field_mill_plotter_OpeningFcn, ...
                   'gui_OutputFcn',  @field_mill_plotter_OutputFcn, ...
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


% --- Executes just before field_mill_plotter is made visible.
function field_mill_plotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to field_mill_plotter (see VARARGIN)

h = varargin{1};
handles.h = h;

p1 = get( h.plotter2,'Position');
set(handles.field_mill_plotter,'Position',[p1(1)+82 p1(2)+25 43.8 16.9231])


% Set values
set_values(handles);

% Choose default command line output for field_mill_plotter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes field_mill_plotter wait for user response (see UIRESUME)
% uiwait(handles.field_mill_plotter);




% --- Outputs from this function are returned to the command line.
function varargout = field_mill_plotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in KSC01.
function KSC01_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2'));    
    h.g.KSC_FM(1) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC02.
function KSC02_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2'));    
    h.g.KSC_FM(2) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC03.
function KSC03_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); ;
    h.g.KSC_FM(3) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC04.
function KSC04_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(4) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC05.
function KSC05_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); ;
    h.g.KSC_FM(5) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC06.
function KSC06_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(6) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC07.
function KSC07_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(7) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC08.
function KSC08_Callback(hObject, eventdata, handles)
    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(8) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC09.
function KSC09_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(9) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC10.
function KSC10_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(10) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC11.
function KSC11_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(11) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC12.
function KSC12_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(12) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)
    


% --- Executes on button press in KSC13.
function KSC13_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(13) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC14.
function KSC14_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(14) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC15.
function KSC15_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(15) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC16.
function KSC16_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(16) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC17.
function KSC17_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(17) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC18.
function KSC18_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(18) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC19.
function KSC19_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(19) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC20.
function KSC20_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(20) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC21.
function KSC21_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(21) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC22.
function KSC22_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(22) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC23.
function KSC23_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(23) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC24.
function KSC24_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(24) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC25.
function KSC25_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(25) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC26.
function KSC26_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(26) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC27.
function KSC27_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(27) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC28.
function KSC28_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(28) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC29.
function KSC29_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(29) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC30.
function KSC30_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(30) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)

% --- Executes on button press in KSC31.
function KSC31_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(31) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC32.
function KSC32_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(32) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC33.
function KSC33_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(33) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in KSC34.
function KSC34_Callback(hObject, eventdata, handles)

    h=guidata(findall(0,'Tag','plotter2')); 
    h.g.KSC_FM(34) = get(hObject,'Value');
    handles.h = h;
    guidata(hObject, handles);
    guidata(findall(0,'Tag','plotter2'), h)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function set_values(handles)

set(handles.KSC01,'Value',logical(handles.h.g.KSC_FM(1)));
set(handles.KSC02,'Value',logical(handles.h.g.KSC_FM(2)));
set(handles.KSC03,'Value',logical(handles.h.g.KSC_FM(3)));
set(handles.KSC04,'Value',logical(handles.h.g.KSC_FM(4)));
set(handles.KSC05,'Value',logical(handles.h.g.KSC_FM(5)));
set(handles.KSC06,'Value',logical(handles.h.g.KSC_FM(6)));
set(handles.KSC07,'Value',logical(handles.h.g.KSC_FM(7)));
set(handles.KSC08,'Value',logical(handles.h.g.KSC_FM(8)));
set(handles.KSC09,'Value',logical(handles.h.g.KSC_FM(9)));
set(handles.KSC10,'Value',logical(handles.h.g.KSC_FM(10)));
set(handles.KSC11,'Value',logical(handles.h.g.KSC_FM(11)));
set(handles.KSC12,'Value',logical(handles.h.g.KSC_FM(12)));
set(handles.KSC13,'Value',logical(handles.h.g.KSC_FM(13)));
set(handles.KSC14,'Value',logical(handles.h.g.KSC_FM(14)));
set(handles.KSC15,'Value',logical(handles.h.g.KSC_FM(15)));
set(handles.KSC16,'Value',logical(handles.h.g.KSC_FM(16)));
set(handles.KSC17,'Value',logical(handles.h.g.KSC_FM(17)));
set(handles.KSC18,'Value',logical(handles.h.g.KSC_FM(18)));
set(handles.KSC19,'Value',logical(handles.h.g.KSC_FM(19)));
set(handles.KSC20,'Value',logical(handles.h.g.KSC_FM(20)));
set(handles.KSC21,'Value',logical(handles.h.g.KSC_FM(21)));
set(handles.KSC22,'Value',logical(handles.h.g.KSC_FM(22)));
set(handles.KSC23,'Value',logical(handles.h.g.KSC_FM(23)));
set(handles.KSC24,'Value',logical(handles.h.g.KSC_FM(24)));
set(handles.KSC25,'Value',logical(handles.h.g.KSC_FM(25)));
set(handles.KSC26,'Value',logical(handles.h.g.KSC_FM(26)));
set(handles.KSC27,'Value',logical(handles.h.g.KSC_FM(27)));
set(handles.KSC28,'Value',logical(handles.h.g.KSC_FM(28)));
set(handles.KSC29,'Value',logical(handles.h.g.KSC_FM(29)));
set(handles.KSC30,'Value',logical(handles.h.g.KSC_FM(30)));
set(handles.KSC31,'Value',logical(handles.h.g.KSC_FM(31)));
set(handles.KSC32,'Value',logical(handles.h.g.KSC_FM(32)));
set(handles.KSC33,'Value',logical(handles.h.g.KSC_FM(33)));
set(handles.KSC34,'Value',logical(handles.h.g.KSC_FM(34)));
