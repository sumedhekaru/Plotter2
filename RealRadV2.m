function varargout = RealRadV2(varargin)
% REALRADV2 MATLAB code for RealRadV2.fig
%      REALRADV2, by itself, creates a new REALRADV2 or raises the existing
%      singleton*.
%
%      H = REALRADV2 returns the handle to a new REALRADV2 or the handle to
%      the existing singleton*.
%
%      REALRADV2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REALRADV2.M with the given input arguments.
%
%      REALRADV2('Property','Value',...) creates a new REALRADV2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RealRadV2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RealRadV2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RealRadV2

% Last Modified by GUIDE v2.5 07-Aug-2014 23:56:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RealRadV2_OpeningFcn, ...
                   'gui_OutputFcn',  @RealRadV2_OutputFcn, ...
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


% --- Executes just before RealRadV2 is made visible.
function RealRadV2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RealRadV2 (see VARARGIN)

% Plot Mississippi map
axes(handles.axes1) 
hold all
data = load('ms_map.mat');
hLine = plot(data.x/1000,data.y/1000,'Color',[0.6 0.6 0.6]);
set(get(get(hLine,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off')
ylim([-500 200])
xlim([-350 350])
daspect([1 1 1])

% Start polling data for the first time
handles.timer = timer(...
    'ExecutionMode', 'fixedRate', ...   % Run timer repeatedly
    'Period', 55, ...                % Initial period is 1 sec.
    'TimerFcn', {@update_display,handles},...
    'StartDelay',5); % Specify callback

start(handles.timer);

% Choose default command line output for RealRadV2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RealRadV2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RealRadV2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in VerticalPlane.
function VerticalPlane_Callback(hObject, eventdata, handles)
% hObject    handle to VerticalPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% Necessary to provide this function to prevent timer callback
% from causing an error after GUI code stops executing.
% Before exiting, if the timer is running, stop it.

try
    if strcmp(get(handles.timer, 'Running'), 'on')
        stop(handles.timer);
    end
    % Destroy timer
    delete(handles.timer)
end
% Destroy figure
delete(hObject);

function update_display(handles)


stationId = 'KNQA';


% Base radar downloading site (obtained from GRAnalysit)
base_site = 'http://mesonet-nexrad.agron.iastate.edu/level2/raw/';

% Base folder for saving data
rbf = 'C:\Users\sumedhe\Desktop\RealTimeRadar\';


%% Start the program
str = urlread([base_site stationId]);

indx = strfind(str,stationId);

% Just before last file (last file might be still writing).
lastFile = str(indx(end-2):indx(end-2)+17);

year = lastFile(6:9);
month = lastFile(10:11);
date = lastFile(12:13);
hour = lastFile(15:16);
mm = lastFile(17:18);

% Download this file

bf1 = sprintf('%s/%s/Compressed/%s/%s/%s/%s/',rbf,stationId,year,month,date,hour);

if ~exist(bf1,'dir')
    mkdir(bf1)
end

% If not already downloaded, let's download it
if ~exist([bf1 lastFile],'file')
    urlwrite([base_site stationId '/' lastFile],[bf1 lastFile]);
end

bf2 = sprintf('%s/%s/%s/%s/%s/%s/',rbf,stationId,year,month,date,hour);

if ~exist(bf2,'dir')
    mkdir(bf2)
end


% If have not Extracted it, let's do it
if ~exist([bf2 lastFile '.netcdf'],'file')
    cmd = sprintf('"C:/Program Files/Java/jre6/bin/java" -classpath "C:/Users/sumedhe/Desktop/Plotter2/toolsUI-4.3.jar" ucar.nc2.FileWriter -in "%s%s" -out "%s%s.netcdf"',...
        bf1,lastFile,bf2,lastFile);
    system(cmd)
end


% All done, let's plot it
sen_set = open('sensor_setting.mat');
sen_set.radarFn = [bf2 lastFile '.netcdf'];
sen_set.radEleAngInd = 1;
axis(handles.axis1)
radar_plot(sen_set)

