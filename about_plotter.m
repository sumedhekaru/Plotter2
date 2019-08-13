function varargout = about_plotter(varargin)
% ABOUT_PLOTTER M-file for about_plotter.fig
%      ABOUT_PLOTTER, by itself, creates a new ABOUT_PLOTTER or raises the existing
%      singleton*.
%
%      H = ABOUT_PLOTTER returns the handle to a new ABOUT_PLOTTER or the handle to
%      the existing singleton*.
%
%      ABOUT_PLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ABOUT_PLOTTER.M with the given input arguments.
%
%      ABOUT_PLOTTER('Property','Value',...) creates a new ABOUT_PLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before about_plotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to about_plotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help about_plotter

% Last Modified by GUIDE v2.5 13-May-2010 14:59:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @about_plotter_OpeningFcn, ...
                   'gui_OutputFcn',  @about_plotter_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = varargin{1};
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before about_plotter is made visible.
function about_plotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to about_plotter (see VARARGIN)

% Choose default command line output for about_plotter
handles.output = hObject;
ver=strcat('Version :',varargin(1));
date=strcat('Date : ',varargin(2));
set(handles.ver,'String',ver)
set(handles.date,'String',date)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes about_plotter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = about_plotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structu

varargout{1} = handles.output;
