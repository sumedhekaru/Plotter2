function varargout = video_analyzer(varargin)
% VIDEO_ANALYZER MATLAB code for video_analyzer.fig
%      VIDEO_ANALYZER, by itself, creates a new VIDEO_ANALYZER or raises the existing
%      singleton*.
%
%      H = VIDEO_ANALYZER returns the handle to a new VIDEO_ANALYZER or the handle to
%      the existing singleton*.
%
%      VIDEO_ANALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIDEO_ANALYZER.M with the given input arguments.
%
%      VIDEO_ANALYZER('Property','Value',...) creates a new VIDEO_ANALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before video_analyzer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to video_analyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help video_analyzer

% Last Modified by GUIDE v2.5 29-Jul-2014 19:46:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @video_analyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @video_analyzer_OutputFcn, ...
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


% --- Executes just before video_analyzer is made visible.
function video_analyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to video_analyzer (see VARARGIN)

% Choose default command line output for video_analyzer
handles.output = hObject;

%Load Phantom DLLS
LoadPhantomLibraries();
RegisterPhantom(true)
clc

set(handles.hsv_analyzer,'Position',[1080 31 63 30])

% Delete previos video windos
try; delte(500); end;

b.plot_e_field  = 1;
b.file_name     = 'HSV file name';
b.dir           = '';
b.begining      = '';
b.end0           = '';
b.start_frame   = '';
b.stop_frame    = '';
b.align_frame   = '';
b.align_time    = '';
b.bg_frame      = '';
b.n_of_e_frames = 11;
b.color_map     = 9;
b.color_map_str = 'Gray';
b.graphs        = 1;
b.cf_sli_posi = 0;
b.FF_RQ_frames  = 100;
b.current_frame = '';
b.save_fn = 0;
b.save_pn = 0;

b.t1 = ''; % Current frame times
b.t2 = '';
b.T1 = ''; % Total time range (from start to stop)
b.T2 = '';

b.ldar = 0;
b.linet = 0;
b.pbfa = 0;
b.hsvp = 0;
b.pre_plot = 0;
b.post_plot = 0;

% Location values
b.cam_x = -20636.6;
b.cam_y = 1569.3;
b.cam_z = 0;
b.ref_x = '';
b.ref_y = '';
b.ref_z = '';
b.scr_x = '';
b.scr_y = '';

b.vh    = 500; % Handle for video
b.gh    = 501; % Handle for plotter graphs


try
    % check whether user asking for new plotter
    if varargin{1}==1
        % Do nothing        
    end
catch
    try
        % If not try to load last GUI data

        b = open([tempdir 'hsv_analyser_last_gui_data.mat']);
        b.save_fn = 0;
    catch
        % Do nothing
    end
end


handles.b       = b;
handles.lines   = [nan nan];

% Video setting file
handles.v_set = open('video_settings.mat');



global c
c.is_paused = false;
handles.ph.cineHandle = NaN;

handles = set_values(handles);


% Get sensor_settings
handles.sen_set = open('sensor_setting.mat');


% Load plotter2 data
try
    h=guidata(findall(0,'Tag','plotter2'));
    handles.g = h.g;
catch
    fprintf('\nE-field plotting parameters are \ncoming from Plotter2. \n\nPlease run plotter2 first.\n')
end

handles.temp.gt1 = handles.g.t1;
handles.temp.gt2 = handles.g.t2;

% Open cine file
handles = set_cine(handles);

try
    allign_time_Callback(hObject, eventdata, handles)
end

% Add aditional tools to figures
% handles = add_tools(handles);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes video_analyzer wait for user response (see UIRESUME)
% uiwait(handles.hsv_analyzer);


% --- Outputs from this function are returned to the command line.
function varargout = video_analyzer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_fn_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fn as text
%        str2double(get(hObject,'String')) returns contents of edit_fn as a double


% --- Executes during object creation, after setting all properties.
function edit_fn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in choose_fn.
function choose_fn_Callback(hObject, eventdata, handles)

[fn,dir] = uigetfile('H:\14Aug2011\*.cine');

if fn~=0
    handles.b.file_name = fn;
    handles.b.dir = dir;
    

    % get cine handle
    [HRES, cineHandle] = PhNewCineFromFile([dir fn]);
    handles.ph.cineHandle = cineHandle;
    
       
    % Get info about the cine
    pFirstIm = libpointer('int32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_FIRSTIMAGENO, pFirstIm)
    handles.b.begining= double(pFirstIm.Value);
    pImCount = libpointer('uint32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_IMAGECOUNT, pImCount);
    handles.b.end0 = (handles.b.begining + double(pImCount.Value) - 1);
    handles.b.start_frame  = handles.b.begining;
    handles.b.stop_frame   = handles.b.end0;
    handles.b.bg_frame     = handles.b.begining;
    handles.b.current_frame = handles.b.begining;
   
    
    
    %get cine image buffer size
    pInfVal = libpointer('uint32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_MAXIMGSIZE, pInfVal);
    imgSizeInBytes = pInfVal.Value;
    %The image flip for GetCineImage function is inhibated.
    pInfVal = libpointer('int32Ptr',false);
    PhSetCineInfo(cineHandle, PhFileConst.GCI_VFLIPVIEWACTIVE, pInfVal);
    %Create the image reange to be readed
    imgRange = get(libstruct('tagIMRANGE'));
    %take one image at imageNo
    imgRange.First = handles.b.start_frame;
    imgRange.Cnt = 1;
    
    handles.ph.cineHandle = cineHandle;
    handles.ph.imgRange = imgRange;
    handles.ph.imgSizeInBytes = imgSizeInBytes;
    %clc
    
    % Get Frame Rate
    pInfVal = libpointer('uint32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_FRAMERATE, pInfVal); 
    handles.ph.frame_rate = double(pInfVal.Value);
    
    
    % Get time
    pInfVal = libpointer('uint32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_TRIGTIMESEC , pInfVal);   
    sec= pInfVal.Value;
    
    sec = double(sec);
    
    pInfVal = libpointer('uint32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_TRIGTIMEFR , pInfVal);    
    secf= pInfVal.Value;
    
    secf = double(secf);
    
    binval = dec2bin(secf);
    
    binval2 = binval(1:end-2);
    binval2 = [binval2 '00']; % Last two digits contains IRIG sycro info
    
    secf = bin2dec(binval2);
    
    sec = sec+secf/2^32;
    
    unix_epoch = datenum(1970,1,1,0,0,0);
    matlab_time = sec./86400 + unix_epoch;
    trig_time = datevec(matlab_time);
    % trig_time = datenum(trig_time - [trig_time(1:3) 0 0 0])
    handles.ph.trig_time = trig_time;
        
%     pInfVal = libpointer('int8Ptr',0);
%     PhGetCineInfo(cineHandle, PhFileConst.GCI_LENSDESCRIPTION   , pInfVal)  
%     secf= double(pInfVal.Value)
%     
     set_values(handles);
     
     % Show the first frame
     handles = show(handles,handles.b.current_frame);
     
     % Set up image (one time)
     handles = setup_image(handles);
end

guidata(hObject, handles);



function allign_frame_Callback(hObject, eventdata, handles)

temp = str2double(get(hObject,'String'));

% If temp is outside the frame range, let's exit
if temp < handles.b.begining || temp > handles.b.end0
    errordlg('Allign Frame given is out of range.','Align Frame Error')
    return
end

handles.b.align_frame = temp;
handles.b.current_frame = temp;

% Let's show the frame
if ~isempty(handles.b.align_time)
    handles = show(handles,temp);
end

% Let's load ldar data
handles = load_ldar(handles);

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function allign_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allign_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function allign_time_Callback(hObject, eventdata, handles)

try
    at= str2double(get(hObject,'String')) ;
    handles.b.align_time = at;
catch
    at = handles.b.align_time;
end

% Let's load the nevest plotter values
try
    h=guidata(findall(0,'Tag','plotter2'));
catch
    errordlg(sprintf('E-field plotting parameters are \ncoming from Plotter2. \n\nPlease run plotter2 first.\n'),'Plotter2 not found')
    return
end

handles.g = h.g;
handles.cal_info = h.cal_info;

% Let's find out whether the allign time is within pltter time range
if at < h.g.t1min || at > h.g.t2max
    errordlg(sprintf('Allign time must be within %is - %is\nIf not choose different date/time from plotter2',h.g.t1min,h.g.t2max),'Allign Time Error')  
    return
end



% Let's plot all the data that can ve within trigger
vtl = 3;    % Vedio trigger length (time in secconds)

% Lower time
if handles.g.t1min < (at - vtl)
    handles.g.t1 = at - vtl;
else
    handles.g.t1 = handles.g.t1min;
end



% Upper time
if handles.g.t2max > (at + vtl)
    handles.g.t2 = at + vtl;
else
    handles.g.t2 = handles.g.t2max;
end

% This will load the LDAR data
handles = load_ldar(handles);


handles = plot_all_for_video(handles);

handles = show(handles,handles.b.current_frame);

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function allign_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allign_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
%clc
i = handles.b.current_frame;
sf = handles.b.start_frame;
total = handles.b.stop_frame - sf;


set(handles.hsv_analyzer,'Units','Pixels');
pos = get(handles.hsv_analyzer,'Outerposition');

wbh = waitbar((i-sf)/total,'Please close this to stop playing',...
    'name','Play','Units','Pixels','Outerposition',[pos(1)+pos(3)/2-183 pos(2)-103 366 103]);


while i < handles.b.stop_frame
    handles = show(handles,i);
    
    try
        waitbar((i-sf)/total,wbh)
    catch
        handles.b.current_frame = i;
        handles = show(handles,i);
        guidata(hObject,handles)
        return
    end
    
    i = i+1;
end

handles.b.current_frame = i;
handles = show(handles,i);
delete(wbh)
guidata(hObject,handles)


% --- Executes on selection change in color_map.
function color_map_Callback(hObject, eventdata, handles)

contents = cellstr(get(hObject,'String'));
handles.b.color_map = get(hObject,'Value');
handles.b.color_map_str = contents{get(hObject,'Value')};
figure(handles.b.vh)
colormap(handles.b.color_map_str)
handles = show(handles, handles.b.current_frame);
handles = setup_image(handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function color_map_CreateFcn(hObject, eventdata, handles)
% hObject    handle to color_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in graph_selection.
function graph_selection_Callback(hObject, eventdata, handles)
    
    handles.b.graphs = get(hObject,'Value');
    handles = show(handles, handles.b.current_frame);
    guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function graph_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to graph_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bg_frame_Callback(hObject, eventdata, handles)






% Hints: get(hObject,'String') returns contents of bg_frame as text
%        str2double(get(hObject,'String')) returns contents of bg_frame as a double


% --- Executes during object creation, after setting all properties.
function bg_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bg_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in next_frame.
function next_frame_Callback(hObject, eventdata, handles)
cf = handles.b.current_frame;

if cf < handles.b.stop_frame
    handles.b.current_frame = cf + 1;
    handles = show(handles,handles.b.current_frame);
    guidata(hObject,handles)
    
end


% --- Executes on button press in frev_frame.
function frev_frame_Callback(hObject, eventdata, handles)
cf = handles.b.current_frame;

if cf > handles.b.start_frame
    handles.b.current_frame =cf - 1;
    handles = show(handles,handles.b.current_frame);
    guidata(hObject,handles)
end

% --- Executes on button press in rewind.
function rewind_Callback(hObject, eventdata, handles)
%clc
i= handles.b.current_frame;
sf = handles.b.start_frame;
total = handles.b.stop_frame -sf;


set(handles.hsv_analyzer,'Units','Pixels');
pos = get(handles.hsv_analyzer,'Outerposition');

wbh = waitbar((i-sf)/total,'Please close this to stop rewinding',...
    'name','Rewind','Units','Pixels','Outerposition',[pos(1)+pos(3)/2-183 pos(2)-103 366 103]);



while i > handles.b.start_frame
    handles = show(handles,i);
    
    try
        waitbar((i-sf)/total,wbh)
    catch
        handles.b.current_frame = i;
        handles = show(handles,i);
        guidata(hObject,handles)
        return
    end
    
    i = i-handles.b.FF_RQ_frames;
end

i = handles.b.start_frame;

handles.b.current_frame = i;
handles = show(handles,i);
delete(wbh)
guidata(hObject,handles)

% --- Executes on button press in fast_foward.
function fast_foward_Callback(hObject, eventdata, handles)
%clc
i = handles.b.current_frame;
sf = handles.b.start_frame;
total = handles.b.stop_frame - sf;

set(handles.hsv_analyzer,'Units','Pixels');
pos = get(handles.hsv_analyzer,'Outerposition');

wbh = waitbar((i-sf)/total,'Please close this to stop playing',...
    'name','FastFoward','Units','Pixels','Outerposition',[pos(1)+pos(3)/2-183 pos(2)-103 366 103]);



while i < handles.b.stop_frame
    handles = show(handles,i);
    
    try
        waitbar((i-sf)/total,wbh)
    catch
        handles.b.current_frame = i;
        handles = show(handles,i);
        guidata(hObject,handles)
        return
    end
    
    i = i+handles.b.FF_RQ_frames;
end

i = handles.b.stop_frame;

handles.b.current_frame = i;
handles = show(handles,i);
delete(wbh)
guidata(hObject,handles)


% --- Executes on button press in reverse.
function reverse_Callback(hObject, eventdata, handles)
%clc
i = handles.b.current_frame;
sf = handles.b.start_frame;
total = handles.b.stop_frame - sf;

set(handles.hsv_analyzer,'Units','Pixels');
pos = get(handles.hsv_analyzer,'Outerposition');

wbh = waitbar((i-sf)/total,'Please close this to stop reverse play',...
    'name','Reverse','Units','Pixels','Outerposition',[pos(1)+pos(3)/2-183 pos(2)-103 366 103]);


while i > handles.b.start_frame
    handles = show(handles,i);
    
    try
        waitbar((i-sf)/total,wbh)
    catch
        handles.b.current_frame = i;
        handles = show(handles,i);
        guidata(hObject,handles)
        return
    end
    
    i = i-1;
end

handles.b.current_frame = i;
handles = show(handles,i);
delete(wbh)
guidata(hObject,handles)




% --- Executes on button press in plot_e_field.
function plot_e_field_Callback(hObject, eventdata, handles)
    handles.b.plot_e_field = get(hObject,'Value');
    
    guidata(hObject,handles)


% --- Executes on slider movement.
function current_frame_slider_Callback(hObject, eventdata, handles)

st=handles.b.start_frame;
en=handles.b.stop_frame;

if ~isempty(st) || ~isempty(en)
    val = get(hObject,'Value');
    handles.b.current_frame = round(st + val*(en - st));
    
    set(handles.curr_frame_text,'String',num2str(handles.b.current_frame))
    
    % if allign parameters are sat, let's plot it
    if ~isempty(handles.b.align_frame) || ~isempty(handles.b.align_time)
        handles = show(handles,handles.b.current_frame);
    end
    
    guidata(hObject,handles)
    
end

    

%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function current_frame_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_frame_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)

b= handles.b;

[fn,pn] = uigetfile([b.dir b.file_name(1:end-4) 'mat']);

filename=[pn fn];

if fn==0
    % File did not selected, do nothing!
else
    b=open(filename);
    
    % this will be necessary if user rename the file
    b.save_fn = fn;
    b.save_pn = pn;
    
    b.isSaved = 1;
    
    handles.b = b;
    
    % Set values
    handles = set_values(handles);
    
 
    % Open cine file
    handles = set_cine(handles);
    
   
    allign_time_Callback(hObject, eventdata, handles)
    
    guidata(hObject,handles)
    
end




% --------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)

b=handles.b;

if b.save_fn(1)==0
    % suggested filename
    
    [b.save_fn,b.save_pn] = uiputfile('*.mat','Save..',[b.dir b.file_name(1:end-4) 'mat']);
    if b.save_fn(1)~=0
        fn=[b.save_pn b.save_fn];
        save(fn,'-struct', 'b')
        set(handles.hsv_analyzer,'Name',['HSV Analyzer -' b.save_fn])
        b.isSaved=1;
    end
else
    fn=[b.save_pn b.save_fn];
    save(fn,'-struct', 'b')
    b.isSaved=1;
end
handles.b=b;
guidata(hObject, handles);


% --------------------------------------------------------------------
function save_as_Callback(hObject, eventdata, handles)
b=handles.b;

[fn,pn] = uiputfile('*.mat','Save As..',[b.dir b.file_name(1:end-4) 'mat']);
if fn(1)~=0
    b.save_fn=fn;
    b.save_pn=pn;
    fn=[pn fn];
    save(fn,'-struct', 'b')
    set(handles.hsv_analyzer,'Name',['HSV Analyzer -' b.save_fn])
    b.isSaved=1;
    handles.sen_set.working_dir = b.save_pn;
    handles.b=b;
    guidata(hObject, handles);
end


% --------------------------------------------------------------------
function tools_Callback(hObject, eventdata, handles)
% hObject    handle to tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function use_meters_Callback(hObject, eventdata, handles)
v_set = handles.v_set;
if strcmp(get(gcbo, 'Checked'),'on')
    set(gcbo, 'Checked', 'off');
    v_set.useMeters = 0;
else
    set(gcbo, 'Checked', 'on');
    v_set.useMeters = 1;
end

save('video_settings.mat','-Struct','v_set')
handles.v_set = v_set;
axisChange(handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function find_bright_frames_Callback(hObject, eventdata, handles)

b=handles.b;

% Let's see whther the file name is given
if isempty(handles.b.dir)
    errordlg('Please choose a cine file first','Cine File Error')
    return
end


% Does user need to crop

[a, bt] = settingsdlg(...
    'Description', 'This will analyze and find bright frames of HSV',... 
    'title'      , 'Bright Frames',...
    'separator'  , 'Number of frames to omit',...
    {'Omit';'N'}, 10,...
    'separator'  , ' ',...
    {'Crop the image'; 'crop'}, true);        %240
clc


if ~strcmp(bt,'ok')
    return
end
%clc
N = a.N;

if a.crop == true
    fg = figure(handles.b.vh);
    [temp temp tem rect] = imcrop();
else
    matlabIm = get_image(handles,handles.b.start_frame);
    [R B] = size(matlabIm);
    rect = [1 1 R B];
end

% Manual rectangle
rect = [181 1 79 120];

wbh = waitbar(0,'Please Wait','name','Loking for bright frames');

tic


i_array = handles.b.start_frame:N:handles.b.stop_frame;
mean_array = i_array-i_array;
leng = length(i_array);

for i=1:leng
    
    try 
        x = i/leng;
        msg = sprintf('Analysing frame %i of %i',i,leng);
        waitbar(x,wbh,msg)
    catch
        fprintf('Finding Return Strock stopped by user\n')
        return
    end

    matlabIm = get_image(handles,i_array(i));
    matlabIm = imcrop(matlabIm, rect);
    %matlabIm = matlabIm(a.L:a.R-1,a.T:a.B-1);
    %x = mean2(matlabIm)
    %figure
    %image(matlabIm,'CDataMapping','scaled')
    %imagesc(matlabIm,[0, 65535])
    mean_array(i) = mean2(matlabIm);

end

try
    delete(wbh)
end

% Normalze data
mean_array = double(mean_array - min(mean_array));
mean_array = mean_array / max(mean_array);

% save_data.intensity = mean_array;
% save_data.frames = i_array;
% 
% save('C:\Users\Sumedhe\Desktop\video_intensity_data2.mat','-Struct','save_data')

toc

[peakLoc ~] = peakfinder(mean_array,0.1);

figure
plot(i_array,mean_array)
assignin('base','frames',i_array)
assignin('base','intensity',mean_array)
hold all
plot(i_array(peakLoc),mean_array(peakLoc),'pr','MarkerFaceColor','r')
xlabel('Frame Number')
ylabel('Normalize Average Intensity')
set(0,'defaulttextinterpreter','none')
title(handles.b.file_name)
set(0,'defaulttextinterpreter','default')

% Display peak info on matlab
fprintf('\nBright Frames are :\n')
for i=1:length(peakLoc)
    fprintf('\t\t%6.0f\n',i_array(peakLoc(i)))
end





function handles = set_values(handles)

set(handles.plot_e_field,'Value',handles.b.plot_e_field)
set(handles.edit_fn,'String',handles.b.file_name)
set(handles.beg,'String',['Begining: ' num2str(handles.b.begining)])
set(handles.endf,'String',['End: ' num2str(handles.b.end0)])
set(handles.start_frame,'String',handles.b.start_frame)
set(handles.stop_frame,'String',handles.b.stop_frame)
set(handles.allign_frame,'String',handles.b.align_frame)
set(handles.allign_time,'String',sprintf('%.7f',handles.b.align_time))
set(handles.bg_frame,'String',handles.b.bg_frame)
set(handles.n_of_e_frames,'String',handles.b.n_of_e_frames)
set(handles.color_map,'Value',handles.b.color_map)
set(handles.graph_selection,'Value',handles.b.graphs)
set(handles.curr_frame_text,'String',num2str(handles.b.current_frame))
set(handles.current_frame_slider,'value',handles.b.cf_sli_posi)
set(handles.hsv_analyzer,'Name',['HSV Analyzer -' handles.b.save_fn])
set(handles.ldar,'Value',handles.b.ldar)
set(handles.linet,'Value',handles.b.linet)
set(handles.pbfa,'Value',handles.b.pbfa)
set(handles.hsvp,'Value',handles.b.hsvp)
set(handles.pre_plot,'Value',handles.b.pre_plot)
set(handles.post_plot,'Value',handles.b.post_plot)

% setings based values
if handles.v_set.useMeters
    set(handles.use_meters, 'Checked', 'on');
else
    set(handles.use_meters, 'Checked', 'off');
end





function matlabIm = get_image(handles,imageNu)

handles.ph.imgRange.First = imageNu;

%Read the cine image into the buffer
[~, unshiftedIm, imgHeader] = PhGetCineImage(handles.ph.cineHandle,...
                                                handles.ph.imgRange, ...
                                                handles.ph.imgSizeInBytes);
pImCount = libpointer('int32Ptr',1);
                                                                                      

% Transform 1D image pixels to 1D/3D image pixels to be used with MATLAB
[unshiftedIm] = ExtractImageMatrixFromImageBuffer(unshiftedIm, imgHeader);

bps = GetEffectiveBitsFromIH(imgHeader);
[matlabIm, ~] = ConstructMatlabImage(unshiftedIm, imgHeader.biWidth, imgHeader.biHeight, 1, bps);

%data.img = matlabIm;
%save('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F18\background.mat','-struct','data')



function handles=plot_all_for_video(handles)

% Le't close if there is a current graph
try
    delete(handles.b.gh)
end


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
if nnz(g.chgraphs)==0 && nnz(g.lpgraphs)==0 && g.ldar==0 && ...
        g.linet==0 && g.pbfa ==0 && g.nldn==0 && g.dedt == 0 && g.dedt_trig == 0
    errordlg('No graphs were selected to plot!','Graph selection Error','modal')
    uiwait
    return
end

% Lets setup calibration values before plot
handles=calibration(handles);

plot_all3(handles.g)
% set(gca,'YLimMode','auto') 
grid off
handles.b.gh = gcf;
set(handles.b.gh,'Position',[600 378 560 420])
% set(handles.b.gh,'visible','off')





function start_frame_Callback(hObject, eventdata, handles)

temp = str2double(get(hObject,'String'));

% If temp is outside the frame range, let's exit
if temp < handles.b.begining || temp > handles.b.end0
    errordlg('Begining Frame given is out of range.','Begining Frame Error')
    set(handles.start_frame,'String',handles.b.begining)
    handles.b.start_frame = handles.b.begining;
    guidata(hObject, handles);
    return
end

handles.b.start_frame = temp;
handles.b.current_frame = temp;
handles = show(handles,temp);

handles = load_ldar(handles);

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function start_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stop_frame_Callback(hObject, eventdata, handles)

temp = str2double(get(hObject,'String'));

% If temp is outside the frame range, let's exit
if temp < handles.b.begining || temp > handles.b.end0
    errordlg('Ending Frame given is out of range.','Ending Frame Error')
    set(handles.stop_frame,'String',handles.b.end0)
    handles.b.stop_frame = handles.b.end0;
    guidata(hObject, handles);
    return
end

handles.b.stop_frame = temp;
handles = load_ldar(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function stop_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stop_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = set_xrange(handles,img)

% f = handles.ph.frame_rate;
% n = handles.b.n_of_e_frames;
% T = handles.b.align_time;


if ~ishandle(handles.b.gh)
    % The graph not found
    disp('Warning: E-field graph not found')
end

% figure out range of the graphs
% 
% t1 = T + (img - handles.b.align_frame - (n+1)/2)/f;
% t2 = t1 + n/f;

try
    delete(handles.lines)
end

figure(handles.b.gh)

axh = gca;
set(axh,'XLim',[handles.b.window_t1 handles.b.window_t2])


% Grids for middle frame
t2 = handles.b.t1;
t1 = handles.b.t2;

x2 = [t2 t2 ] ;
y = get(axh,'YLim') ;
x1 = [t1 t1];

hold all
handles.lines(1) = plot(x1,y,'Color','r','LineStyle','--','LineWidth',2);
%handles.lines(2) = plot(x1,y,'Color','b','LineStyle',':','LineWidth',2);

handles.lines(2) = plot(x2,y,'Color','r','LineStyle','--','LineWidth',2);
%handles.lines(4) = plot(x2,y,'Color','b','LineStyle',':','LineWidth',2);

% Set xlebel to show time range of middle frame
    str = sprintf('Time {%0.7f - %0.7f}s',t1,t2);
xlabel(str)

% % Let's save t1 and t2 for future use
% handles.b.t1 = t1;
% handles.b.t2 = t2;


    

function handles = show(handles,N)
% This will display N th image
matlabIm = get_image(handles,N);

fg = figure(handles.b.vh);


cla
switch handles.b.graphs
    case 1
         %image(double(matlabIm),'CDataMapping','scaled');
         %matlabIm = imresize(matlabIm,5);
         %imagesc(matlabIm,[1697/16128*65535 9219/(16128-2353)*65535])
         imagesc(matlabIm,[0, 65535])
         %image(double(matlabIm));
         
     case 2
        matlabIm = imcomplement(matlabIm);
        image(matlabIm,'CDataMapping','scaled');

        %imagesc(matlabIm,[0, 65535])
    case 3
        bg = get_image(handles,handles.b.bg_frame);
        matlabIm = double(matlabIm) - double(bg);
        %image(matlabIm,'CDataMapping','scaled');
        imagesc(matlabIm)
    case 4
        bg = get_image(handles,handles.b.bg_frame);
        matlabIm = double(matlabIm) - double(bg);
        matlabIm = imcomplement(matlabIm);
        image(matlabIm,'CDataMapping','scaled'); 
        %imagesc(matlabIm,[0, 65535])
end
        
daspect([1 1 1])
title(num2str(N))
hold all


% Calculate middle frame times
f = handles.ph.frame_rate;
n = handles.b.n_of_e_frames;
T = handles.b.align_time;


% Time range for the veiw
t1 = T + (N - handles.b.align_frame - (n+1)/2)/f;
t2 = t1 + n/f;

handles.b.window_t1 = t1;
handles.b.window_t2 = t2;


% Grids for middle frame
t2 = t2-(n-1)/2/f;
t1 = t2-1/f;

handles.b.t1 = t1;
handles.b.t2 = t2;

% Plotting LDAR/PBFA/LINET
if handles.b.ldar; handles = plot_ldar(handles); end;
if handles.b.pbfa; handles = plot_pbfa(handles); end;
if handles.b.hsvp; handles = plot_hsvp(handles); end;

% xlim([0 320])
% ylim([-150 240])
% zticks = 0:1000:7000;
% [scrX scrY] = ldar2scr(zeros(size(zticks))-12816,zeros(size(zticks))+1254,zticks,handles)
% [x y z] = scr2ldar(sscr,scrY,handles);
% set(gca,'YTickLabel',round(z))




% E-field graph arrange
if handles.b.plot_e_field && ishandle(handles.b.gh)
    handles = set_xrange(handles,N);
end

% Update slider
update_slider(handles,N)






function n_of_e_frames_Callback(hObject, eventdata, handles)
temp = str2double(get(hObject,'String'));

if mod(temp,2) == 0
    errordlg('Number of frames must be an ODD number','Error')
else
    handles.b.n_of_e_frames = temp;
end

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function n_of_e_frames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_of_e_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function update_slider(handles,cf)
    st = handles.b.start_frame;
    en = handles.b.stop_frame;
    % cf = handles.b.current_frame;
    set(handles.curr_frame_text,'String',num2str(cf))
    val = (cf-st)/(en-st);
    if val <= 1 && val >= 0
        set(handles.current_frame_slider,'Val',val)
    else
        disp(val)
    end
 
% --------------------------------------------------------------------
function crop_current_image_Callback(hObject, eventdata, handles)

if ishandle(handles.b.vh)
    [I,rect]= imcrop(handles.b.vh);
    % Show Cropped Image
    if ~isempty(I)
        fg = figure;
        set(fg,'Name','Cropped Image')
        %colormap(handles.b.color_map_str)
        colormap(gray(2^8))
        image(I,'CDataMapping','scaled')
            %'XData',[rect(1) rect(3)],'YData',[rect(3) rect(4)]); 
        title(['Cropped Image of image #: ' num2str(handles.b.current_frame)])
        daspect([1 1 1])
        
    end
else
    errordlg('Image window is not found','Crop Error')
end

function handles = setup_image(handles)

set(handles.b.vh,'Name','Video Window')
set(handles.b.vh,'Position',[20 378 560 420])
colormap(handles.b.color_map_str)

    
% --------------------------------------------------------------------
function vert_int_Callback(hObject, eventdata, handles)

if ishandle(handles.b.vh)
    %clc
    [I,rect]= imcrop(handles.b.vh);
    % Show Cropped Image
    if ~isempty(I)
        fg = figure;
        subplot(1,2,1)
        set(fg,'Name','Cropped Image')
        colormap(handles.b.color_map_str)
        image(I,'CDataMapping','scaled')
        
            %'XData',[rect(1) rect(3)],'YData',[rect(3) rect(4)]); 
        title(['Cropped Image of image #: ' num2str(handles.b.current_frame)])
        daspect([1 1 1])
        
        [h l]=size(I);
        
        inten = nan(1,h);
        
        for i = 1:h
            inten(i) = mean2(I(i,:));
        end
        
        subplot(1,2,2)
        ma = max(inten);
        mi = min(inten);
        plot((inten-mi)/(ma - mi),1:h)
        set(gca,'YDir','reverse')
        xlabel('Normalized Intensity')
            
        
        
        
    end
else
    errordlg('Image window is not found','Crop Error')
end



function curr_frame_text_Callback(hObject, eventdata, handles)

st=handles.b.start_frame;
en=handles.b.stop_frame;

if ~isempty(st) || ~isempty(en)
    
    cf = str2double(get(hObject,'String'));
    
    if cf < st || cf > en
        errordlg('Frame entered outof range','Current frame error')
        return
    end
    
    handles.b.current_frame = cf;
    
    update_slider(handles,cf)
    
    % if allign parameters are sat, let's plot it
    if ~isempty(handles.b.align_frame) || ~isempty(handles.b.align_time)
        handles = show(handles,handles.b.current_frame);
    end
    
    guidata(hObject,handles)
    
end



% hObject    handle to curr_frame_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of curr_frame_text as text
%        str2double(get(hObject,'String')) returns contents of curr_frame_text as a double


% --- Executes during object creation, after setting all properties.
function curr_frame_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curr_frame_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Settings_Callback(hObject, eventdata, handles)
% hObject    handle to Settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function location_reference_Callback(hObject, eventdata, handles)

b = handles.b;

[settings, button] = settingsdlg(...
    'title'      , 'Location Settings',...
    'separator'  , 'Camera Location',...
    {'X'; 'cam_x'}, b.cam_x,...
    {'Y'; 'cam_y'}, b.cam_y,...
    {'Z'; 'cam_z'}, b.cam_z,...
 'separator'  , 'Reference Point',...
    {'X'; 'ref_x'}, b.ref_x,...
    {'Y'; 'ref_y'}, b.ref_y,...
    {'Z'; 'ref_z'}, b.ref_z,...
     'separator'  , 'Screen Reference',...
    {'X'; 'scr_x'}, b.scr_x,...
    {'Y'; 'scr_y'}, b.scr_y);

if strcmp(button,'ok')
    b.cam_x = settings.cam_x;
    b.cam_y = settings.cam_y;
    b.cam_z = settings.cam_z;
    b.ref_x = settings.ref_x;
    b.ref_y = settings.ref_y;
    b.ref_z = settings.ref_z;
    b.scr_x = settings.scr_x;
    b.scr_y = settings.scr_y;
    
    handles.b = b;
    guidata(hObject,handles)
    
end







function handles = set_cine(handles)
%clc

fn = [handles.b.dir handles.b.file_name];

if ~exist(fn,'file')
    % file name is not there, no point to continue
    if ~strcmp('HSV file name',fn)
        errordlg(['File "' fn '" not found'])
    end
    return
end  
    


    fn = handles.b.file_name ;
    dir = handles.b.dir;
    

    % get cine handle
    [HRES, cineHandle] = PhNewCineFromFile([dir fn]);
    handles.ph.cineHandle = cineHandle;
    
       
    % Get info about the cine
    pFirstIm = libpointer('int32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_FIRSTIMAGENO, pFirstIm)
    handles.b.begining= double(pFirstIm.Value);
    pImCount = libpointer('uint32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_IMAGECOUNT, pImCount);
    handles.b.end0 = (handles.b.begining + double(pImCount.Value) - 1);
    %handles.b.start_frame  = handles.b.begining;
    %handles.b.stop_frame   = handles.b.end0;
    handles.b.bg_frame     = handles.b.begining;
    %handles.b.current_frame = handles.b.begining;
   
    
    
    %get cine image buffer size
    pInfVal = libpointer('uint32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_MAXIMGSIZE, pInfVal);
    imgSizeInBytes = pInfVal.Value;
    %The image flip for GetCineImage function is inhibated.
    pInfVal = libpointer('int32Ptr',false);
    PhSetCineInfo(cineHandle, PhFileConst.GCI_VFLIPVIEWACTIVE, pInfVal);
    %Create the image reange to be readed
    imgRange = get(libstruct('tagIMRANGE'));
    %take one image at imageNo
    imgRange.First = handles.b.start_frame;
    imgRange.Cnt = 1;
    
    handles.ph.cineHandle = cineHandle;
    handles.ph.imgRange = imgRange;
    handles.ph.imgSizeInBytes = imgSizeInBytes;
    %clc
    
    % Get Frame Rate
    pInfVal = libpointer('uint32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_FRAMERATE, pInfVal); 
    handles.ph.frame_rate = double(pInfVal.Value);
    
    
    % Get time
    pInfVal = libpointer('uint32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_TRIGTIMESEC , pInfVal);   
    sec= pInfVal.Value;
    
    sec = double(sec);
    
    pInfVal = libpointer('uint32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_TRIGTIMEFR , pInfVal);    
    secf= pInfVal.Value;
    
    secf = double(secf);
    
    binval = dec2bin(secf);
    
    binval2 = binval(1:end-2);
    binval2 = [binval2 '00']; % Last two digits contains IRIG sycro info
    
    secf = bin2dec(binval2);
    
    sec = sec+secf/2^32;
    
    unix_epoch = datenum(1970,1,1,0,0,0);
    matlab_time = sec./86400 + unix_epoch;
    trig_time = datevec(matlab_time);
    % trig_time = datenum(trig_time - [trig_time(1:3) 0 0 0])
    handles.ph.trig_time = trig_time;
        
%     pInfVal = libpointer('int8Ptr',0);
%     PhGetCineInfo(cineHandle, PhFileConst.GCI_LENSDESCRIPTION   , pInfVal)  
%     secf= double(pInfVal.Value)
%     


     set_values(handles);
     
    % This will load the LDAR data
    handles = load_ldar(handles);
    
    
    %%%%

    
    
    
     
     % Show the first frame
     handles = show(handles,handles.b.current_frame);
     handles = add_tools(handles);
     
     % Set up image (one time)
     handles = setup_image(handles);
     
     if handles.v_set.useMeters
        axisChange(handles)
     end


% --- Executes on button press in ldar.
function ldar_Callback(hObject, eventdata, handles)
    handles.b.ldar = get(hObject,'Value');
    figure(handles.b.vh)
    handles = show(handles,handles.b.current_frame);
    guidata(hObject,handles)


% --- Executes on button press in linet.
function linet_Callback(hObject, eventdata, handles)
    handles.b.linet = get(hObject,'Value');
    
    guidata(hObject,handles)


% --- Executes on button press in pbfa.
function pbfa_Callback(hObject, eventdata, handles)
    handles.b.pbfa = get(hObject,'Value');
    figure(handles.b.vh)
    handles = show(handles,handles.b.current_frame);    
    guidata(hObject,handles)


% --- Executes on button press in hsvp.
function hsvp_Callback(hObject, eventdata, handles)
    handles.b.hsvp = get(hObject,'Value');
    [x y z] = scr2ldar(handles.b.scr_x,handles.b.scr_y,handles);
    handles = show(handles,handles.b.current_frame); 
    guidata(hObject,handles)


% --- Executes on button press in pre_plot.
function pre_plot_Callback(hObject, eventdata, handles)
    handles.b.pre_plot = get(hObject,'Value');
    figure(handles.b.vh)
    handles = show(handles,handles.b.current_frame);
    guidata(hObject,handles)


% --- Executes on button press in post_plot.
function post_plot_Callback(hObject, eventdata, handles)
    handles.b.post_plot = get(hObject,'Value');
    figure(handles.b.vh)
    handles = show(handles,handles.b.current_frame);
    guidata(hObject,handles)


% --- Executes when user attempts to close hsv_analyzer.
function hsv_analyzer_CloseRequestFcn(hObject, eventdata, handles)
    b = handles.b;
    save([tempdir 'hsv_analyser_last_gui_data.mat'],'-struct', 'b')
    try; delete(handles.b.vh); end;
    try; delete(handles.b.gh); end;
    delete(handles.hsv_analyzer)


% --------------------------------------------------------------------
function new_Callback(hObject, eventdata, handles)
    close(handles.hsv_analyzer)
    pause(1)
    video_analyzer(1)
    
    
function handles = load_ldar(handles)
       
    b = handles.b;
    g = handles.g;
    settings = handles.sen_set;
    
    oriLt = handles.temp.gt1;
    oriHt = handles.temp.gt2;
    
    % Do nothing if any of ldar, linet, pbfa, hsvp is not selected
    if ~b.ldar && ~b.linet && ~b.pbfa && ~b.hsvp
        %disp('Nothing to do')
        return
    end
    
    
    % Do nothing if allign_time, allign_fame, start_frame or stop_frame
    
    if isempty(b.align_time) || isempty(b.align_frame) || ...
            isempty(b.start_frame) || isempty(b.stop_frame)
        disp('Do nothing')
        return
    end
    
    if isnan(b.align_time) || isnan(b.align_frame) || ...
            isnan(b.start_frame) || isnan(b.stop_frame)
        %disp('Do nothing')
        return
    end
    
    % Let's calculate minimum and maximu times for plots
     
    f = handles.ph.frame_rate;
    T = b.align_time;


    % figure out range of the graphs
    b.T1 = T - (b.align_frame - b.start_frame)/f;
    b.T2 = b.T1 + (b.stop_frame - b.start_frame)/f;
    
    absant_fn = '';
    handles.ldard =[];
    handles.pbfad =[];
    
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
    
    % Let's load LDAR2 data
    if b.ldar
        
        dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
            -datenum(str2double(g.YYYY),0,0));
        
        ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
            settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);
        
        if exist(ldar_fn,'file')
            % If file exsist let's load ldar data
            [CG,CAL,DLS]=ldarExtract2(ldar_fn,oriLt,oriHt,str2double(settings.ldar_r),...
                    b.cam_x,b.cam_y,b.cam_z,0);
                                
            handles.ldard(:,1) = DLS(:,10);
            handles.ldard(:,2) = DLS(:,6);
            handles.ldard(:,3) = DLS(:,7);
            handles.ldard(:,4) = DLS(:,8);
            
            
            handles.cglssd(:,1) = CG(:,10);
            handles.cglssd(:,2) = CG(:,6);
            handles.cglssd(:,3) = CG(:,7);
            handles.cglssd(:,4) = CG(:,8);
            
            
            [N_x, N_y] = ldar2scr(handles.ldard(:,2),handles.ldard(:,3),handles.ldard(:,4),handles);
          
            handles.ldard(:,5) = N_x;
            handles.ldard(:,6) = N_y;
            
            [N_x, N_y] = ldar2scr(handles.cglssd(:,2),handles.cglssd(:,3),handles.cglssd(:,4),handles);
            
            handles.cglssd(:,5) = N_x;
            handles.cglssd(:,6) = N_y;
            
            % Lets sort raws considering time
            handles.ldard = sortrows(handles.ldard,1);
            handles.cglssd = sortrows(handles.cglssd,1);
            
        else
            % If file is not exist don't store the file name
            absant_fn=[absant_fn ldar_fn];
        end
    end
    
    % Let's load pbfa data
    if b.pbfa
        
        pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
            settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
            g.MM{:},g.DD{:},g.hh,ext);
        
        if exist(pbfa_fn,'file')
            % If the file exisit, let's load data
            %t1 = g.t1
            %t2 = g.t2
            
            
            PBFA=pbfaExtract(pbfa_fn,oriLt,oriHt,str2double(settings.ldar_r),...
                b.cam_x,b.cam_y,b.cam_z);
            
            handles.pbfad(:,1) = PBFA(:,6);
            handles.pbfad(:,2) = PBFA(:,3);
            handles.pbfad(:,3) = PBFA(:,4);
            handles.pbfad(:,4) = PBFA(:,5);
            
            
            [N_x, N_y] = ldar2scr(handles.pbfad(:,2),handles.pbfad(:,3),handles.pbfad(:,4),handles);
                        
            handles.pbfad(:,5) = N_x;
            handles.pbfad(:,6) = N_y;
            
                        
            % Lets sort raws considering time
            handles.pbfad = sortrows(handles.pbfad,1);
                          
        else
            % If file is not exist let store it to let user know
            absant_fn=[absant_fn pbfa_fn];
        end
    end
   
    %sdf
    
    
        % Let's load pbfa data
        %b.hsvp = 0
    if b.hsvp
        
        hsvp_fn = [b.dir b.file_name(1:end-4) 'txt'];
        
        if exist(hsvp_fn,'file')
            % If the file exisit, let's load data
            data = hsvp_extract(hsvp_fn);
            
            handles.hsvp_data = data;
            handles.hsvpd(:,1) = data.frameN;
            handles.hsvpd(:,2) = data.N_x;
            handles.hsvpd(:,3) = data.N_y;
                        
            % Lets sort raws considering time
            handles.hsvpd = sortrows(handles.hsvpd,1);
            
                          
        else
            % If file is not exist let store it to let user know
            absant_fn=[absant_fn hsvp_fn];
        end
    end
    
    
    
    % Should add LINET
    
    % Let user know about missing files
    if ~isempty(absant_fn)
        errordlg(absant_fn,'Files not found!')
    end
  
    
    
   
    
    
function handles = plot_ldar(handles)
% This function will plot all the location points (eventhough it named as
% PLOT_LDAR

% If ldar data not available, just return
try
    temp = handles.ldard;
catch
    return
end;

b = handles.b;

% Let's load LDAR data if has not been loaded yet
try 
    temp = handles.ldard(:,1);
catch
    handles = load_ldar(handles);
end


% LDAR points in the current frame
lol = sum(handles.ldard(:,1) < b.t1)+1;
ul  = sum(handles.ldard(:,1) <= b.t2);

lol2 = sum(handles.cglssd(:,1) < b.t1)+1;
ul2  = sum(handles.cglssd(:,1) <= b.t2);

% Plot Current ldar point
if b.ldar
    plot(handles.ldard(lol:ul,5),handles.ldard(lol:ul,6),'co','MarkerFaceColor','c','MarkerSize',2)
    plot(handles.cglssd(lol2:ul2,5),handles.cglssd(lol2:ul2,6),'gs','MarkerFaceColor','g','MarkerSize',4)
end

% Plot early ldar points
if b.pre_plot
    plot(handles.ldard(1:lol-1,5),handles.ldard(1:lol-1,6),'co','MarkerFaceColor','c','MarkerSize',2)
    plot(handles.cglssd(1:lol2-1,5),handles.cglssd(1:lol2-1,6),'gs','MarkerFaceColor','g','MarkerSize',4)
end

% Plot post ldar points
if b.post_plot
    plot(handles.ldard(ul+1:end,5),handles.ldard(ul+1:end,6),'co','MarkerFaceColor','c','MarkerSize',2)
    plot(handles.cglssd(ul2+1:end,5),handles.cglssd(ul2+1:end,6),'gs','MarkerFaceColor','g','MarkerSize',4)
end



function handles = plot_pbfa(handles)
% This function will plot all the location points (eventhough it named as
% PLOT_LDAR

% If ldar data not available, just return
try
    temp = handles.pbfad;
catch
    return
end;

b = handles.b;

% Let's load LDAR data if has not been loaded yet
try 
    temp = handles.pbfad(:,1);
catch
    handles = load_ldar(handles);
end

% LDAR points in the current frame
lol = sum(handles.pbfad(:,1) < b.t1)+1;
ul  = sum(handles.pbfad(:,1) <= b.t2);

% Plot Current ldar point
if b.pbfa
    plot(handles.pbfad(lol:ul,5),handles.pbfad(lol:ul,6),'ro','MarkerFaceColor','r','MarkerSize',2)
end

% Plot early ldar points
if b.pre_plot
    plot(handles.pbfad(1:lol-1,5),handles.pbfad(1:lol-1,6),'ro','MarkerFaceColor','r','MarkerSize',2)
end

% Plot post ldar points
if b.post_plot
    plot(handles.pbfad(ul+1:end,5),handles.pbfad(ul+1:end,6),'ro','MarkerFaceColor','r','MarkerSize',2)
end

       
function [N_x,N_y]=ldar2scr(x,y,z,handles)
b=handles.b;

f = 8.0e-3;  %Focal length
p = 20.0e-6; %Fixel Size

N_x=b.scr_x-(f*tan(atan((y-b.cam_y)./(x-b.cam_x ))- atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x )) ))/p;

N_y=b.scr_y - (z - b.cam_z).*f./(p.*sqrt((x-b.cam_x ).^2+(y-b.cam_y ).^2 ))...
    ./sqrt(1+(tan(atan((y-b.cam_y)./(x-b.cam_x ))- atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x )))).^2);

function [x y z] = scr2ldar(N_x,N_y,handles)

b=handles.b;

f = 8.0e-3;  %Focal length
p = 20.0e-6; %Fixel Size


r = sqrt((((b.ref_x-b.cam_x)).^2+(b.ref_y-b.cam_y ).^2 ).*(1+(p^2.*(b.scr_x-N_x ).^2)./f^2 ));

theta = atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x ))+atan(p*(b.scr_x-N_x)/f);

x = b.cam_x + r.*cos(theta);
y = b.cam_y + r.*sin(theta);
z = b.cam_z + r.*((b.scr_y-N_y).*p)./f;


function handles = add_tools(handles)

% tbh = findall(handles.b.vh,'Type','uitoolbar');
% 
% tth = uitoggletool(tbh,'CData',rand(20,20,3),...
%     'ClickedCallback',@curs,...
%     'HandleVisibility','off');

figure(handles.b.vh)
f = uimenu('Label','HSV tools');
uimenu(f,'Label','x,y,z data curser','Callback','hsv_tools(1)','Accelerator','L');
uimenu(f,'Label','Generate HSVP file','Callback','hsv_tools(2)','Accelerator','G');
uimenu(f,'Label','Find distance','Callback','hsv_tools(3)','Accelerator','D')



% --------------------------------------------------------------------
function time2frame_Callback(hObject, eventdata, handles)

time = inputdlg('Enter Time','Time to frame',1,{'100'});

if isempty(time) || strcmp(time , '')
    return
end

try
    time = str2double(time);
    if time < 0 || time > 86400
        errordlg('Input time out of the range','Time to frame')
        return
    end
catch
    errordlg('Error occured','Time to frame')
    return
end

b = handles.b;

frame = round( b.align_frame + ...
    (time - b.align_time)*handles.ph.frame_rate);

clipboard('copy',frame)
msg = sprintf('%0.7f s = %i',time,frame);
msgbox(msg,'Time to frame')






% --------------------------------------------------------------------
function frame2time_Callback(hObject, eventdata, handles)


frame = inputdlg('Enter Frame Nu:','Frame to time',1,{'100'});

if isempty(frame) || strcmp(frame , '')
    return
end

try
    frame = round(str2double(frame));
catch
    errordlg('Error occured','Frame to time')
    return
end

b = handles.b;

time = b.align_time + (frame - b.align_frame)/handles.ph.frame_rate;

%clipboard('copy',time)
msg = sprintf('%i = %0.7f s',frame,time);
msgbox(msg,'Frame to time')


% --------------------------------------------------------------------
function line_figure_creater_Callback(hObject, eventdata, handles)
clc

% [a, bt] = settingsdlg(...
%     'Description', 'This will created a line of images. ',... 
%     'title'      , 'Line Image Creater',...
%     'separator'  , 'Number of frames to omit',...
%     {'Omit';'images'}, '',...
%     'separator'  , ' ',...
%     {'Crop the image'; 'crop'}, true);      %240
% clc
% 
% % xxx = a.images
% 
% if ~strcmp(bt,'ok')
%     return
% end
% 
% % try
%     images = str2num(a.images);
% % catch
% %     images = a.images;
% % end
% 
% 
% if a.crop == true
%     fg = figure(handles.b.vh);
%     [temp temp tem rect] = imcrop();
% else
%     matlabIm = get_image(handles,handles.b.start_frame);
%     [R B] = size(matlabIm);
%     rect = [1 1 R B];
% end

%rect = [48 1 64 120]; % for fig 4,6,8 IB light paper
rect = [48 1 128 240]; % for fig 4,6,8 IB light paper

% Number of vertical fixels between figures
vth = 3;
blank = zeros(rect(4),vth)+65400;
%blank = zeros(rect(4)+1,vth);

images = handles.b.start_frame:2:handles.b.stop_frame;
%rect = [48 1 64 120];
lineIm = [];

%bgf = get_image(handles,-62988);

for i = 1:length(images)
    I = get_image(handles,images(i));    
    I = imcrop(I,rect);
    size(I)
    size(blank)
    lineIm = [lineIm I blank  ];
end 
figure
%colormap(handles.b.color_map_str)
colormap(gray(2^8))
%image(lineIm,'CDataMapping','scaled')
%imagesc(lineIm,[0, 30535])
%lineIm= imcomplement(lineIm);
imagesc(lineIm,[0, 50400])
%imagesc(lineIm,[0, 30535]);

daspect([1 1 1])
set(gca,'ytick',[])
set(gca,'xtick',rect(3)+1:rect(3)+vth+1:rect(3)*length(images))
set(gca,'xticklabel',[])
xl = xlim;
xlim([xl(1) xl(2)-vth])


% Figure numbers
for i = 1:length(images)
    text((i-1)*(rect(3)+vth+1)+5,rect(2)+rect(4)-8,num2str(-images(i)-60000),'color',[1 1 1],'fontsize',12)
end

tools2fig

%title(a.images)


% --------------------------------------------------------------------
function figure_range_balance_Callback(hObject, eventdata, handles)

figure_light_balance(handles);


function handles = plot_hsvp(handles)

% If ldar data not available, just return
b = handles.b;

% if b.hsvp
%     return
% end

try
    temp = handles.hsvpd;
catch
    return
end;



% HSVP for current frame index;
lol = sum(handles.hsvpd(:,1) <  b.current_frame)+1;
ul  = sum(handles.hsvpd(:,1) <=  b.current_frame);


% Plot post hsvp points
if b.post_plot && b.pre_plot
    indx = sum(~isnan(handles.hsvpd(:,2)));
    t = handles.hsvpd(1:indx,1);
    map=colormap(jet);
    miv=min(t);
    mav=max(t);
    %miv = -53804;
    %m1v = -54000;
    %mav = -53497;
    clrstep = (mav-miv)/size(map,1);    
    marker='o';
    
   for nc=1:size(map,1)
        iv = find(t>=miv+(nc-1)*clrstep & t<=miv+nc*clrstep);
        plot(handles.hsvpd(iv',2),handles.hsvpd(iv',3),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',4)
   end
    colormap(gray)
    
    
    
else
    
    % Plot early hsvp points
    if b.pre_plot
        plot(handles.hsvpd(1:lol-1,2),handles.hsvpd(1:lol-1,3),'bo','MarkerFaceColor','b','MarkerSize',5)
    end
    
    if b.post_plot
        plot(handles.hsvpd(ul+1:end,2),handles.hsvpd(ul+1:end,3),'ro','MarkerFaceColor','r','MarkerSize',5)
    end
    
end

% Plot Current hsvp point
if b.hsvp
    plot(handles.hsvpd(lol:ul,2),handles.hsvpd(lol:ul,3),'g+','MarkerFaceColor','g','MarkerSize',5)
end



function axisChange(handles)
fg = figure(handles.b.vh);

if handles.v_set.useMeters    
    
    % Change xy axis to meters
   % xlim([0.5 320])
    %ylim([-150 240])
    
    zticks = 0:1000:12000;
    [scrX scrY] = ldar2scr(zeros(size(zticks))+handles.b.ref_x,...
        zeros(size(zticks))+handles.b.ref_y,zticks,handles);
    [x y z] = scr2ldar(scrX,scrY,handles);
    set(gca,'YTick',fliplr(scrY),'YTickLabel',round(fliplr(z)))
    ylabel('Altitude (m)')
    
    
    [x1 y1] = scr2ldar(0,handles.b.scr_y,handles);
    [x2 y2] = scr2ldar(320,handles.b.scr_y,handles);
    hdist = sqrt((x1-x2)^2+(y1-y2)^2);
    n = floor(hdist/1000);
    xticklabs = 0:1000:n*1000;
    scrX  = 320/hdist*xticklabs;
    set(gca,'XTick',scrX,'XTickLabel',xticklabs);
    xlabel('Distance (m)')
else
    set(gca, 'XTickMode', 'auto', 'YTickMode', 'auto')
    set(gca, 'XTickLabelMode', 'auto', 'YTickLabelMode','auto')
end


% --------------------------------------------------------------------
function create_movie_Callback(hObject, eventdata, handles)
clc
NAr = 0; % Number of arrows
i = handles.b.start_frame;
sf = handles.b.start_frame;
total = handles.b.stop_frame - sf;


% set(handles.hsv_analyzer,'Units','Pixels');
pos = get(handles.hsv_analyzer,'Outerposition');

wbh = waitbar((i-sf)/total,'Please close this to stop playing',...
    'name','Play','Units','Pixels','Outerposition',[pos(1)+pos(3)/2-183 pos(2)-103 366 103]);

fg = figure;
set(fg,'Position',[50 50 1000 740])

subplot(1,2,1)
% set(gca,'nextplot','replacechildren');
axis image
daspect([ 1 1 1])

subplot(1,2,2)
axis image
daspect([ 1 1 1])

colormap(handles.b.color_map_str)

% Crom rectangle
%rect = [165 1 50 70]; % IBP
rect = [40 10 135 230];
%rect = [140 90 50 70]; % Leader
% rect = [164 31 82 61]; % Horizontal leader
%rect = [1 140 30 40];  % Figure 8 IBP

% Image range
range = [5000, 35535];

F =  moviein(handles.stop_frame - handles.start_frame);
b = handles.b;

% Load Arrow data
for i = 1:NAr
    fn = sprintf('%s%s-Arrow%i.txt',b.dir,b.file_name(1:end-5),i);
    fIDA = fopen(fn,'r');
    % Header
    hdr=textscan(fIDA,'%s %s',7,'HeaderLines',2);
    % data
    tmp = textscan(fIDA,'%f %f %f %f','HeaderLines',5);
    fclose(fIDA);
    arrow(i).color      = hdr{2}(1);
    arrow(i).HeadLength = hdr{2}(2);
    arrow(i).HeadWidth  = hdr{2}(3);
    arrow(i).LineWidth  = hdr{2}(4);
    arrow(i).LineLength = hdr{2}(5);
    arrow(i).Dirrection = hdr{2}(6);
    arrow(i).Smooth     = hdr{2}(7);
    
    tmp = cell2mat(tmp);
    tmp = sortrows(tmp,1);
    
    arrow(i).frameN = tmp(:,1);
    arrow(i).N_x    = tmp(:,2) - rect(1)+1;
    arrow(i).N_y = tmp(:,3) - rect(2)+1;
end

i = sf;
% bgim = get_image(handles,sf);
% bgim = imcrop(bgim,rect);

while i < handles.b.stop_frame
    
    matlabIm = get_image(handles,i);
    matlabIm = imcrop(matlabIm,rect);
    
   
    % Show image    
    %cla
    
    %daspect([1 1 1])
    switch handles.b.graphs
        case 1
            %image(double(matlabIm),'CDataMapping','scaled');
            %matlabIm = imresize(matlabIm,5);
            %imagesc(matlabIm,[1697/16128*65535 9219/(16128-2353)*65535])
            imagesc(matlabIm,range)
            %image(double(matlabIm));
            
        case 2
            I = matlabIm;
            matlabIm = imcomplement(matlabIm);
            %image(matlabIm,'CDataMapping','scaled');
            
            subplot(1,2,1)
            %imagesc(matlabIm,range)
            % imagesc(I,[0 10535]) % Fll leader
            imagesc(I,[0 25565])
            %image(I-bgim,'CDataMapping','scaled');
            %image(I,'CDataMapping','scaled');
            time = b.align_time + (i - b.align_frame)/handles.ph.frame_rate;
            title(sprintf('Frame Num: %i Time: %0.6f',i,time))
            
            subplot(1,2,2)
            %I = imresize(matlabIm,10,'lanczos3');
            %I = imresize(I,10,'lanczos3');
            %imagesc(I,range)
            % imagesc(matlabIm,[55000 65535]) %for F11 leader
            imagesc(matlabIm,[40000  65535])
            %image(matlabIm,'CDataMapping','scaled');
            %title('Enhanced')
            
        case 3
            bg = get_image(handles,handles.b.bg_frame);
            matlabIm = double(matlabIm) - double(bg);
            %image(matlabIm,'CDataMapping','scaled');
            imagesc(matlabIm,range)
        case 4
            bg = get_image(handles,handles.b.bg_frame);
            matlabIm = double(matlabIm) - double(bg);
            matlabIm = imcomplement(matlabIm);
            %image(matlabIm,'CDataMapping','scaled');
            imagesc(matlabIm,range)
    end
    
    subplot(1,2,1)
    
  
    % Plotting arrows
    for k = 1:NAr
        ind = sum(arrow(k).frameN <=i);
        
        aFn = arrow(k).frameN(ind);
        aN_x = arrow(k).N_x(ind);
        aN_y = arrow(k).N_y(ind);
        
        try
            aFn2 = arrow(k).frameN(ind+1);
            aN_x2 = arrow(k).N_x(ind+1);
            aN_y2 = arrow(k).N_y(ind+1);
        catch
            aN_x2 = NaN;
            aN_y2 = NaN;
        end
        
        
        
        if str2double(arrow(k).Smooth) && ~isnan(aN_x2)
            aN_x = aN_x + (aN_x2 - aN_x)*(i-aFn)/(aFn2 - aFn);
            aN_y = aN_y + (aN_y2 - aN_y)*(i-aFn)/(aFn2 - aFn);
        end
        
        if ind > 0
            
            switch arrow(k).Dirrection{:}
                case 'right'
                    dx = str2double(arrow(k).LineLength{:});
                    posi = [aN_x-dx, aN_y, dx , 0];
                case 'left'
                    dx = str2double(arrow(k).LineLength{:});
                    posi = [aN_x+dx, aN_y, -dx , 0 ];
                case 'up'
                    dy = str2double(arrow(k).LineLength{:});
                    posi = [aN_x, aN_y+dy, 0, -dy ];
                case 'down'
                    dy = str2double(arrow(k).LineLength{:});
                    posi = [aN_x, aN_y-dy, 0 , dy ];
            end
            
            % [arrow(k).N_x(ind)-dx, arrow(k).N_y(ind)-dy, dx , dy ]

            ah = annotation('arrow');
            set(ah,'parent',gca)
            set(ah,'position',posi,...
                'color',arrow(k).color{:},...
                'HeadLength',str2double(arrow(k).HeadLength{:}),...
                'HeadWidth',str2double(arrow(k).HeadWidth{:}),...
                'LineWidth',str2double(arrow(k).LineWidth{:}))
        end
    end
    
    subplot(1,2,2)
    % Plotting arrows
    for k = 1:NAr
        ind = sum(arrow(k).frameN <=i);
        
        aFn = arrow(k).frameN(ind);
        aN_x = arrow(k).N_x(ind);
        aN_y = arrow(k).N_y(ind);
        
        try
            aFn2 = arrow(k).frameN(ind+1);
            aN_x2 = arrow(k).N_x(ind+1);
            aN_y2 = arrow(k).N_y(ind+1);
        catch
            aN_x2 = NaN;
            aN_y2 = NaN;
        end
        
        
        
        if str2double(arrow(k).Smooth) && ~isnan(aN_x2)
            aN_x = aN_x + (aN_x2 - aN_x)*(i-aFn)/(aFn2 - aFn);
            aN_y = aN_y + (aN_y2 - aN_y)*(i-aFn)/(aFn2 - aFn);
        end
        
        if ind > 0
            
            switch arrow(k).Dirrection{:}
                case 'right'
                    dx = str2double(arrow(k).LineLength{:});
                    posi = [aN_x-dx, aN_y, dx , 0];
                case 'left'
                    dx = str2double(arrow(k).LineLength{:});
                    posi = [aN_x+dx, aN_y, -dx , 0 ];
                case 'up'
                    dy = str2double(arrow(k).LineLength{:});
                    posi = [aN_x, aN_y+dy, 0, -dy ];
                case 'down'
                    dy = str2double(arrow(k).LineLength{:});
                    posi = [aN_x, aN_y-dy, 0 , dy ];
            end
            
            % [arrow(k).N_x(ind)-dx, arrow(k).N_y(ind)-dy, dx , dy ]

            ah = annotation('arrow');
            set(ah,'parent',gca)
            set(ah,'position',posi,...
                'color',arrow(k).color{:},...
                'HeadLength',str2double(arrow(k).HeadLength{:}),...
                'HeadWidth',str2double(arrow(k).HeadWidth{:}),...
                'LineWidth',str2double(arrow(k).LineWidth{:}))
        end
    end
    
    %axis image
    
    subplot(1,2,1)
    % set(gca,'nextplot','replacechildren');
    axis image
    %daspect([ 1 1 1])
    
    subplot(1,2,2)
    axis image
    %daspect([ 1 1 1])
    
    F(i-sf+1) = getframe(gcf);
    
    
    try
        waitbar((i-sf)/total,wbh)
    catch
        %handles.b.current_frame = i;
        %handles = show(handles,i);
        %guidata(hObject,handles)
        return
    end
    
    i = i+1;
    

end
close(fg)
%clf
%M = movie(F);
movie2avi(F,'C:\Users\Sumedhe\Desktop\test3_05fps1.avi',...
    'compression','none',...
    'fps',5, ...
    'quality',100 ...
   )

% handles.b.current_frame = i;
%handles = show(handles,i);
delete(wbh)
guidata(hObject,handles)


% --------------------------------------------------------------------
function create_movie_2_Callback(hObject, eventdata, handles)
clc
i = handles.b.start_frame;
sf = handles.b.start_frame;
total = handles.b.stop_frame - sf

F =  moviein(handles.stop_frame - handles.start_frame);


%fg = figure(handles.b.vh);

set(handles.hsv_analyzer,'Units','Pixels');
pos = get(handles.hsv_analyzer,'Outerposition');

wbh = waitbar((i-sf)/total,'Please close this to stop playing',...
  'name','Play','Units','Pixels','Outerposition',[pos(1)+pos(3)/2-183 pos(2)-103 366 103]);

% try
%     delete(2000)
% end

figure(2000)

handles.b.vh = subplot(2,2,[2 4]);
set(gca,'Color',[0.70 0.70 0.70]);
axis ij
xlim([150 250])
ylim([-160 240])
colormap gray



box on
daspect([1 1 1])
set(gca,'YTick',[-148.7 -97.6 -46.5 4.6 55.7 106.8 157.9 203])
set(gca,'YTickLabel',fliplr([0 1000 2000 3000 4000 5000 6000 7000]))
set(gca,'XTick',[166.9 214.61])
set(gca,'XTickLabel',[3500 4500])
xlabel('Distance (m)')
ylabel('Altitude (m)')




while i < handles.b.stop_frame
    handles = show2(handles,i);
    
    try
        waitbar((i-sf)/total,wbh)
    catch
        handles.b.current_frame = i;
        handles = show2(handles,i);
        return
    end
  
    F(i-sf+1) = getframe(2000);
    i = i+1;
end

movie2avi(F(1:end-1),'C:\Users\Sumedhe\Desktop\Leader2-test.avi',...
    'compression', 'Cinepak',...
    'fps',25, ...
    'quality',20 ...
   )

handles.b.current_frame = i;
handles = show2(handles,i);
delete(wbh)






function handles = show2(handles,N)
% This will display N th image
matlabIm = get_image(handles,N);

subplot(handles.b.vh);


cla
switch handles.b.graphs
    case 1
         %image(double(matlabIm),'CDataMapping','scaled');
         %matlabIm = imresize(matlabIm,5);
         %imagesc(matlabIm,[1697/16128*65535 9219/(16128-2353)*65535])
         imagesc(matlabIm,[0, 65535])
         %image(double(matlabIm));
         
     case 2
        matlabIm = imcomplement(matlabIm);
        image(matlabIm,'CDataMapping','scaled');
        imagesc(matlabIm,[0, 55535])
        %pause(0.1)
    case 3
        bg = get_image(handles,handles.b.bg_frame);
        matlabIm = double(matlabIm) - double(bg);
        %image(matlabIm,'CDataMapping','scaled');
        imagesc(matlabIm)
    case 4
        bg = get_image(handles,handles.b.bg_frame);
        matlabIm = double(matlabIm) - double(bg);
        matlabIm = imcomplement(matlabIm);
        image(matlabIm,'CDataMapping','scaled'); 
        %imagesc(matlabIm,[0, 65535])
end
        
daspect([1 1 1])
title(num2str(N))
hold all


% Calculate middle frame times
f = handles.ph.frame_rate;
n = handles.b.n_of_e_frames;
T = handles.b.align_time;


% Time range for the veiw
t1 = T + (N - handles.b.align_frame - (n+1)/2)/f;
t2 = t1 + n/f;

handles.b.window_t1 = t1;
handles.b.window_t2 = t2;


% Grids for middle frame
t2 = t2-(n-1)/2/f;
t1 = t2-1/f;

handles.b.t1 = t1;
handles.b.t2 = t2;
% 
% % Plotting LDAR/PBFA/LINET
if handles.b.ldar; handles = plot_ldar(handles); end;
if handles.b.pbfa; handles = plot_pbfa(handles); end;
if handles.b.hsvp; handles = plot_hsvp(handles); end;


% 
% % xlim([0 320])
% % ylim([-150 240])
% % zticks = 0:1000:7000;
% % [scrX scrY] = ldar2scr(zeros(size(zticks))-12816,zeros(size(zticks))+1254,zticks,handles)
% % [x y z] = scr2ldar(sscr,scrY,handles);
% % set(gca,'YTickLabel',round(z))
% 
% 
% 
% 
% % E-field graph arrange
% if handles.b.plot_e_field && ishandle(handles.b.gh)
%     handles = set_xrange(handles,N);
% end
% 
% % Update slider
% update_slider(handles,N)


% --------------------------------------------------------------------
function image_integrator_Callback(hObject, eventdata, handles)

% Get image number to intergrate
clc

prompt = {'Iimage numbers:','Image number(s) to keep the backgroud:','Black threshold','Dim'};
dlg_title = 'Image integrator';
num_lines = 1;
def = {'-55136,-52768,-50894,-43133,-38710,-34842,-46309',...
       '-46309',...
       '0.025, 0.035,0.01,0.2,0.2,0.01,0.07', ...
       '0.75,0.75,0.75,0,0.0,0.55,0'};

%
answer = inputdlg(prompt,dlg_title,num_lines,def);

if isempty(answer)
    return
end

images = str2num(answer{1});
keep_bg = str2num(answer{2});
bg_trsld = str2num(answer{3});
dim = str2num(answer{4});

bg = get_image(handles,handles.b.bg_frame);

mm1 = max(bg(:));

im = bg - bg;
% Get images
for i = 1:length(images)
    if ismember(images(i),keep_bg)
        im = im + get_image(handles,images(i));
    else
        temp = get_image(handles,images(i)) - bg;
        temp(temp < bg_trsld(i)*mm1) = 0;
        temp = temp*dim(i);
        %temp = medfilt2(temp);
        temp = wiener2(temp,[5 5]);
        im = im + temp;
        
    end
end

figure
colormap(gray)
image(im,'CDataMapping','scaled');
axis('image')
hold all
if handles.b.ldar; handles = plot_ldar(handles); end
if handles.b.pbfa; handles = plot_pbfa(handles); end;
if handles.b.hsvp; handles = plot_hsvp(handles); end;


% --------------------------------------------------------------------
function plot_pbhsv_on_current_zt_Callback(hObject, eventdata, handles)

d= handles.hsvpd;

b = handles.b;

t = b.align_time + (d(:,1) - b.align_frame)./handles.ph.frame_rate;

[x y z] = scr2ldar(d(:,2),d(:,3),handles);

lg = get(findobj(gcf,'Type','axes','Tag','legend'),'string')


figure(gcf)
h=findobj(gcf,'Type','axes');
lg = get(findobj(gcf,'Type','axes','Tag','legend'),'string');
%[AX,H1,H2]=plotyy(nan,nan,nan,nan);
%plot (AX(2),DLS(:,10),DLS(:,8),'ko','MarkerFaceColor','k','MarkerSize',mz)
plot(h(2),t,z,'go','markerfacecolor','g','markersize',1)
legend([lg 'PBHSV'])


% --------------------------------------------------------------------
function follow_intensity_of_pixels_Callback(hObject, eventdata, handles)

% This function will plot individual pixel intensities (normalized by
% deviding 2^16) with frame numbers. The function can be accesed via
% "Tools" of viedo_analyzer.m program.

%% User inputs

% Pixels to follow in (x,y)
pixels = [190 111
          100 100
          20  20];
%% Program     
% Number of pixels requestd     
[L1, ~] = size(pixels);

% Let's see whther the file name is given
if isempty(handles.b.dir)
    errordlg('Please choose a cine file first','Cine File Error')
    return
end


wbh = waitbar(0,'Please Wait','name','Loking for bright frames');

i_array = handles.b.start_frame:handles.b.stop_frame;
leng = length(i_array);

% variable to save intencity data

int_data = nan(leng,L1);

for i=1:leng
    
    try 
        x = i/leng;
        msg = sprintf('Analysing frame %i of %i',i,leng);
        waitbar(x,wbh,msg)
    catch
        fprintf('Finding Return Strock stopped by user\n')
        return
    end

    matlabIm = get_image(handles,i_array(i));
    inds = sub2ind(size(matlabIm),pixels(:,1),pixels(:,2));
    % Get intensities
    Is = matlabIm(inds);
    int_data(i,:) = Is';
end

% Plot data
figure
hold all
box on
ylabel('Normalized light intensity')
xlabel(sprintf('Frame number from the frame %i',i_array(1)))
for i = 1:L1
    %plot(int_data(:,i)/65536)
    %stairs(i_array,int_data(:,i)/65536)
    stairs(int_data(:,i)/65536)
end

try
    delete(wbh)
end
