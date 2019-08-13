function varargout = figure_light_balance(varargin)
% FIGURE_LIGHT_BALANCE MATLAB code for figure_light_balance.fig
%      FIGURE_LIGHT_BALANCE, by itself, creates a new FIGURE_LIGHT_BALANCE or raises the existing
%      singleton*.
%
%      H = FIGURE_LIGHT_BALANCE returns the handle to a new FIGURE_LIGHT_BALANCE or the handle to
%      the existing singleton*.
%
%      FIGURE_LIGHT_BALANCE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIGURE_LIGHT_BALANCE.M with the given input arguments.
%
%      FIGURE_LIGHT_BALANCE('Property','Value',...) creates a new FIGURE_LIGHT_BALANCE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before figure_light_balance_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to figure_light_balance_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help figure_light_balance

% Last Modified by GUIDE v2.5 30-Jun-2012 16:10:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @figure_light_balance_OpeningFcn, ...
                   'gui_OutputFcn',  @figure_light_balance_OutputFcn, ...
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


% --- Executes just before figure_light_balance is made visible.
function figure_light_balance_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to figure_light_balance (see VARARGIN)
clc
handles.mainH = varargin{1};
handles.figure = 7500;
handles.high = 65535;
handles.low = 0;
handles.I = get_image(handles.mainH,handles.mainH.b.current_frame);
handles.rect = [1 1 240 320];
handles.tit = num2str(handles.mainH.b.current_frame);
handles.pbfad.x = [];
handles.pbfad.y = [];
handles.ldard.x = [];
handles.ldard.y = [];

set(handles.frames,'string',num2str(handles.mainH.b.current_frame));
handles.exp_val = 1;
show(handles);

% Choose default command line output for figure_light_balance
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes figure_light_balance wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = figure_light_balance_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function low_val_Callback(hObject, eventdata, handles)
handles.low = get(hObject,'Value');
set(handles.bottomtext,'String',num2str(handles.low));
guidata(hObject, handles);
show(handles)


% --- Executes during object creation, after setting all properties.
function low_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to low_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function high_val_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.high = get(hObject,'Value');
set(handles.toptext,'String',num2str(handles.high));
guidata(hObject, handles);
show(handles)


% --- Executes during object creation, after setting all properties.
function high_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to high_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
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


function show(handles)
if handles.high < handles.low
    disp('Top slider must have a higher value')
    return
end


figure(handles.figure);
clf
%colormap(handles.mainH.b.color_map_str)
colormap(gray(2^8))
I = handles.I;
%I = wiener2(I,[2 2]);
%I = medfilt2(I);
%I = ordfilt2(I,4,true(2));
%I = imresize(I,handles.exp_val,'lanczos3');

% Getting negative
%I = imcomplement(I);

imagesc(I,[handles.low, handles.high])

daspect([1 1 1])
hold all

% plot(handles.pbfad.x,handles.pbfad.y,'ro',...
%     'markerfacecolor','r','markersize',2)
% 
% plot(handles.ldard.x,handles.ldard.y,'go',...
%     'markerfacecolor','g','markersize',2)

title(handles.tit)



function frames_Callback(hObject, eventdata, handles)
clc
handles.tit = get(hObject,'String');
images = str2num(handles.tit);
lineIm = [];
rect = handles.rect - [0 0 1 1];

sumIm = [];
backgro = get_image(handles.mainH,handles.mainH.b.bg_frame);
backgro = imcrop(backgro,rect);
avgval = mean2(backgro);



handles.pbfad.x = [];
handles.pbfad.y = [];
handles.ldard.x = [];
handles.ldard.y = [];

for i = 1:length(images)
    I = get_image(handles.mainH,images(i));
    I = imcrop(I,rect);
    lineIm = [lineIm I ];
    handles = get_pbfa(handles,images(i),i);
    handles = get_ldar(handles,images(i),i);
   
    %handles = get_ldar(handles);
    if i == 1;
        sumIm = I;
        mean2(I-backgro)
    elseif images(i) == -31619
        newI = I - backgro;
        newI(find(newI < avgval*0.2)) = 0;
        sumIm = sumIm + newI;     
    else
        newI = I - backgro;
        newI(find(newI < avgval*0.08)) = 0;
        sumIm = sumIm + newI;
     
    end  
    
end

% figure
% imagesc(sumIm,[handles.low, handles.high])
% daspect([1 1 1])
% colormap(gray(2^8));
% return

%x = handles.pbfad.x
%y = handles.pbfad.y
handles.I = lineIm;

%handles.I = lineIm;
guidata(hObject, handles);
show(handles);



% --- Executes during object creation, after setting all properties.
function frames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function topleftx_Callback(hObject, eventdata, handles)
handles.rect(1) = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function topleftx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to topleftx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in crop.
function crop_Callback(hObject, eventdata, handles)
fg = figure(handles.figure);
[temp temp tem rect] = imcrop();
handles.rect = rect;
set(handles.topleftx,'String',rect(1))
set(handles.toplefty,'String',rect(2))
set(handles.xrange,'String',rect(3))
set(handles.yrange,'String',rect(4))
guidata(hObject, handles);
show(handles)




function toplefty_Callback(hObject, eventdata, handles)
handles.rect(2) = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function toplefty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to toplefty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xrange_Callback(hObject, eventdata, handles)
handles.rect(3) = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function xrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yrange_Callback(hObject, eventdata, handles)
handles.rect(4) = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function yrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function toptext_Callback(hObject, eventdata, handles)

handles.high =   str2double(get(hObject,'String'));
set(handles.high_val,'value',handles.high)
guidata(hObject, handles);
show(handles)




% --- Executes during object creation, after setting all properties.
function toptext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to toptext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bottomtext_Callback(hObject, eventdata, handles)
handles.low =   str2double(get(hObject,'String'));
set(handles.low_val,'value',handles.low)
guidata(hObject, handles);
show(handles)


% --- Executes during object creation, after setting all properties.
function bottomtext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bottomtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function expand_slider_Callback(hObject, eventdata, handles)
handles.exp_val = get(hObject,'Value');
set(handles.expand_val,'String',num2str(handles.exp_val));
guidata(hObject, handles);
show(handles)



% --- Executes during object creation, after setting all properties.
function expand_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expand_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function expand_val_Callback(hObject, eventdata, handles)
handles.exp_val =   str2double(get(hObject,'String'));
set(handles.low_val,'value',handles.exp_val)
guidata(hObject, handles);
show(handles)

% --- Executes during object creation, after setting all properties.
function expand_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expand_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function orgH = get_ldar(handles,image,i)
% This function will plot all the location points (eventhough it named as
% PLOT_LDAR

orgH = handles;
rect = handles.rect;
handles = handles.mainH;

% If ldar data not available, just return
try
    temp = handles.ldard;
catch
    return
end;

b = handles.b;

% Let's load LDAR data if has not been loaded yet
% try 
%     temp = handles.ldard(:,1);
% catch
%     disp('LDAR data has not loaded yet')
%     return
% end

% Calculate middle frame times
f = handles.ph.frame_rate;
%n = handles.b.n_of_e_frames;
T = handles.b.align_time;

% Grids for middle frame
t1 = T + (image - handles.b.align_frame - 2)/f;
t2 = t1 + 3/f;



% LDAR points in the current frame
lol = sum(handles.ldard(:,1) < t1)+1;
ul  = sum(handles.ldard(:,1) <= t2);

% Plot Current ldar point
if b.ldar
    orgH.ldard.x = [orgH.ldard.x; handles.ldard(lol:ul,5)-rect(1)+(i-1)*rect(3)];
    orgH.ldard.y = [orgH.ldard.y; handles.ldard(lol:ul,6)-rect(2)];
       
%    plot(handles.ldard(lol:ul,5)-rect(1),handles.ldard(lol:ul,6)--rect(2),'go','MarkerFaceColor','g','MarkerSize',2)
end

% % Plot early ldar points
% if b.pre_plot
%     plot(handles.ldard(1:lol-1,5)-rect(1),handles.ldard(1:lol-1,6)-rect(2),'go','MarkerFaceColor','g','MarkerSize',2)
% end
% 
% % Plot post ldar points
% if b.post_plot
%     plot(handles.ldard(ul+1:end,5)-rect(1),handles.ldard(ul+1:end,6)-rect(2),'go','MarkerFaceColor','g','MarkerSize',2)
% end



function orgH = get_pbfa(handles,image,i)
% This function will plot all the location points (eventhough it named as
% PLOT_LDAR

% original handle
orgH = handles;
rect = handles.rect;
handles = handles.mainH;


% If ldar data not available, just return
try
    temp = handles.pbfad;
catch
    return
end;

b = handles.b;

% Let's load LDAR data if has not been loaded yet
%try 
%    temp = handles.pbfad(:,1);
%catch
%    handles = load_ldar(handles);
%end

% Calculate middle frame times
f = handles.ph.frame_rate;
%n = handles.b.n_of_e_frames;
T = handles.b.align_time;

% Grids for middle frame
t1 = T + (image - handles.b.align_frame - 2)/f;
t2 = t1 + 3/f;

% LDAR points in the current frame
lol = sum(handles.pbfad(:,1) < t1)+1;
ul  = sum(handles.pbfad(:,1) <= t2);

% Plot Current ldar point
if b.pbfa
    orgH.pbfad.x = [orgH.pbfad.x; handles.pbfad(lol:ul,5)-rect(1)+(i-1)*rect(3)];
    orgH.pbfad.y = [orgH.pbfad.y; handles.pbfad(lol:ul,6)-rect(2)];
    
    %plot(handles.pbfad(lol:ul,5)-rect(1),handles.pbfad(lol:ul,6)-rect(2),'ro','MarkerFaceColor','r','MarkerSize',2)
end

%x = orgH.pbfad.x
%y = orgH.pbfad.y


% % Plot early ldar points
% if b.pre_plot
%     plot(handles.pbfad(1:lol-1,5)-rect(1),handles.pbfad(1:lol-1,6)-rect(2),'ro','MarkerFaceColor','r','MarkerSize',2)
% end
% 
% % Plot post ldar points
% if b.post_plot
%     plot(handles.pbfad(ul+1:end,5)-rect(1),handles.pbfad(ul+1:end,6)-rect(2),'ro','MarkerFaceColor','r','MarkerSize',2)
% end
