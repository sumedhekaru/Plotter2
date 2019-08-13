function varargout = anotation(varargin)
% ANOTATION MATLAB code for anotation.fig
%      ANOTATION, by itself, creates a new ANOTATION or raises the existing
%      singleton*.
%
%      H = ANOTATION returns the handle to a new ANOTATION or the handle to
%      the existing singleton*.
%
%      ANOTATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANOTATION.M with the given input arguments.
%
%      ANOTATION('Property','Value',...) creates a new ANOTATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before anotation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to anotation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help anotation

% Last Modified by GUIDE v2.5 28-Mar-2015 02:25:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @anotation_OpeningFcn, ...
                   'gui_OutputFcn',  @anotation_OutputFcn, ...
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


% --- Executes just before anotation is made visible.
function anotation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to anotation (see VARARGIN)

% Choose default command line output for anotation
handles.output = hObject;


% Default values
fig = figure;
hold all
a.fig = fig;
a.ax = gca;
xL = get(a.ax,'Xlim');
yL = get(a.ax,'Ylim');
a.textstr = 'Enter text';
a.color = 'k';
a.textArrow = 1;
a.objectType = 'line';
a.objectHandle = NaN;
a.solid_line = 1;
a.dtheta = 1;
a.line_ang = 0;
a.dx = 10^(round(log10(range(xL)/10)));
a.dy = 10^(round(log10(range(yL)/10)));
set(handles.dx,'String',a.dx)
set(handles.dy,'String',a.dy)
a.bold = 0;
a.italic = 0;

a.x1 = xL(1);
a.x2 = xL(2);
a.y1 = yL(1) + range(yL)/2;
a.y2 = a.y1;
a.thickness = 1;

a.fontSize = 10;

a.rotate_from = 'middle';
a.x0  = (a.x1+a.x2)/2;
a.y0 = (a.y1 +a.y2)/2;

handles.a = a;
handles = draw_line(handles);


      

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes anotation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = anotation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in left.
function left_Callback(hObject, eventdata, handles)
% hObject    handle to left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.a;

if strcmp(a.objectType,'line')
    handles = deleteObject(handles);
    handles.a.x1 = handles.a.x1 - handles.a.dx;
    handles.a.x2 = handles.a.x2 - handles.a.dx;
    handles = draw_line(handles);
elseif sum(strcmp(a.objectType,{'arrow','textarrow'}))
    handles = deleteObject(handles);
    handles.a.x1 = handles.a.x1 - handles.a.dx;
    handles.a.x2 = handles.a.x2 - handles.a.dx;
    handles = draw_arrow(handles);
elseif sum(strcmp(a.objectType,{'rectangle','ellipse'}))
    a.x1 = a.x1 - a.dx;
    set(a.objectHandle,'Position',[a.x1 a.y1 a.x2 a.y2])
    handles.a = a;
elseif strcmp(a.objectType,'text')
    a.x1 = a.x1 - a.dx;
    set(a.objectHandle,'Position',[a.x1,a.y1])
    handles.a = a;
    
end

guidata(hObject,handles)


% --- Executes on button press in right.
function right_Callback(hObject, eventdata, handles)
% hObject    handle to right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

if strcmp(a.objectType,'line')
    handles = deleteObject(handles);
    handles.a.x1 = handles.a.x1 + handles.a.dx;
    handles.a.x2 = handles.a.x2 + handles.a.dx;
    handles = draw_line(handles);
elseif sum(strcmp(a.objectType,{'arrow','textarrow'}))
    handles = deleteObject(handles);
    handles.a.x1 = handles.a.x1 + handles.a.dx;
    handles.a.x2 = handles.a.x2 + handles.a.dx;
    handles = draw_arrow(handles);
elseif sum(strcmp(a.objectType,{'rectangle','ellipse'}))
    a.x1 = a.x1 + a.dx;
    set(a.objectHandle,'Position',[a.x1 a.y1 a.x2 a.y2])
    handles.a = a;
elseif strcmp(a.objectType,'text')
    a.x1 = a.x1 + a.dx;
    set(a.objectHandle,'Position',[a.x1,a.y1])
    handles.a = a;    
end
        

guidata(hObject,handles)


% --- Executes on button press in up.
function up_Callback(hObject, eventdata, handles)
% hObject    handle to up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.a;

if strcmp(a.objectType,'line')
    handles = deleteObject(handles);
    handles.a.y1 = handles.a.y1 + handles.a.dy;
    handles.a.y2 = handles.a.y2 + handles.a.dy;
    handles = draw_line(handles);
elseif sum(strcmp(a.objectType,{'arrow','textarrow'}))
    handles = deleteObject(handles);
    handles.a.y1 = handles.a.y1 + handles.a.dy;
    handles.a.y2 = handles.a.y2 + handles.a.dy;
    handles = draw_arrow(handles);
elseif sum(strcmp(a.objectType,{'rectangle','ellipse'}))
    a.y1 = a.y1 + a.dy;
    set(a.objectHandle,'Position',[a.x1 a.y1 a.x2 a.y2])
    handles.a = a;
elseif strcmp(a.objectType,'text')
    a.y1 = a.y1 + a.dy;
    set(a.objectHandle,'Position',[a.x1,a.y1])
    handles.a = a;
    
end

guidata(hObject,handles)


% --- Executes on button press in down.
function down_Callback(hObject, eventdata, handles)
% hObject    handle to down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

if strcmp(a.objectType,'line')
    handles = deleteObject(handles);
    handles.a.y1 = handles.a.y1 - handles.a.dy;
    handles.a.y2 = handles.a.y2 - handles.a.dy;
    handles = draw_line(handles);
elseif sum(strcmp(a.objectType,{'arrow','textarrow'}))
    handles = deleteObject(handles);
    handles.a.y1 = handles.a.y1 - handles.a.dy;
    handles.a.y2 = handles.a.y2 - handles.a.dy;
    handles = draw_arrow(handles);
elseif sum(strcmp(a.objectType,{'rectangle','ellipse'}))
    a.y1 = a.y1 - a.dy;
    set(a.objectHandle,'Position',[a.x1 a.y1 a.x2 a.y2])
    handles.a = a;
elseif strcmp(a.objectType,'text')
    a.y1 = a.y1 - a.dy;
    set(a.objectHandle,'Position',[a.x1,a.y1])
    handles.a = a;
    
end

guidata(hObject,handles)

% --- Executes on button press in rotate_ccw.
function rotate_ccw_Callback(hObject, eventdata, handles)
% hObject    handle to rotate_ccw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

if sum(strcmp(a.objectType,{'line','arrow','textarrow'}))
    
%     % Lets convert coordinates to dataspace coordinates.
%     [x0,y0] = nfu2ds(a.ax,a.x0,a.y0);
%     [x1,y1] = nfu2ds(a.ax,a.x1,a.y1);
%     [x2,y2] = nfu2ds(a.ax,a.x2,a.y2);
%     
%     [THETA,R] = cart2pol(x1-x0,y1-y0); %Convert to polar coordinates
%     THETA=THETA+a.dtheta*pi/180; %Add a_rad to theta
%     [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
%     x1 = xr + x0;
%     y1 = yr + y0;
%     
%     
%     [THETA,R] = cart2pol(x2-x0,y2-y0); %Convert to polar coordinates
%     THETA=THETA+a.dtheta*pi/180; %Add a_rad to theta
%     [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
%     x2 = xr + x0;
%     y2 = yr + y0;
%     
%     % Convert back to figure coordinates from data space coordinates
%     [a.x1, a.y1] = ds2nfu(a.ax,x1,y1);
%     [a.x2, a.y2] = ds2nfu(a.ax,x2,y2);

    
    [THETA,R] = cart2pol(a.x1-a.x0,a.y1-a.y0); %Convert to polar coordinates
    THETA=THETA+a.dtheta*pi/180; %Add a_rad to theta
    [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
    a.x1 = xr + a.x0;
    a.y1 = yr + a.y0;
    
    [THETA,R] = cart2pol(a.x2-a.x0,a.y2-a.y0); %Convert to polar coordinates
    THETA=THETA+a.dtheta*pi/180; %Add a_rad to theta
    [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
    a.x2 = xr + a.x0;
    a.y2 = yr + a.y0;

    fprintf('Angle = %0.1f degrees\n',THETA*180/pi)
    a.line_ang = THETA*180/pi;
    handles.a = a;
    handles = deleteObject(handles);
    if strcmp(a.objectType,'line')
        handles = draw_line(handles);
    else
        handles = draw_arrow(handles);
    end
    
elseif strcmp(a.objectType,'text')
       a.line_ang = a.line_ang + a.dtheta;
       set(a.objectHandle,'Rotation',a.line_ang)  
       handles.a = a;
       
end


guidata(hObject,handles)



    



% --- Executes on button press in rotate_cw.
function rotate_cw_Callback(hObject, eventdata, handles)
% hObject    handle to rotate_cw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

if sum(strcmp(a.objectType,{'line','arrow','textarrow'}))
        
        % Convert to data space units
%         [x
%         [xn,yn] = nfu2ds(hAx,x,y)
       
    
    
        [THETA,R] = cart2pol(a.x1-a.x0,a.y1-a.y0); %Convert to polar coordinates
        THETA=THETA-a.dtheta*pi/180; %Add a_rad to theta
        [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
        a.x1 = xr + a.x0;
        a.y1 = yr + a.y0;
        
        [THETA,R] = cart2pol(a.x2-a.x0,a.y2-a.y0); %Convert to polar coordinates
        THETA=THETA-a.dtheta*pi/180; %Add a_rad to theta
        [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
        a.x2 = xr + a.x0;
        a.y2 = yr + a.y0;
        
        fprintf('Angle = %0.1f degrees\n',THETA*180/pi)
        a.line_ang = THETA*180/pi;
        
        handles.a = a;
        handles = deleteObject(handles);
        if strcmp(a.objectType,'line')
            handles = draw_line(handles);
        else
            handles = draw_arrow(handles);
        end
elseif strcmp(a.objectType, 'text')
        a.line_ang = a.line_ang - a.dtheta;        
        set(a.objectHandle,'Rotation',a.line_ang)  
        handles.a = a;
        
end
        
guidata(hObject,handles)

function dtheta_Callback(hObject, eventdata, handles)
% hObject    handle to dtheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dtheta as text
%        str2double(get(hObject,'String')) returns contents of dtheta as a double
handles.a.dtheta = str2double(get(hObject,'String'));

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function dtheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dtheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dy_Callback(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dy as text
%        str2double(get(hObject,'String')) returns contents of dy as a double
handles.a.dy = str2double(get(hObject,'String'));

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function dy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx_Callback(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx as text
%        str2double(get(hObject,'String')) returns contents of dx as a double

handles.a.dx = str2double(get(hObject,'String'));

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function dx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in newline_h.
function newline_h_Callback(hObject, eventdata, handles)
% hObject    handle to newline_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles. a;
xL = get(a.ax,'Xlim');
yL = get(a.ax,'Ylim');

a.x1 = xL(1);
a.x2 = xL(2);
a.y1 = yL(1) + range(yL)/2;
a.y2 = a.y1;

a.line_ang = 0;

handles.a = a;
handles = get_pivot(handles);
handles = draw_line(handles);
guidata(hObject,handles)



% --- Executes on button press in newline_v.
function newline_v_Callback(hObject, eventdata, handles)
% hObject    handle to newline_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles. a;
xL = get(a.ax,'Xlim');
yL = get(a.ax,'Ylim');

a.x1 = xL(1) + range(xL)/2;
a.x2 = a.x1;
a.y1 = yL(1);
a.y2 = yL(2);

a.line_ang = 90;

handles.a = a;
handles = get_pivot(handles);
handles = draw_line(handles);
guidata(hObject,handles)


% --- Executes on button press in newArrow.
function newArrow_Callback(hObject, eventdata, handles)
% hObject    handle to newArrow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.a;
Lh = imline;
pos = Lh.getPosition;
delete(Lh);
a.x1 = pos(1);
a.x2 = pos(2);
a.y1 = pos(3);
a.y2 = pos(4);

a.line_ang = cart2pol(a.x2-a.x1,a.y2-a.y1); %Convert to polar coordinates

handles.a = a;
handles = get_pivot(handles);
handles = draw_arrow(handles);
guidata(hObject,handles)


function h = draw_arrow(h)
a = h.a;

[x1,y1] =  ds2nfu(a.ax,a.x1,a.y1);
[x2,y2] =  ds2nfu(a.ax,a.x2,a.y2);

% Looking for a text arrow
if h.a.textArrow
   a.objectHandle = annotation('textarrow',[x1 x2],[y1 y2],...
        'String',a.textstr);
    a.objectType = 'textarrow';
    
    if a.bold
        set(a.objectHandle,'FontWeight','bold')
    else
        set(a.objectHandle,'FontWeight','normal')
    end
    
    if a.italic
        set(a.objectHandle,'FontAngle','italic')
    else
        set(a.objectHandle,'FontAngle','normal')
    end
else
    a.objectHandle = annotation('arrow',[x1 x2],[y1 y2]);
    a.objectType = 'arrow';
end

set(a.objectHandle,'LineWidth',a.thickness,'color',a.color);

if a.solid_line
    set(a.objectHandle,'LineStyle','-')
else
    set(a.objectHandle,'LineStyle','--')
end



% Pin to axis
hA = handle(a.objectHandle);
hA.pinAtAffordance(1);
hA.pinAtAffordance(2);

h.a = a;





% --- Executes on button press in textArrow.
function textArrow_Callback(hObject, eventdata, handles)
% hObject    handle to textArrow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of textArrow
handles.a.textArrow = get(hObject,'Value');
if sum(strcmp(handles.a.objectType,{'arrow','textarrow'}))
    handles = deleteObject(handles);
    handles = draw_arrow(handles);
end


guidata(hObject,handles)



function textstr_Callback(hObject, eventdata, handles)
% hObject    handle to textstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.a;
a.textstr = get(hObject,'String');

if strcmp(a.objectType,'text') || strcmp(a.objectType,'textarrow')
    set(a.objectHandle,'String',a.textstr)
end

handles.a = a;

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function textstr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in newtext.
function newtext_Callback(hObject, eventdata, handles)
% hObject    handle to newtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

pos = ginput(1);
%pos = get(a.objectHandle,'Position');
a.x1 = pos(1);
a.y1 = pos(2);
a.x2 = NaN; 
a.y2 = NaN;
a.objectType = 'text';


a.line_ang = 0;
handles.a = a;
handles = add_new_text(handles);
guidata(hObject,handles)


function handles = add_new_text(handles)

a = handles.a;
a.objectType = 'text';
a.objectHandle = text(a.x1,a.y1,a.textstr);
set(a.objectHandle,'color',a.color,'FontSize',a.fontSize,'Rotation',a.line_ang);
%get(a.objectHandle)

if a.bold
    set(a.objectHandle,'FontWeight','bold')
else
    set(a.objectHandle,'FontWeight','normal')
end

if a.italic
    set(a.objectHandle,'FontAngle','italic')
else
    set(a.objectHandle,'FontAngle','normal')
end


handles.a = a;




% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)
% hObject    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
deleteObject(handles);

function h = draw_line(h)

a = h.a;
a.objectType = 'line';
if a.solid_line
    lineSpec = a.color;
else
    lineSpec = [a.color '--'];
end

a.objectHandle = plot(a.ax,[a.x1,a.x2],[a.y1,a.y2],lineSpec,...
    'LineWidth',a.thickness);
a.objectType = 'line';
h.a = a;



function h = deleteObject(h)

try
delete(h.a.objectHandle)

catch
   disp('Object delete is not found')
end


% --- Executes during object creation, after setting all properties.
function rotate_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotate_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes during object creation, after setting all properties.
function rotate_le_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotate_le (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function rotate_panel_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to rotate_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function rotate_panel_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to rotate_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in rotate_panel.
function rotate_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in rotate_panel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

handles = get_pivot(handles);
guidata(hObject,handles)

function handles = get_pivot(handles)

val = get(handles.rotate_le,'Value');
if val; 
    handles.a.rotate_from = 'leftend'; 
    handles.a.x0  = handles.a.x1;
    handles.a.y0 =  handles.a.y1;
end

val = get(handles.rotate_mid,'Value');
if val; 
    handles.a.rotate_from = 'middle'; 
    handles.a.x0  = (handles.a.x1+handles.a.x2)/2;
    handles.a.y0 = (handles.a.y1 +handles.a.y2)/2;
end

val = get(handles.rotate_re,'Value');
if val; 
    handles.a.rotate_from ='rightend'; 
    handles.a.x0  = handles.a.x2;
    handles.a.y0 = handles.a.y2;
end





% --- Executes when selected object is changed in line_specs.
function line_specs_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in line_specs 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

handles.a.solid_line = get(handles.solid,'Value');

if sum(strcmp(handles.a.objectType,...
        {'line','arrow','textarrow','rectangle','ellipse'}))
    if handles.a.solid_line
        set(handles.a.objectHandle,'LineStyle','-')
    else
        set(handles.a.objectHandle,'LineStyle','--')
    end
end


guidata(hObject,handles)


% --- Executes when selected object is changed in color_panel.
function color_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in color_panel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

val = get(handles.red,'Value');
if val; handles.a.color = 'r'; end

val = get(handles.green,'Value');
if val; handles.a.color = 'g'; end

val = get(handles.blue,'Value');
if val; handles.a.color = 'b'; end

val = get(handles.black,'Value');
if val; handles.a.color = 'k'; end


if sum(strcmp(handles.a.objectType,{'line','arrow','textarrow'}))
    set(handles.a.objectHandle,'Color',handles.a.color)
elseif sum(strcmp(handles.a.objectType,{'rectangle','ellipse'}))
    set(handles.a.objectHandle,'EdgeColor',handles.a.color)
end
    
guidata(hObject,handles)




function thickness_Callback(hObject, eventdata, handles)
% hObject    handle to thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thickness as text
%        str2double(get(hObject,'String')) returns contents of thickness as a double

handles.a.thickness = str2double(get(hObject,'String'));

if strcmp(handles.a.objectType,'text')
        deleteObject(handles);
        handles = add_new_text(handles);
else
        set(handles.a.objectHandle,'LineWidth',handles.a.thickness)
end


guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in shorten.
function shorten_Callback(hObject, eventdata, handles)
% hObject    handle to shorten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.a;

if sum(strcmp(a.objectType,{'line','arrow','textarrow'}))

    if a.line_ang < 45 || a.line_ang > 315 || (a.line_ang > 135 && a.line_ang < 225)
        newy = interp1([a.x1 a.x2],[a.y1 a.y2],[a.x1 + a.dx/2, a.x2 - a.dx/2],'line','extrap');
        a.x1 = a.x1 + a.dx/2;
        a.x2 = a.x2 - a.dx/2;
        a.y1 = newy(1);
        a.y2 = newy(2);
    else
        newx = interp1([a.y1 a.y2],[a.x1 a.x2],[a.y1 + a.dy/2, a.y2 - a.dy/2],'line','extrap');
        a.y1 = a.y1 + a.dy/2;
        a.y2 = a.y2 - a.dy/2;
        a.x1 = newx(1);
        a.x2 = newx(2);
    end
    
    
    handles.a = a;
    handles = deleteObject(handles);
    
    if strcmp(a.objectType,'line')
        handles = draw_line(handles);
    else
        handles = draw_arrow(handles);
    end
  
end

guidata(hObject,handles)


% --- Executes on button press in lengthen.
function lengthen_Callback(hObject, eventdata, handles)
% hObject    handle to lengthen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

if sum(strcmp(a.objectType,{'line','arrow','textarrow'}))

    if a.line_ang < 45 || a.line_ang > 315 || (a.line_ang > 135 && a.line_ang < 225)
        newy = interp1([a.x1 a.x2],[a.y1 a.y2],[a.x1 - a.dx/2, a.x2 + a.dx/2],'line','extrap');
        a.x1 = a.x1 - a.dx/2;
        a.x2 = a.x2 + a.dx/2;
        a.y1 = newy(1);
        a.y2 = newy(2);
    else
        newx = interp1([a.y1 a.y2],[a.x1 a.x2],[a.y1 - a.dy/2, a.y2 + a.dy/2],'line','extrap');
        a.y1 = a.y1 - a.dy/2;
        a.y2 = a.y2 + a.dy/2;
        a.x1 = newx(1);
        a.x2 = newx(2);
    end
    
    
    handles.a = a;
    handles = deleteObject(handles);
    
    if strcmp(a.objectType,'line')
        handles = draw_line(handles);
    else
        handles = draw_arrow(handles);
    end
  
end
guidata(hObject,handles)



function dL_Callback(hObject, eventdata, handles)
% hObject    handle to dL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dL as text
handles.a.dL = str2double(get(hObject,'String')); % returns contents of dL as a double
add_data_to_figure(a,'update')
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function dL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hand_draw.
function hand_draw_Callback(hObject, eventdata, handles)
% hObject    handle to hand_draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;
Lh = imline;
pos = Lh.getPosition;
delete(Lh);
a.x1 = pos(1);
a.x2 = pos(2);
a.y1 = pos(3);
a.y2 = pos(4);

a.line_ang = cart2pol(a.x2-a.x1,a.y2-a.y1); %Convert to polar coordinates


handles.a = a;
handles = get_pivot(handles);
handles = draw_line(handles);
guidata(hObject,handles)


% --- Executes on button press in txt_update.
function txt_update_Callback(hObject, eventdata, handles)
% hObject    handle to txt_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




function fontSize_Callback(hObject, eventdata, handles)
% hObject    handle to fontSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;
a.fontSize = str2double(get(hObject,'String'));
handles.a = a;

if strcmp(a.objectType,'text')
    handles = deleteObject(handles);
    handles = add_new_text(handles);
elseif strcmp(a.objectType,'textarrow')
    set(handles.a.objectHandle,'FontSize',a.fontSize)
end

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function fontSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fontSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bold.
function bold_Callback(hObject, eventdata, handles)
% hObject    handle to bold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.a.bold = get(hObject,'Value');

if sum(strcmp(handles.a.objectType,{'text','textarrow'}))
    if handles.a.bold
        set(handles.a.objectHandle,'FontWeight','bold')
    else
        set(handles.a.objectHandle,'FontWeight','normal')
    end    
end


guidata(hObject,handles)



% --- Executes on button press in italic.
function italic_Callback(hObject, eventdata, handles)
% hObject    handle to italic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of italic

handles.a.italic = get(hObject,'Value');

if sum(strcmp(handles.a.objectType,{'text','textarrow'}))
    if handles.a.italic
        set(handles.a.objectHandle,'FontAngle','italic')
    else
        set(handles.a.objectHandle,'FontAngle','normal')
    end    
end

guidata(hObject,handles)


function [xn,yn] = nfu2ds(hAx,x,y)

%% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');
axpos = get(hAx,'Position');
axlim = axis(hAx);
axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));


%% Transform data

%xn = (x-axlim(1))*axpos(3)/axwidth + axpos(1);
xn = axlim(1) + (x-axpos(1))*axwidth/axpos(3);
%yn = (y-axlim(3))*axpos(4)/axheight + axpos(2);
yn = axlim(3) + (y-axpos(2))*axheight/axpos(4);

%% Restore axes units
set(hAx,'Units',axun)



   


% --- Executes on button press in ellipse.
function ellipse_Callback(hObject, eventdata, handles)
% hObject    handle to ellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
a = handles.a;

msg = sprintf('(1) Please draw the eliipse on the figure. \n(2) Do the changes you want. \n(3) Double click in the ellipse to finish.');

th = annotation('textbox',[0.15 0.75 0.6 0.15],'string',msg);

h = imellipse(a.ax); 
pos = wait(h);
pos = getPosition(h);
% Delete massage and imrect
delete(h)
delete(th)

if isempty(pos)
    return
end

a.x1 = pos(1);
a.y1 = pos(2);
a.x2 = pos(3);
a.y2 = pos(4);

handles.a = a;
handles = draw_ellipse(handles);
guidata(hObject,handles)

function handles = draw_ellipse(handles)

a = handles.a;
a.objectType = 'ellipse';
a.objectHandle = rectangle('Position',[a.x1,a.y1,a.x2,a.y2],...
    'LineWidth',a.thickness,'edgecolor',a.color,'curvature',[1,1]);
   
if a.solid_line
    set(a.objectHandle,'LineStyle','-')
else
    set(a.objectHandle,'LineStyle','--')
end

handles.a = a;


% --- Executes on button press in rectangle.
function rectangle_Callback(hObject, eventdata, handles)
% hObject    handle to rectangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

msg = sprintf('(1) Please draw the rectangle on the figure. \n(2) Do the changes you want. \n(3) Double click in the rectangle to finish.');

th = annotation('textbox',[0.15 0.75 0.6 0.15],'string',msg);

h = imrect(a.ax); 
pos = wait(h);
% Delete massage and imrect
delete(h)
delete(th)

if isempty(pos)
    return
end

a.x1 = pos(1);
a.y1 = pos(2);
a.x2 = pos(3);
a.y2 = pos(4);

handles.a = a;
handles = draw_rect(handles);
guidata(hObject,handles)

function handles = draw_rect(handles)
a = handles.a;
a.objectType = 'rectangle';
a.objectHandle = rectangle('Position',[a.x1,a.y1,a.x2,a.y2],...
    'LineWidth',a.thickness,'edgecolor',a.color);
   
if a.solid_line
    set(a.objectHandle,'LineStyle','-')
else
    set(a.objectHandle,'LineStyle','--')
end

handles.a = a;

function varargout = anotation(varargin)
% ANOTATION MATLAB code for anotation.fig
%      ANOTATION, by itself, creates a new ANOTATION or raises the existing
%      singleton*.
%
%      H = ANOTATION returns the handle to a new ANOTATION or the handle to
%      the existing singleton*.
%
%      ANOTATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANOTATION.M with the given input arguments.
%
%      ANOTATION('Property','Value',...) creates a new ANOTATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before anotation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to anotation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help anotation

% Last Modified by GUIDE v2.5 28-Mar-2015 02:25:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @anotation_OpeningFcn, ...
                   'gui_OutputFcn',  @anotation_OutputFcn, ...
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


% --- Executes just before anotation is made visible.
function anotation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to anotation (see VARARGIN)

% Choose default command line output for anotation
handles.output = hObject;


% Default values
fig = figure;
hold all
a.fig = fig;
a.ax = gca;
xL = get(a.ax,'Xlim');
yL = get(a.ax,'Ylim');
a.textstr = 'Enter text';
a.color = 'k';
a.textArrow = 1;
a.objectType = 'line';
a.objectHandle = NaN;
a.solid_line = 1;
a.dtheta = 1;
a.line_ang = 0;
a.dx = 10^(round(log10(range(xL)/10)));
a.dy = 10^(round(log10(range(yL)/10)));
set(handles.dx,'String',a.dx)
set(handles.dy,'String',a.dy)
a.bold = 0;
a.italic = 0;

a.x1 = xL(1);
a.x2 = xL(2);
a.y1 = yL(1) + range(yL)/2;
a.y2 = a.y1;
a.thickness = 1;

a.fontSize = 10;

a.rotate_from = 'middle';
a.x0  = (a.x1+a.x2)/2;
a.y0 = (a.y1 +a.y2)/2;

handles.a = a;
handles = draw_line(handles);


      

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes anotation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = anotation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in left.
function left_Callback(hObject, eventdata, handles)
% hObject    handle to left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.a;

if strcmp(a.objectType,'line')
    handles = deleteObject(handles);
    handles.a.x1 = handles.a.x1 - handles.a.dx;
    handles.a.x2 = handles.a.x2 - handles.a.dx;
    handles = draw_line(handles);
elseif sum(strcmp(a.objectType,{'arrow','textarrow'}))
    handles = deleteObject(handles);
    handles.a.x1 = handles.a.x1 - handles.a.dx;
    handles.a.x2 = handles.a.x2 - handles.a.dx;
    handles = draw_arrow(handles);
elseif sum(strcmp(a.objectType,{'rectangle','ellipse'}))
    a.x1 = a.x1 - a.dx;
    set(a.objectHandle,'Position',[a.x1 a.y1 a.x2 a.y2])
    handles.a = a;
elseif strcmp(a.objectType,'text')
    a.x1 = a.x1 - a.dx;
    set(a.objectHandle,'Position',[a.x1,a.y1])
    handles.a = a;
    
end

guidata(hObject,handles)


% --- Executes on button press in right.
function right_Callback(hObject, eventdata, handles)
% hObject    handle to right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

if strcmp(a.objectType,'line')
    handles = deleteObject(handles);
    handles.a.x1 = handles.a.x1 + handles.a.dx;
    handles.a.x2 = handles.a.x2 + handles.a.dx;
    handles = draw_line(handles);
elseif sum(strcmp(a.objectType,{'arrow','textarrow'}))
    handles = deleteObject(handles);
    handles.a.x1 = handles.a.x1 + handles.a.dx;
    handles.a.x2 = handles.a.x2 + handles.a.dx;
    handles = draw_arrow(handles);
elseif sum(strcmp(a.objectType,{'rectangle','ellipse'}))
    a.x1 = a.x1 + a.dx;
    set(a.objectHandle,'Position',[a.x1 a.y1 a.x2 a.y2])
    handles.a = a;
elseif strcmp(a.objectType,'text')
    a.x1 = a.x1 + a.dx;
    set(a.objectHandle,'Position',[a.x1,a.y1])
    handles.a = a;    
end
        

guidata(hObject,handles)


% --- Executes on button press in up.
function up_Callback(hObject, eventdata, handles)
% hObject    handle to up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.a;

if strcmp(a.objectType,'line')
    handles = deleteObject(handles);
    handles.a.y1 = handles.a.y1 + handles.a.dy;
    handles.a.y2 = handles.a.y2 + handles.a.dy;
    handles = draw_line(handles);
elseif sum(strcmp(a.objectType,{'arrow','textarrow'}))
    handles = deleteObject(handles);
    handles.a.y1 = handles.a.y1 + handles.a.dy;
    handles.a.y2 = handles.a.y2 + handles.a.dy;
    handles = draw_arrow(handles);
elseif sum(strcmp(a.objectType,{'rectangle','ellipse'}))
    a.y1 = a.y1 + a.dy;
    set(a.objectHandle,'Position',[a.x1 a.y1 a.x2 a.y2])
    handles.a = a;
elseif strcmp(a.objectType,'text')
    a.y1 = a.y1 + a.dy;
    set(a.objectHandle,'Position',[a.x1,a.y1])
    handles.a = a;
    
end

guidata(hObject,handles)


% --- Executes on button press in down.
function down_Callback(hObject, eventdata, handles)
% hObject    handle to down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

if strcmp(a.objectType,'line')
    handles = deleteObject(handles);
    handles.a.y1 = handles.a.y1 - handles.a.dy;
    handles.a.y2 = handles.a.y2 - handles.a.dy;
    handles = draw_line(handles);
elseif sum(strcmp(a.objectType,{'arrow','textarrow'}))
    handles = deleteObject(handles);
    handles.a.y1 = handles.a.y1 - handles.a.dy;
    handles.a.y2 = handles.a.y2 - handles.a.dy;
    handles = draw_arrow(handles);
elseif sum(strcmp(a.objectType,{'rectangle','ellipse'}))
    a.y1 = a.y1 - a.dy;
    set(a.objectHandle,'Position',[a.x1 a.y1 a.x2 a.y2])
    handles.a = a;
elseif strcmp(a.objectType,'text')
    a.y1 = a.y1 - a.dy;
    set(a.objectHandle,'Position',[a.x1,a.y1])
    handles.a = a;
    
end

guidata(hObject,handles)

% --- Executes on button press in rotate_ccw.
function rotate_ccw_Callback(hObject, eventdata, handles)
% hObject    handle to rotate_ccw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

if sum(strcmp(a.objectType,{'line','arrow','textarrow'}))
    
%     % Lets convert coordinates to dataspace coordinates.
%     [x0,y0] = nfu2ds(a.ax,a.x0,a.y0);
%     [x1,y1] = nfu2ds(a.ax,a.x1,a.y1);
%     [x2,y2] = nfu2ds(a.ax,a.x2,a.y2);
%     
%     [THETA,R] = cart2pol(x1-x0,y1-y0); %Convert to polar coordinates
%     THETA=THETA+a.dtheta*pi/180; %Add a_rad to theta
%     [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
%     x1 = xr + x0;
%     y1 = yr + y0;
%     
%     
%     [THETA,R] = cart2pol(x2-x0,y2-y0); %Convert to polar coordinates
%     THETA=THETA+a.dtheta*pi/180; %Add a_rad to theta
%     [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
%     x2 = xr + x0;
%     y2 = yr + y0;
%     
%     % Convert back to figure coordinates from data space coordinates
%     [a.x1, a.y1] = ds2nfu(a.ax,x1,y1);
%     [a.x2, a.y2] = ds2nfu(a.ax,x2,y2);

    
    [THETA,R] = cart2pol(a.x1-a.x0,a.y1-a.y0); %Convert to polar coordinates
    THETA=THETA+a.dtheta*pi/180; %Add a_rad to theta
    [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
    a.x1 = xr + a.x0;
    a.y1 = yr + a.y0;
    
    [THETA,R] = cart2pol(a.x2-a.x0,a.y2-a.y0); %Convert to polar coordinates
    THETA=THETA+a.dtheta*pi/180; %Add a_rad to theta
    [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
    a.x2 = xr + a.x0;
    a.y2 = yr + a.y0;

    fprintf('Angle = %0.1f degrees\n',THETA*180/pi)
    a.line_ang = THETA*180/pi;
    handles.a = a;
    handles = deleteObject(handles);
    if strcmp(a.objectType,'line')
        handles = draw_line(handles);
    else
        handles = draw_arrow(handles);
    end
    
elseif strcmp(a.objectType,'text')
       a.line_ang = a.line_ang + a.dtheta;
       set(a.objectHandle,'Rotation',a.line_ang)  
       handles.a = a;
       
end


guidata(hObject,handles)



    



% --- Executes on button press in rotate_cw.
function rotate_cw_Callback(hObject, eventdata, handles)
% hObject    handle to rotate_cw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

if sum(strcmp(a.objectType,{'line','arrow','textarrow'}))
        
        % Convert to data space units
%         [x
%         [xn,yn] = nfu2ds(hAx,x,y)
       
    
    
        [THETA,R] = cart2pol(a.x1-a.x0,a.y1-a.y0); %Convert to polar coordinates
        THETA=THETA-a.dtheta*pi/180; %Add a_rad to theta
        [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
        a.x1 = xr + a.x0;
        a.y1 = yr + a.y0;
        
        [THETA,R] = cart2pol(a.x2-a.x0,a.y2-a.y0); %Convert to polar coordinates
        THETA=THETA-a.dtheta*pi/180; %Add a_rad to theta
        [xr,yr] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
        a.x2 = xr + a.x0;
        a.y2 = yr + a.y0;
        
        fprintf('Angle = %0.1f degrees\n',THETA*180/pi)
        a.line_ang = THETA*180/pi;
        
        handles.a = a;
        handles = deleteObject(handles);
        if strcmp(a.objectType,'line')
            handles = draw_line(handles);
        else
            handles = draw_arrow(handles);
        end
elseif strcmp(a.objectType, 'text')
        a.line_ang = a.line_ang - a.dtheta;        
        set(a.objectHandle,'Rotation',a.line_ang)  
        handles.a = a;
        
end
        
guidata(hObject,handles)

function dtheta_Callback(hObject, eventdata, handles)
% hObject    handle to dtheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dtheta as text
%        str2double(get(hObject,'String')) returns contents of dtheta as a double
handles.a.dtheta = str2double(get(hObject,'String'));

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function dtheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dtheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dy_Callback(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dy as text
%        str2double(get(hObject,'String')) returns contents of dy as a double
handles.a.dy = str2double(get(hObject,'String'));

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function dy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx_Callback(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx as text
%        str2double(get(hObject,'String')) returns contents of dx as a double

handles.a.dx = str2double(get(hObject,'String'));

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function dx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in newline_h.
function newline_h_Callback(hObject, eventdata, handles)
% hObject    handle to newline_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles. a;
xL = get(a.ax,'Xlim');
yL = get(a.ax,'Ylim');

a.x1 = xL(1);
a.x2 = xL(2);
a.y1 = yL(1) + range(yL)/2;
a.y2 = a.y1;

a.line_ang = 0;

handles.a = a;
handles = get_pivot(handles);
handles = draw_line(handles);
guidata(hObject,handles)



% --- Executes on button press in newline_v.
function newline_v_Callback(hObject, eventdata, handles)
% hObject    handle to newline_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles. a;
xL = get(a.ax,'Xlim');
yL = get(a.ax,'Ylim');

a.x1 = xL(1) + range(xL)/2;
a.x2 = a.x1;
a.y1 = yL(1);
a.y2 = yL(2);

a.line_ang = 90;

handles.a = a;
handles = get_pivot(handles);
handles = draw_line(handles);
guidata(hObject,handles)


% --- Executes on button press in newArrow.
function newArrow_Callback(hObject, eventdata, handles)
% hObject    handle to newArrow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.a;
Lh = imline;
pos = Lh.getPosition;
delete(Lh);
a.x1 = pos(1);
a.x2 = pos(2);
a.y1 = pos(3);
a.y2 = pos(4);

a.line_ang = cart2pol(a.x2-a.x1,a.y2-a.y1); %Convert to polar coordinates

handles.a = a;
handles = get_pivot(handles);
handles = draw_arrow(handles);
guidata(hObject,handles)


function h = draw_arrow(h)
a = h.a;

[x1,y1] =  ds2nfu(a.ax,a.x1,a.y1);
[x2,y2] =  ds2nfu(a.ax,a.x2,a.y2);

% Looking for a text arrow
if h.a.textArrow
   a.objectHandle = annotation('textarrow',[x1 x2],[y1 y2],...
        'String',a.textstr);
    a.objectType = 'textarrow';
    
    if a.bold
        set(a.objectHandle,'FontWeight','bold')
    else
        set(a.objectHandle,'FontWeight','normal')
    end
    
    if a.italic
        set(a.objectHandle,'FontAngle','italic')
    else
        set(a.objectHandle,'FontAngle','normal')
    end
else
    a.objectHandle = annotation('arrow',[x1 x2],[y1 y2]);
    a.objectType = 'arrow';
end

set(a.objectHandle,'LineWidth',a.thickness,'color',a.color);

if a.solid_line
    set(a.objectHandle,'LineStyle','-')
else
    set(a.objectHandle,'LineStyle','--')
end



% Pin to axis
hA = handle(a.objectHandle);
hA.pinAtAffordance(1);
hA.pinAtAffordance(2);

h.a = a;





% --- Executes on button press in textArrow.
function textArrow_Callback(hObject, eventdata, handles)
% hObject    handle to textArrow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of textArrow
handles.a.textArrow = get(hObject,'Value');
if sum(strcmp(handles.a.objectType,{'arrow','textarrow'}))
    handles = deleteObject(handles);
    handles = draw_arrow(handles);
end


guidata(hObject,handles)



function textstr_Callback(hObject, eventdata, handles)
% hObject    handle to textstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.a;
a.textstr = get(hObject,'String');

if strcmp(a.objectType,'text') || strcmp(a.objectType,'textarrow')
    set(a.objectHandle,'String',a.textstr)
end

handles.a = a;

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function textstr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in newtext.
function newtext_Callback(hObject, eventdata, handles)
% hObject    handle to newtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

pos = ginput(1);
%pos = get(a.objectHandle,'Position');
a.x1 = pos(1);
a.y1 = pos(2);
a.x2 = NaN; 
a.y2 = NaN;
a.objectType = 'text';


a.line_ang = 0;
handles.a = a;
handles = add_new_text(handles);
guidata(hObject,handles)


function handles = add_new_text(handles)

a = handles.a;
a.objectType = 'text';
a.objectHandle = text(a.x1,a.y1,a.textstr);
set(a.objectHandle,'color',a.color,'FontSize',a.fontSize,'Rotation',a.line_ang);
%get(a.objectHandle)

if a.bold
    set(a.objectHandle,'FontWeight','bold')
else
    set(a.objectHandle,'FontWeight','normal')
end

if a.italic
    set(a.objectHandle,'FontAngle','italic')
else
    set(a.objectHandle,'FontAngle','normal')
end


handles.a = a;




% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)
% hObject    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
deleteObject(handles);

function h = draw_line(h)

a = h.a;
a.objectType = 'line';
if a.solid_line
    lineSpec = a.color;
else
    lineSpec = [a.color '--'];
end

a.objectHandle = plot(a.ax,[a.x1,a.x2],[a.y1,a.y2],lineSpec,...
    'LineWidth',a.thickness);
a.objectType = 'line';
h.a = a;



function h = deleteObject(h)

try
delete(h.a.objectHandle)

catch
   disp('Object delete is not found')
end


% --- Executes during object creation, after setting all properties.
function rotate_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotate_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes during object creation, after setting all properties.
function rotate_le_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotate_le (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function rotate_panel_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to rotate_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function rotate_panel_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to rotate_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in rotate_panel.
function rotate_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in rotate_panel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

handles = get_pivot(handles);
guidata(hObject,handles)

function handles = get_pivot(handles)

val = get(handles.rotate_le,'Value');
if val; 
    handles.a.rotate_from = 'leftend'; 
    handles.a.x0  = handles.a.x1;
    handles.a.y0 =  handles.a.y1;
end

val = get(handles.rotate_mid,'Value');
if val; 
    handles.a.rotate_from = 'middle'; 
    handles.a.x0  = (handles.a.x1+handles.a.x2)/2;
    handles.a.y0 = (handles.a.y1 +handles.a.y2)/2;
end

val = get(handles.rotate_re,'Value');
if val; 
    handles.a.rotate_from ='rightend'; 
    handles.a.x0  = handles.a.x2;
    handles.a.y0 = handles.a.y2;
end





% --- Executes when selected object is changed in line_specs.
function line_specs_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in line_specs 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

handles.a.solid_line = get(handles.solid,'Value');

if sum(strcmp(handles.a.objectType,...
        {'line','arrow','textarrow','rectangle','ellipse'}))
    if handles.a.solid_line
        set(handles.a.objectHandle,'LineStyle','-')
    else
        set(handles.a.objectHandle,'LineStyle','--')
    end
end


guidata(hObject,handles)


% --- Executes when selected object is changed in color_panel.
function color_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in color_panel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

val = get(handles.red,'Value');
if val; handles.a.color = 'r'; end

val = get(handles.green,'Value');
if val; handles.a.color = 'g'; end

val = get(handles.blue,'Value');
if val; handles.a.color = 'b'; end

val = get(handles.black,'Value');
if val; handles.a.color = 'k'; end


if sum(strcmp(handles.a.objectType,{'line','arrow','textarrow'}))
    set(handles.a.objectHandle,'Color',handles.a.color)
elseif sum(strcmp(handles.a.objectType,{'rectangle','ellipse'}))
    set(handles.a.objectHandle,'EdgeColor',handles.a.color)
end
    
guidata(hObject,handles)




function thickness_Callback(hObject, eventdata, handles)
% hObject    handle to thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thickness as text
%        str2double(get(hObject,'String')) returns contents of thickness as a double

handles.a.thickness = str2double(get(hObject,'String'));

if strcmp(handles.a.objectType,'text')
        deleteObject(handles);
        handles = add_new_text(handles);
else
        set(handles.a.objectHandle,'LineWidth',handles.a.thickness)
end


guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in shorten.
function shorten_Callback(hObject, eventdata, handles)
% hObject    handle to shorten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.a;

if sum(strcmp(a.objectType,{'line','arrow','textarrow'}))

    if a.line_ang < 45 || a.line_ang > 315 || (a.line_ang > 135 && a.line_ang < 225)
        newy = interp1([a.x1 a.x2],[a.y1 a.y2],[a.x1 + a.dx/2, a.x2 - a.dx/2],'line','extrap');
        a.x1 = a.x1 + a.dx/2;
        a.x2 = a.x2 - a.dx/2;
        a.y1 = newy(1);
        a.y2 = newy(2);
    else
        newx = interp1([a.y1 a.y2],[a.x1 a.x2],[a.y1 + a.dy/2, a.y2 - a.dy/2],'line','extrap');
        a.y1 = a.y1 + a.dy/2;
        a.y2 = a.y2 - a.dy/2;
        a.x1 = newx(1);
        a.x2 = newx(2);
    end
    
    
    handles.a = a;
    handles = deleteObject(handles);
    
    if strcmp(a.objectType,'line')
        handles = draw_line(handles);
    else
        handles = draw_arrow(handles);
    end
  
end

guidata(hObject,handles)


% --- Executes on button press in lengthen.
function lengthen_Callback(hObject, eventdata, handles)
% hObject    handle to lengthen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

if sum(strcmp(a.objectType,{'line','arrow','textarrow'}))

    if a.line_ang < 45 || a.line_ang > 315 || (a.line_ang > 135 && a.line_ang < 225)
        newy = interp1([a.x1 a.x2],[a.y1 a.y2],[a.x1 - a.dx/2, a.x2 + a.dx/2],'line','extrap');
        a.x1 = a.x1 - a.dx/2;
        a.x2 = a.x2 + a.dx/2;
        a.y1 = newy(1);
        a.y2 = newy(2);
    else
        newx = interp1([a.y1 a.y2],[a.x1 a.x2],[a.y1 - a.dy/2, a.y2 + a.dy/2],'line','extrap');
        a.y1 = a.y1 - a.dy/2;
        a.y2 = a.y2 + a.dy/2;
        a.x1 = newx(1);
        a.x2 = newx(2);
    end
    
    
    handles.a = a;
    handles = deleteObject(handles);
    
    if strcmp(a.objectType,'line')
        handles = draw_line(handles);
    else
        handles = draw_arrow(handles);
    end
  
end
guidata(hObject,handles)



function dL_Callback(hObject, eventdata, handles)
% hObject    handle to dL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dL as text
handles.a.dL = str2double(get(hObject,'String')); % returns contents of dL as a double
add_data_to_figure(a,'update')
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function dL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hand_draw.
function hand_draw_Callback(hObject, eventdata, handles)
% hObject    handle to hand_draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;
Lh = imline;
pos = Lh.getPosition;
delete(Lh);
a.x1 = pos(1);
a.x2 = pos(2);
a.y1 = pos(3);
a.y2 = pos(4);

a.line_ang = cart2pol(a.x2-a.x1,a.y2-a.y1); %Convert to polar coordinates


handles.a = a;
handles = get_pivot(handles);
handles = draw_line(handles);
guidata(hObject,handles)


% --- Executes on button press in txt_update.
function txt_update_Callback(hObject, eventdata, handles)
% hObject    handle to txt_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




function fontSize_Callback(hObject, eventdata, handles)
% hObject    handle to fontSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;
a.fontSize = str2double(get(hObject,'String'));
handles.a = a;

if strcmp(a.objectType,'text')
    handles = deleteObject(handles);
    handles = add_new_text(handles);
elseif strcmp(a.objectType,'textarrow')
    set(handles.a.objectHandle,'FontSize',a.fontSize)
end

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function fontSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fontSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bold.
function bold_Callback(hObject, eventdata, handles)
% hObject    handle to bold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.a.bold = get(hObject,'Value');

if sum(strcmp(handles.a.objectType,{'text','textarrow'}))
    if handles.a.bold
        set(handles.a.objectHandle,'FontWeight','bold')
    else
        set(handles.a.objectHandle,'FontWeight','normal')
    end    
end


guidata(hObject,handles)



% --- Executes on button press in italic.
function italic_Callback(hObject, eventdata, handles)
% hObject    handle to italic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of italic

handles.a.italic = get(hObject,'Value');

if sum(strcmp(handles.a.objectType,{'text','textarrow'}))
    if handles.a.italic
        set(handles.a.objectHandle,'FontAngle','italic')
    else
        set(handles.a.objectHandle,'FontAngle','normal')
    end    
end

guidata(hObject,handles)


function [xn,yn] = nfu2ds(hAx,x,y)

%% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');
axpos = get(hAx,'Position');
axlim = axis(hAx);
axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));


%% Transform data

%xn = (x-axlim(1))*axpos(3)/axwidth + axpos(1);
xn = axlim(1) + (x-axpos(1))*axwidth/axpos(3);
%yn = (y-axlim(3))*axpos(4)/axheight + axpos(2);
yn = axlim(3) + (y-axpos(2))*axheight/axpos(4);

%% Restore axes units
set(hAx,'Units',axun)



   


% --- Executes on button press in ellipse.
function ellipse_Callback(hObject, eventdata, handles)
% hObject    handle to ellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
a = handles.a;

msg = sprintf('(1) Please draw the eliipse on the figure. \n(2) Do the changes you want. \n(3) Double click in the ellipse to finish.');

th = annotation('textbox',[0.15 0.75 0.6 0.15],'string',msg);

h = imellipse(a.ax); 
pos = wait(h);
pos = getPosition(h);
% Delete massage and imrect
delete(h)
delete(th)

if isempty(pos)
    return
end

a.x1 = pos(1);
a.y1 = pos(2);
a.x2 = pos(3);
a.y2 = pos(4);

handles.a = a;
handles = draw_ellipse(handles);
guidata(hObject,handles)

function handles = draw_ellipse(handles)

a = handles.a;
a.objectType = 'ellipse';
a.objectHandle = rectangle('Position',[a.x1,a.y1,a.x2,a.y2],...
    'LineWidth',a.thickness,'edgecolor',a.color,'curvature',[1,1]);
   
if a.solid_line
    set(a.objectHandle,'LineStyle','-')
else
    set(a.objectHandle,'LineStyle','--')
end

handles.a = a;


% --- Executes on button press in rectangle.
function rectangle_Callback(hObject, eventdata, handles)
% hObject    handle to rectangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = handles.a;

msg = sprintf('(1) Please draw the rectangle on the figure. \n(2) Do the changes you want. \n(3) Double click in the rectangle to finish.');

th = annotation('textbox',[0.15 0.75 0.6 0.15],'string',msg);

h = imrect(a.ax); 
pos = wait(h);
% Delete massage and imrect
delete(h)
delete(th)

if isempty(pos)
    return
end

a.x1 = pos(1);
a.y1 = pos(2);
a.x2 = pos(3);
a.y2 = pos(4);

handles.a = a;
handles = draw_rect(handles);
guidata(hObject,handles)

function handles = draw_rect(handles)
a = handles.a;
a.objectType = 'rectangle';
a.objectHandle = rectangle('Position',[a.x1,a.y1,a.x2,a.y2],...
    'LineWidth',a.thickness,'edgecolor',a.color);
   
if a.solid_line
    set(a.objectHandle,'LineStyle','-')
else
    set(a.objectHandle,'LineStyle','--')
end

handles.a = a;

