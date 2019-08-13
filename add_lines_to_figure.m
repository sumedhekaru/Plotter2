function add_lines_to_figure(figH)

if nargin < 1
    figure
    figH = gcf;
    axeH = gca;
else
    figH = figH;
    axeH = gca;
end

hold all

% Create the GUI figure
fig = figure('MenuBar','none');
bg = get(fig,'color');
pos = get(figH,'position');
L = 240;
H = 180;
set(fig,'position',[pos(1)+pos(3)+20, pos(2)+pos(4)/2-H/2,L,H])



% Create the handles structure
h = guihandles(fig);
h.figH = figH;
h.axeH = axeH;

% Add horizontal line to the figure
xL = get(h.axeH,'Xlim');
yL = get(h.axeH,'Ylim');
h.x1 = xL(1);
h.x2 = xL(2);
h.y1 = yL(1) + range(yL)/2;
h.y2 = h.y1;

h = plot_line(h);

% Line type
h.solid = 1;


% Add two buttons
h.bt_up = uicontrol(fig,'Style','pushbutton','String','Up',...
                'Position',[80 10 70 20]);

h.bt_down = uicontrol(fig,'Style','pushbutton','String','Down',...
                'Position',[10 10 70 20]);
            
h.bt_left = uicontrol(fig,'Style','pushbutton','String','Left',...
                'Position',[10 40 70 20]);
            
h.bt_right = uicontrol(fig,'Style','pushbutton','String','Right',...
                'Position',[80 40 70 20]);

h.bt_addnewH = uicontrol(fig,'Style','pushbutton','String','New H Line',...
                'Position',[10 120 70 20]); 

h.bt_addnewV = uicontrol(fig,'Style','pushbutton','String','New V Line',...
                'Position',[80 120 70 20]); 
            
h.bt_delete = uicontrol(fig,'Style','pushbutton','String','Delete',...
                'Position',[10 150 70 20]); 
            
h.bt_addArr = uicontrol(fig,'Style','pushbutton','String','New Arrow',...
                'Position',[80 150 70 20]); 
            
% Add textboxes
h.tx_vshift = uicontrol(fig,'Style','edit',...
                'String',10^(round(log10(range(xL)/10))),...
                'Position',[160 10 70 20]);
            
h.tx_hshift = uicontrol(fig,'Style','edit',...
                'String',10^(round(log10(range(yL)/10))),...
                'Position',[160 40 70 20]);

% Add button group for line selection
h.butg = uibuttongroup(fig,'Title','Line type',...
            'units','pixels',...
            'Position',[10 70 220 40],'background',bg);
        
h.line_solid = uicontrol(h.butg,'Style','radiobutton',...
                'String','Solid','background',bg,...
                'Units','normalized',...
                'Value',1,'Position',[0.1 0.5 0.25 0.5]);

h.line_broken = uicontrol(h.butg,'Style','radiobutton',...
                'String','Broken','background',bg,...
                'Units','normalized',...
                'Value',0,'Position',[0.6 0.5 0.25 0.5]);
 

            
            
% add static massage
% h.t1 = uicontrol(fig,'Style','text',...
%                 'String','Vertical E-shift:',...
%                 'Position',[10 40 130 20],'backgroundcolor',bg);



% assign functions
set(h.bt_down,'callback',@moveDown)
set(h.bt_up,'callback',@moveUp)
set(h.bt_left,'callback',@moveLeft)
set(h.bt_right,'callback',@moveRight)
set(h.bt_addnewH,'callback',@newLineH)
set(h.bt_addnewV,'callback',@newLineV)
set(h.bt_addArr,'callback',@newArrow)
set(h.bt_delete,'callback',@delete_curr)
set(h.line_solid,'callback',@test)
set(h.line_broken,'callback',@test)


guidata(fig,h)

function test(hObject,eventdata)
h = guidata(gcbo);
h.solid = get(h.line_solid,'Value');
delete_curr(hObject,eventdata)
h = plot_line(h);
guidata(hObject,h)


function h = moveDown(hObject,eventdata)

h = guidata(gcbo);


% get the vertical change
h.vshift = str2double(get(h.tx_vshift, 'String'));

h.y1 = h.y1 - h.vshift;
h.y2 = h.y2 - h.vshift;

% delete current line
delete_curr(hObject,eventdata)

% draw a new line
h = plot_line(h);

guidata(hObject,h)

function h = moveUp(hObject,eventdata)

h = guidata(gcbo);

% get the vertical change
h.vshift = str2double(get(h.tx_vshift, 'String'));

h.y1 = h.y1 + h.vshift;
h.y2 = h.y2 + h.vshift;

% delete current line
delete_curr(hObject,eventdata)

% draw a new line
h = plot_line(h);


guidata(hObject,h)

function h = moveLeft(hObject,eventdata)

h = guidata(gcbo);

% get the vertical change
h.vshift = str2double(get(h.tx_vshift, 'String'));

h.x1 = h.x1 - h.vshift;
h.x2 = h.x2 - h.vshift;

% delete current line
delete_curr(hObject,eventdata)

% draw a new line
h = plot_line(h);


guidata(hObject,h)

function h = moveRight(hObject,eventdata)

h = guidata(gcbo);

% get the vertical change
h.vshift = str2double(get(h.tx_vshift, 'String'));

h.x1 = h.x1 + h.vshift;
h.x2 = h.x2 + h.vshift;

% delete current line
delete_curr(hObject,eventdata)

% draw a new line
h = plot_line(h);


guidata(hObject,h)


function h = newLineH(hObject,eventdata)

h = guidata(gcbo);

% Creat a new line
xL = get(h.axeH,'Xlim');
yL = get(h.axeH,'Ylim');

h.x1 = xL(1);
h.x2 = xL(2);
h.y1 = yL(1) + range(yL)/2;
h.y2 = h.y1;

h = plot_line(h);

guidata(hObject,h)


function h = newLineV(hObject,eventdata)

h = guidata(gcbo);

% Creat a new line
xL = get(h.axeH,'Xlim');
yL = get(h.axeH,'Ylim');

h.x1 = xL(1)+range(xL)/2;
h.x2 = h.x1;
h.y1 = yL(1);
h.y2 = yL(2);

h = plot_line(h);

guidata(hObject,h)


function h = newArrow(hObject,eventdata)

h = guidata(gcbo);

guidata(hObject,h)



function delete_curr(hObject,eventdata)
h = guidata(gcbo);
% delete current line
try
    delete(h.hLine)
end
guidata(hObject,h)

function h = plot_line(h)

if h.solid
    Lt = '';
else
    Lt = '--'; 
end

h.hLine = plot(h.axeH,[h.x1 h.x2],[h.y1 h.y2],Lt);

