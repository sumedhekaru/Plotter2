function extendticklabel(varargin)
%EXTENDTICKLABEL extends the number of significant figures that an axes' label displays
%  EXTENDTICKLABEL(H,AX,POINT) takes the axes with handle H, and extends the 
%  labels of AX to POINT significant figures.  AX should be a string with 
%  value 'x', 'y', or 'z'.  POINT should be an integer.
%
%  EXTENDTICKLABEL(AX,POINT) will assume H is the current axes, GCA.
%
%  The motivation behind this tool is the fact that MATLAB defaults to 
%  five decimal places when creating axes ticklabels.
%
%  Example:
%       x = 1:10;  y = 5000:0.000001:5000+9*0.000001;
%       plot(x,y)
%       extendticklabel(gca,'y',10)

%  Greg Aloe
%  Revision 0.2
%  Copyright 2002 The Mathworks, Inc.

% Parse inputs
error(nargchk(2,3,nargin))

point = varargin{end};

if floor(point) ~= point
    point = floor(point);
    warning(sprintf('Non-integer input POINT was rounded down to %i.',point))
end

ax    = varargin{end-1};
ax    = ax(1);
if ~strcmp(ax,'x') & ~strcmp(ax,'y') & ~strcmp(ax,'z')
    error('Input AX must take the char value ''x'', ''y'', or ''z''.')
end

if nargin == 3
    h = varargin{1};
else
    h = gca;
end

curformat = get(0,'format');
curvals = get(h,[ax 'tick']);

format long
newticklabels = num2str(curvals',point);
set(0,'format',curformat)
set(h,[ax 'ticklabel'],newticklabels);

% Manipulate axes position to fit the new labels
shift = point - 6;
if (shift > 0) & strcmp(ax,'y')
    curaxunits = get(h,'units');
    set(h,'units','char')
    newaxpos = get(h,'position');
    newaxpos(1) = newaxpos(1)+shift;
    newaxpos(3) = newaxpos(3)-shift;
    set(h,'position',newaxpos);
    set(h,'units',curaxunits);
end

if strcmp(get(get(h,'parent'),'renderer'),'OpenGL')
    warning(sprintf(['If an exponent is provided inadvertently for the ' ...
        'ticklabel,\n         it can be removed by changing the figure''s' ...
        ' renderer to ''zbuffer''.\n         ' ...
        'This is an OpenGL bug in MATLAB 6.1 (R12.1) and previous.']))
end