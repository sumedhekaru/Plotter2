function set_figure_size(figh)

% Figur handle
if nargin < 1
    figh = gcf;
end

% get guidata
d = guidata(figh);

try
     % figure size
    fz = d.figProp.fig_position;
catch
    disp('No previous figure properties found')
    return
end
    

% screen size
sz = get(0,'ScreenSize');

% scale the figure if it is longer
if fz(3) > sz(3)
    ratio = sz(3)/fz(3);
    fz(3) = round(fz(3)*ratio*0.9);
    fz(4) = round(fz(4)*ratio*0.9);
end

% scale the figure if it is toller
if fz(4) > sz(4)
    ratio = sz(4)/fz(4);
    fz(3) = round(fz(3)*ratio*0.9);
    fz(4) = round(fz(4)*ratio*0.9);
end

% move the figure to the left if it is out of the area
if (fz(1)+fz(3)) > sz(3)
    fz(1) = 1;
    fz(2) = 1;
end

if (fz(2)+fz(4)) > sz(4)
    fz(2) = 1;
    fz(1) = 1;
end

% set figure size
set(figh,'position',fz)

