function test_save_fig(figh)

if nargin < 2
    figh = gcf;
end

% get all the axis handles in the figure
hAllAxes = findobj(figh,'type','axes');
hLeg = findobj(hAllAxes,'tag','legend');
hAxes = setdiff(hAllAxes,hLeg);

% get position of the figure
figProp.fig_position = get(figh,'position');

figProp.hAxes = hAxes;
% get locations of all the axes
for i = 1:length(hAxes);
    figProp.position(i).position = get(hAxes(i),'position');
end

% is there current data in this figure

data = guidata(figh)

% add figure properties to the figure
data.figProp = figProp;
guidata(figh,data)

% save the figure
saveas(figh,'C:\Users\Sumedhe\Desktop\testfig.fig')
disp('done')
    