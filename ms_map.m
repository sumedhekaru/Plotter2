function ms_map
% Draw the ms county map

x_limit = xlim;
y_limit = ylim;

data = load('ms_map.mat');
hold all
hLine = plot(data.x/1000,data.y/1000,'Color',[0.6 0.6 0.6]);
set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                    
axis([x_limit y_limit])