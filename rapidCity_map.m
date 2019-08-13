function rapidCity_map

x_limit = xlim;
y_limit = ylim;

d = open('rapidCityMap2.mat');


hold all
hLine = plot(d.x/1000,d.y/1000,'Color',[0.6 0.6 0.6]);
set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
%hold all
% plot(0,0,'ro')
% plot(x1,y1,'go')
% plot(x2,y2,'ko')
%pbaspect([1 1 1])
axis([x_limit y_limit])
%disp('test')