function florida_map

x_limit = xlim;
y_limit = ylim;

%% Early method (load csv file)
% Calculating x and y
% ldar_lat=28.538486111;
% ldar_lon=80.642633333;
% earth_R=6371e3;

% 
% x = csvread('fl_map_data.csv') ;
% 
% 
% % x
% %data{6}=111200*cos((ldar_lat+data{7})./2.*pi./360).*(ldar_lon-(-data{6}));
% x0 = 111319.491*cos((ldar_lat).*2.*pi./360).*(ldar_lon-(88.0831));
% % y
% %data{7}=111000*(-ldar_lat+data{7});
% y0 = 111319.491*(-ldar_lat+23.9787);
% 
% % x1 = 111319.491*cos((ldar_lat).*2.*pi./360).*(ldar_lon-(80.80838084220886));
% % y1 = 111319.491*(-ldar_lat+28.62142660839247);
% % 
% % x2 = 111319.491*cos((ldar_lat).*2.*pi./360).*(ldar_lon-(80.1552951335907));
% % y2 = 111319.491*(-ldar_lat+27.167840646558165);
% % data.x = (x(:,2)+x0)/1000;
% % data.y = (x(:,1)+y0)/1000;
% % save('florida_map.mat','-Struct','data')
% hold all
% hLine = plot((x(:,2)+x0)/1000,(x(:,1)+y0)/1000,'Color',[0.6 0.6 0.6]);
% set(get(get(hLine,'Annotation'),'LegendInformation'),...
%                         'IconDisplayStyle','off')
% %hold all
% % plot(0,0,'ro')
% % plot(x1,y1,'go')
% % plot(x2,y2,'ko')
% %pbaspect([1 1 1])
% axis([x_limit y_limit])
% 
% %disp('test')

%% New method (load a mat file)
sen_set = open('sensor_setting.mat');

switch sen_set.map_quality
    case 1; data = load('florida_map2.mat');
    case 2; data = load('florida_map3.mat');
    case 3; data = load('florida_map4.mat');
end

hold all
hLine = plot(data.x,data.y,'Color',[0.6 0.6 0.6]);
set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                    
axis([x_limit y_limit])