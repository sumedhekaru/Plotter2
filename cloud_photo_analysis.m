function cloud_photo_analysis
% This function is written to analyze cloud photo obtained by Dr. Marshall
% in Aug 2014

%% User inputs
%Photo file name
a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\IMG_2609.jpeg';
%a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\photos_2015_06_04\IMG_3849.jpg';

% Camera info
a.f = 6e-3;  %Focal length
a.p = 1.8e-6; %Fixel Size

% Ground level photo reference (in pixels)
a.scr_x = 479;
a.scr_y = 1181;

% Real location of the (above) reference (in meters)
a.ref_x = -15840;
a.ref_y = -6925;

% Location of camera in meters
a.cam_x = -12110;
a.cam_y = 11603;
a.cam_z = 0;

% Time info
a.yyyy = 2011;
a.mm = 7;
a.dd = 22;
a.t1 = 61558.0;
a.t2 = 61558.2;


% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\IMG_0465.jpg';
% 
% % Camera info
% a.f = 35e-3*1.62;  %Focal length
% a.p = 1.8e-6; %Fixel Size
% 
% % Ground level photo reference (in pixels)
% a.scr_x = 479;
% a.scr_y = 1181;
% 
% % Real location of the (above) reference (in meters)
% a.ref_x = -15840;
% a.ref_y = -6925;
% 
% % Location of camera in meters
% a.cam_x = -12110;
% a.cam_y = 11603;
% a.cam_z = 0;
% 
% % Time info
% a.yyyy = 2011;
% a.mm = 7;
% a.dd = 22;
% a.t1 = 61558.0;
% a.t2 = 61558.2;



%% Start the program
% Load the image
%info = imfinfo(a.pfn);
I = imread(a.pfn);


% Plot the image
figure
imagesc(I)
daspect([1 1 1])
hold all

% Plot ground level
plot(xlim,[a.scr_y a.scr_y],'w')
text(100,a.scr_y+25,'Ground Level','color','w')

% Plot other levels
levels = 1000:1000:12000;
D = sqrt((a.ref_x - a.cam_x).^2 + (a.ref_y - a.cam_y).^2)

% Location of camera in meters
a.cam_x = -12110;
a.cam_y = 11603;

n = (levels*a.f)/(a.p*D);

for i = 1:length(n)
    plot(xlim,[a.scr_y a.scr_y]-n(i),'k')
    text(100,a.scr_y+25-n(i),sprintf('%i km',levels(i)/1000));
end
    
tools2fig


%% Plot LDAR
% Get a copy of plotter 2 data
try
    h=guidata(findall(0,'Tag','plotter2'));
catch
    disp('Run plotter2 first')
    return
end

% Setup timing
g = h.g;
g.YYYY = {sprintf('%4.4i',a.yyyy)};
g.mm = {sprintf('%4.4i',a.mm)};
g.dd = {sprintf('%4.4i',a.dd)};
g.t1 = a.t1;
g.t2 = a.t2;
g.hh = floor(a.t1/3600);
g.mm = floor((a.t1 - g.hh*3600)/60/5)*5;

% Turn off sensors and unwanted stuffs
g.lpgraphs = zeros(1,60);
g.chgraphs = zeros(1,60);
g.linet = 0; g.ldar = 1; g.cglss = 1; g.pbfa = 0; g.nldn = 0; g.pbfa_old = 0;

% Setting file
sen_set = open('sensor_setting.mat');
xlim('auto')
ylim('auto')

%% Get data
data = load_data(g,sen_set);
DLS = data.DLS;
% DLS info
% 1 - occuring time
% 6 - x
% 7 - y
% 8 - z
% 9 - Occuring time
% 10 - Detection time
% 11 - horizontal Distance to each pulse
x = DLS(:,6);
y = DLS(:,7);
z = DLS(:,8);
t = DLS(:,1);

[N_x,N_y] = ldar2scr(x,y,z,a);

% color code
cmap = colormap;
Lcmap = length(cmap);
cl = round(1 + (Lcmap - 1).*(t-t(1))./(t(end)-t(1)));

for i = 1:length(t)
    plot(N_x(i),N_y(i),'o','markersize',4,'color',cmap(cl(i),:),'markerfacecolor',cmap(cl(i),:))
end



function [N_x,N_y]=ldar2scr(x,y,z,a)

N_x=a.scr_x-(a.f*tan(atan((y-a.cam_y)./(x-a.cam_x ))- atan((a.ref_y-a.cam_y)./(a.ref_x-a.cam_x )) ))/a.p;

N_y=a.scr_y - (z - a.cam_z).*a.f./(a.p.*sqrt((x-a.cam_x ).^2+(y-a.cam_y ).^2 ))...
    ./sqrt(1+(tan(atan((y-a.cam_y)./(x-a.cam_x ))- atan((a.ref_y-a.cam_y)./(a.ref_x-a.cam_x )))).^2);


