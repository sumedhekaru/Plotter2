% finding position from LDAR points to sensors
%clc
a=open('sensor_setting.mat');

% Sensor number 
sn2 = 11;

% Tolerance percentage
tp =5;

h=guidata(findall(0,'Tag','plotter2'));
g = h.g;

a.t1=g.t1;
a.t2=g.t2;

% Overwrite t1 and t2?
t_over = 1;

if t_over
    a.t1 =  82800;
    a.t2 =  a.t1+1800;
end

%% Generating file name for LDAR data
if g.mm < 30
    ext=0;
else
    ext=30;
end

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
    a.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);

if exist(ldar_fn,'file')~=0
    % Check whether the file is exists
    a.ldar_fn=ldar_fn;
else
    % If file is not exist don't store the file name
    a.ldar_fn='';
    % Absant file
   % a.absant_fn=[a.absant_fn ldar_fn];
end

ld_lin='';

% For time shift
x0=0; y0=0; z0=0;


%waitbar(0.7,wbh,'Loading LDAR data','Name','Plotter Busy')
[CG,CAL,DLS]=ldarExtract2(ldar_fn,a.t1,a.t2,str2double(a.ldar_r),...
    x0,y0,z0,0);

% t1=[DLS(:,10);CG(:,10)];
% x1=[DLS(:,6);CG(:,6)];
% y1=[DLS(:,7);CG(:,7)];
% z1=[DLS(:,8);CG(:,8)];

t1=CG(:,10);
x1=CG(:,6);
y1=CG(:,7);
z1=CG(:,8);


%     %For future use to extract ldar data in to the workspace
%     assignin('base', 'lt', x1)
%     assignin('base', 'lx', [DLS(:,6);CG(:,6)])
%     assignin('base', 'ly', [DLS(:,7);CG(:,7)])
%     assignin('base', 'lz', y1)

% Arrays to store the same distance results
xxk02=[];
yyk02=[];

xxk14=[];
yyk14=[];

xxk17=[];
yyk17=[];

xxk24=[];
yyk24=[];

% counter for the same distance results
i14 = 1;
i24 = 1;
i02 = 1;
i17 = 1;

% legend 
lg = {};

fprintf('Sensor\tTime\t\tx\t\ty\tz\tTo sensor\tTo %s',a.sen_IDs{sn2})
for i = 1:length(x1)
   lk02 = sqrt((x1(i)-a.x(1))^2+(y1(i)-a.y(1))^2);
   lk14 = sqrt((x1(i)-a.x(2))^2+(y1(i)-a.y(2))^2);
   lk17 = sqrt((x1(i)-a.x(6))^2+(y1(i)-a.y(6))^2);
   lk24 = sqrt((x1(i)-a.x(3))^2+(y1(i)-a.y(3))^2);
   l2 = sqrt((x1(i)-a.x(sn2))^2+(y1(i)-a.y(sn2))^2);
   
   
   if abs(lk02-l2) <= (lk02+l2)/2*tp/100
       fprintf('\nK02:\t%fs\t%.1fm\t%.1fm\t%.1fm\t%.1fm\t\t%.1fm',t1(i),x1(i),y1(i),z1(i),lk02,l2)
       xxk02(i02) = x1(i);
       yyk02(i02) = y1(i);       
       i02= i02+1;
   end
   
   if abs(lk14-l2) <= (lk14+l2)/2*tp/100
       fprintf('\nK14:\t%fs\t%.1fm\t%.1fm\t%.1fm\t%.1fm\t\t%.1fm',t1(i),x1(i),y1(i),z1(i),lk14,l2)
       xxk14(i14) = x1(i);
       yyk14(i14) = y1(i);
       i14= i14+1;
   end
   
    if abs(lk17-l2) <= (lk17+l2)/2*tp/100
       fprintf('\nK17:\t%fs\t%.1fm\t%.1fm\t%.1fm\t%.1fm\t\t%.1fm',t1(i),x1(i),y1(i),z1(i),lk17,l2)
       xxk17(i17) = x1(i);
       yyk17(i17) = y1(i);
       i17= i17+1;
   end
   
   
   if abs(lk24-l2) <= (lk24+l2)/2*tp/100
       fprintf('\nK24:\t%fs\t%.1fm\t%.1fm\t%.1fm\t%.1fm\t\t%.1fm',t1(i),x1(i),y1(i),z1(i),lk24,l2)
       xxk24(i24) = x1(i);
       yyk24(i24) = y1(i);
       i24= i24+1;
   end    
end

fprintf('\n')

figure
ind = find(a.x ~= 0);
plot(a.x(ind)/1000,a.y(ind)/1000,'rp','MarkerFaceColor','r')
text(a.x(ind)/1000,a.y(ind)/1000,a.sen_IDs(ind))
lg = [lg 'Sensors'];
hold all

if  isempty(x1) == 0
    plot(x1/1000,y1/1000,'go','MarkerFaceColor','g','MarkerSize',2)
    lg = [lg ,'CG flshes'];
end
hold all
if isempty(xxk02) == 0
    plot(xxk02/1000,yyk02/1000,'bo','MarkerFaceColor','b','MarkerSize',2)
    lg = [lg ,'with K02'];      
end
hold all

if isempty(xxk14) == 0
    plot(xxk14/1000,yyk14/1000,'ko','MarkerFaceColor','k','MarkerSize',2)
    lg = [lg ,'with K14'];      
end
hold all

if isempty(xxk17) == 0
    plot(xxk17/1000,yyk17/1000,'mo','MarkerFaceColor','m','MarkerSize',2)
    lg = [lg ,'with K17'];      
end
hold all

if isempty(xxk24) == 0
    plot(xxk24/1000,yyk24/1000,'ro','MarkerFaceColor','r','MarkerSize',2)
    lg = [lg ,'with K24'];      
end
hold all



set(gca,'DataAspectRatio',[1 1 1])
box on
xlabel('East(km)')
ylabel('North(km)')
tit = sprintf('Flashes that have same distance from %s',a.sen_IDs{sn2});
title(tit)
legend(lg)

florida_map


