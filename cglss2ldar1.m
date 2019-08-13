function cglss2ldar1

clc

wbh = waitbar(0,'Please Wait');

% This function convert the CGLSS files to normal files we were using

% Input file name
fn = 'C:\Users\Sumedhe\Desktop\KSCCGLSS20100801.dat';

% Output file folder
out_dir='C:\Users\Sumedhe\Desktop\ldar2/2010/08/01';

if exist(out_dir,'dir')==0    
    mkdir(out_dir)
end



waitbar(0,wbh,'Reading File');

% Open and read all data in the file
fid = fopen(fn);
data=textscan(fid,'%s %f:%f:%f %f:%f:%f %f:%f:%f %f %f %s %f %f %f %s','HeaderLines',0);
% Data contains the following
% 1 - date
% 2 - Hours
% 3 - minutes
% 4 - seconds
% 5,6,7 - Latitude
% 8,9,10 - Logitude
% 11 - Z
% 12 everything else
fclose(fid);

% Total number of points
tot1=length(data{1});

% find occuring time in seconds - store at 1st column
data1{1}=data{2}*3600+data{3}*60+data{4};

% Calculating x and y
ldar_lat=28.538486111;
ldar_lon=80.642633333;
earth_R=6371e3;

% find the x data
data1{2}=111319.491*cos((ldar_lat).*2.*pi./360).*(ldar_lon-(-(data{8}+data{9}/60+data{10}/3600)));
%data{6}=111319.491*cos((ldar_lat).*2.*pi./360).*(ldar_lon-(-data{6}));

% find y data
data1{3}=111319.491*(-ldar_lat+(data{5}+data{6}/60+data{7}/3600));


% Convert data in to mat format
data_new=cell2mat(data(2:11));
data1 = cell2mat(data1);

date=cell2mat(data{1}(1,1));

date = sprintf('%s%s%s',date(7:10),date(1:2),date(4:5));
daten = datenum(str2double(date(1:4)),str2double(date(5:6)),str2double(date(7:8))) ...
    - datenum(str2double(date(1:4)),0,0);




%% Calculating x,y,z

lat = (data{5}+data{6}/60+data{7}/3600)*pi/180;
lon = (data{8}-data{9}/60-data{10}/3600)*pi/180;
h = data_new(:,10); 

% data{8}(1,1)
% data{9}(1,1)
% data{10}(1,1)


h0 = 0;
lat0 = ldar_lat*pi/180;
lon0 = -ldar_lon*pi/180;


a = 6378137.0;
e = 0.081819190842613;

ki0 = sqrt(1-e^2*sin(lat0));
ki = sqrt(1-e^2*sin(lat));


X0 = (a/ki0+h0)*cos(lat0)*cos(lon0);
Y0 = (a/ki0+h0)*cos(lat0)*sin(lon0);
Z0 = (a*(1-e^2)/ki0+h0)*sin(lat0);

X = (a./ki+h).*cos(lat).*cos(lon);
Y = (a./ki+h).*cos(lat).*sin(lon);
Z = (a.*(1-e^2)./ki+h).*sin(lat);


x = -sin(lon0)*(X-X0) + cos(lon0)*(Y-Y0);
y = -sin(lat0)*cos(lon0)*(X-X0) - sin(lat0)*sin(lon0)*(Y-Y0) + cos(lat0)*(Z-Z0);
z =  cos(lat0)*cos(lon0)*(X-X0) + cos(lat0)*sin(lon0)*(Y-Y0) + sin(lat0)*(Z-Z0);


sec_1 = floor(data_new(:,3));
sec_2 = data_new(:,3) - sec_1;



%Lower limit
ll=0;
lln=1;
wbn=0;
%Write File names







for i = 0:47;   
    
    hh = floor(i/2);
    mm = rem(i,2)*30;
    
    wfn=sprintf('%s/ldar2_%s%3.3i%2.2i%2.2i.txt',out_dir,date(1:4),daten,hh,mm);
    fid1 = fopen(wfn,'w');
    
    ul=ll+1800;    
    uln = sum(data1(:,1)<=ul);
    
    str= sprintf('Sumedhe\nMade\nThe\nFollowing\nLdar\nFile\nFrom\nCGLSS\ndata using cglss2ldar1.m matlab program \n\nJDAY\tTime(UTC)\tX(M)\t(Y)\tZ(M)\tEVENT\tTYPE\n\n');
    
    % Writing the header 
    fprintf(fid1,str);
    

    
    for j=lln:uln
        fprintf(fid1,'%i\t%2.2i:%2.2i:%2.2i:%0.6f\t%.0f\t%.0f\t%.0f\tCGLSS\tEVENT\n',...
            daten,data_new(j,1),data_new(j,2),sec_1(j),sec_2(j),x(j),y(j),z(j));        
    end
    
    fclose(fid1);
    
    if uln < lln
        delete(wfn)
    else
        lln=uln;
    end
    ll=ul;
    
    
    
    str=sprintf('Now writing...\nlinet-%.0f-%2.2i%2.2i.txt',date,hh,mm);
        
    try
        waitbar(i/48,wbh,str)
    catch
        fprintf('\nOperation Terminated by User\nV1.0 - Sumedhe \n')
        return
    end  
end



delete(wbh)