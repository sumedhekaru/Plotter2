function data=linetExtract2(fn,t1,t2,rc,x0,y0,z0,disp)
%
%THIS PROGRAM INCLUDES A 19 MICROSEC OFFSET IN TIME TO ACCOUNT FOR A DELAY
%IN THE LINET DATA.  Oct.27, 2010, T. Marshall.  See data{3} in line 57
%below.
%
% clc
% 
% fn='F:\data/LINET/2010/07/11/linet_20100711_1930.txt';
% t1=0;
% t2=80000;
% 
% rc=500000;
% x0=0;
% 
% y0=0;
% z0=0;
% 2013-09-16 Updated from linetextract.m to handle peak currents

if nargin < 8
    disp = 1;
end

% Speed of light
c=299792458.0 ;


% Open and read all data in the file
fid = fopen(fn);
data=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f','HeaderLines',3);
fclose(fid);


% Total number of points
tot1=length(data{1});

% find occuring time in seconds - store at 1st column
data{1}=data{1}*3600+data{2}*60+data{3}+data{4}*(10^(-3))+data{5}*(10^(-6));

% Calculating x and y
ldar_lat=28.538486111;
ldar_lon=80.642633333;
earth_R=6371e3;

%10/27/10--Testing LINET way of calculating x,y from LINET lat,lon;
%previous way is commented out.  I note that the real problem was in the
%mistaken factor of pi./360 in the calculation of data{6}
% backup latitude longitude
data{12} = data{6};
data{13} = data{7};

% x
%data{6}=111200*cos((ldar_lat+data{7})./2.*pi./360).*(ldar_lon-(-data{6}));
data{6}=111319.491*cos((ldar_lat).*2.*pi./360).*(ldar_lon-(-data{6}));
% y
%data{7}=111000*(-ldar_lat+data{7});
data{7}=111319.491*(-ldar_lat+data{7});
% z 
data{8}=data{8}*1000; % was in km and converted in to m


% find distance to each pulse - store at 2nd column
data{2}=sqrt((data{6}-x0).^2+(data{7}-y0).^2+(data{8}-z0).^2);

% Detection time at the sensor - store at 3rd column
if x0==0 && y0==0 && z0 == 0
    data{3}=data{1}-(19e-6);
else
    data{3}=data{1}+data{2}/c-(19e-6);
end


try
    data=cell2mat(data);
catch
    % If above command failed, that mean the input file only contains 8
    % columns intead of 11 columns. Let's fix it.
    tmp = data{8}-data{8};
    data{9} = tmp;
    data{10} = tmp;
    data{11} = tmp;
    data=cell2mat(data);
end
    

%% Choose only data between t1 and t2
%lol=length(data(:,1))- nnz(data(:,1)>t1)+1 ;      
lol=sum(data(:,1)<t1)+1;                       % Index of the lower matrix element
ul=sum(data(:,1)<t2);                         % Index of the upper matrix element  

data=data(lol:ul,:);

% totol number of points in the time range
tot2=length(data(:,1));

%% Choose only data inside the radius rc
% sort data according to the distacne
data=sortrows(data,2);

ul=nnz(data(:,2)<=rc);
data=data(1:ul,:);

% Total number of points in the radius
tot3=length(data(:,1));

% finally convert all to array and sort by detecting time

if isempty(data)==0
    data=sortrows(data,1);
else
    data=[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
end





%% Printing Infomation
if disp
    fprintf('\nTOTAL NUMBER OF LINET POINTS :')
    fprintf('\n\tIn the file \t\t\t= %d\n',tot1)
    fprintf('\tIn the time range \t\t= %d (%.1f%%)\n',tot2,tot2*100/tot1)
    fprintf('\tIn the the radius %d km \t= %d (%.1f%%)\n',rc/1000,tot3,tot3*100/tot1)
end

%% LINET Column DISCRIPTION
% 1 - time
% 2 - Distance from the time correcting sensor
% 3 - time shifted time
% 6 - x distance
% 7 - y distance
% 8 - z distance
% 9 - Type
% 10 - Ip
% 11 - Err
% 12 - lat
% 13 - lon
