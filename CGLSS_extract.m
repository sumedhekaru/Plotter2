function data=CGLSS_extract(fn,t1,t2,rc,x0,y0,z0)
% fn = 'H:\data\cglss\2011\08\KSCCGLSS20110814.dat';
% t1 = 79781.574702;
% t2 = 79785.575860;
% rc = 100000;
% x0 = 1000;
% y0 = 1000;
% z0 = 0;


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
% 11 - Current
% 12 - Z
% 13:16 - don't know
% 17 - sensor's used


fclose(fid);


% find occuring time in seconds - store at 1st column
data{1}=data{2}*3600+data{3}*60+data{4};


% Calculating x and y
ldar_lat=28.538486111;
ldar_lon=80.642633333;
% earth_R=6371e3;
c = 299704764; % Spped of light in air

% find the x data
data{2}=111319.491*cos((ldar_lat).*2.*pi./360).*(ldar_lon-(-(data{8}-data{9}/60-data{10}/3600)));
%data{6}=111319.491*cos((ldar_lat).*2.*pi./360).*(ldar_lon-(-data{6}));

% find y data
data{3}=111319.491*(-ldar_lat+(data{5}+data{6}/60+data{7}/3600));

%% Calculating detection time
%Distance from x0,y0,z0
data{4} = ((x0 - data{2}).^2 + (y0 - data{3}).^2 + (z0 - data{12}).^2).^0.5;

if x0 == 0 && y0 == 0 && z0 == 0
    data{13} = data{1};
else
    data{13} = data{1} + data{4}./c;
end
    

%% I added following line to avoid possible data range miss identifications
lol= nnz(data{1} < t1 - 100) + 1 ;      % Index of the lower matrix element
ul = nnz(data{1} < t2 + 100)  ;                         % Index of the upper matrix element  

for i=1:17
    data{i}=data{i}(lol:ul);
end



%% Choose only data between t1 and t2
lol= nnz(data{1} < t1) + 1 ;      % Index of the lower matrix element
ul = nnz(data{1} < t2)  ;                         % Index of the upper matrix element  

temp = data;
temp{17} = zeros(size(temp{17}));

for i=1:17
    data{i}=data{i}(lol:ul);
end

%% Choose only data within radius r
data{15} = (data{2}.^2 + data{3}.^2).^0.5;
[ix,ix] = sortrows(data{15});

for i=1:17
    data{i} = data{i}(ix,:);
end

ul=nnz(data{13}<=rc);

for i=1:17
    data{i}=data{i}(1:ul);
end

% Sort back according to time
[ix,ix] = sortrows(data{1});

for i=1:17
    data{i} = data{i}(ix,:);
end



% Convert data in to mat format
last_column = cellstr(data{17}); % Backup of column 17
data{17} = zeros(size(data{17}));
data = cell2mat(data);


% Column 17 will contain the number of sensors used
for i = 1:length(last_column)
    data(i,17)= ceil(length(char(last_column(i)))/2);
end

if isempty(data)    
    data = nan(1,17);
end
 
    
%     temp = cell2mat(temp);
%     try
%         for  i = lol-1:ul+1
%             fprintf('%0.6f\t%0.1f\n',temp(i,1),temp(i,11))
%         end
%     catch
%         do nothing
%     end
    

%% Important column definition at the end
% Total number of colums = 17
% 1 : Occuring time
% 2 : x
% 3 : y
% 12: z
% 11: Current
% 4 : distance to x0,y0,z0 given
% 5,6,7 - Latitude
% 8,9,10 - Logitude
% 13 - detection time
% 17 - Total number of sensors used

