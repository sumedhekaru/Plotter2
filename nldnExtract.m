function [NLDNc NLDNg]=nldnExtract(fn,lat0,lon0,t1,t2,rc,x0,y0,z0,output)

% fn = 'C:\Users\Sumedhe\Desktop\NLDN2\2011\08\NLDN2_20110801.txt';
% t1 = 0;
% t2 = 86400;
% rc = 100000000000;
% x0 = 0;
% y0 = 0;
% z0 = 0;
% output = 0;

% west 35
lat0 = 44.033642;
lon0 = 103.285475;

% E35 
%lat0 = 44.061394;
%lon0 = 103.168130;



%%%%%% Column Discription
% 1 - occuring time
% 6 - x
% 7 - y
% 8 - z
% 9 - Occuring time
% 10 - Detection time
% 11 - horizontal Distance to each pulse

try 
    temp = output;
catch
    output = 1;
end

% Speed of light
c=299792458.0 ;


% Open and read all data in the file
fid = fopen(fn);
data=textscan(fid,'%s %f:%f:%f %f %f %f %s %s','HeaderLines',12);
fclose(fid);

% Total number of pints
tot1=length(data{1});

% find occuring time in seconds - store at 1st column
data{1}=data{2}*3600+data{3}*60+data{4};

% Calculating x and y and z
data{2} = 111319.491*cos((lat0).*2.*pi./360).*(lon0-(-data{6}));
data{3} = 111319.491*(-lat0+data{5});
data{4} = strcmp(data{9},'C')*5000;  % all the cloud points would be 5000 m

% find distance to each pulse 
temp =sqrt((data{2}-x0).^2+(data{3}-y0).^2+(data{4}-z0).^2);

% find horizontal distance from the origin
data{9} = sqrt((data{2}).^2+(data{3}).^2);

% Detection time at the sensor - store at 6th column
if x0 == 0 && y0 == 0 && z0 == 0
    data{8}=data{1};
else
    data{8}=data{1}+temp/c;
end

%% Choose only data between t1 and t2
lol=length(data{1})- nnz(data{1}>t1)+1 ;      % Index of the lower matrix element
ul=nnz(data{1}<t2)  ;                         % Index of the upper matrix element  

for i=1:9
    data{i}=data{i}(lol:ul);
end

% totol number of points in the time range
tot2=length(data{1});

%% Choose only data inside the radius rc
% sort data according to the distacne

[ix,ix] = sort(data{9});

for i=1:9
    data{i} = data{i}(ix,:);
end

ul=nnz(data{9}<=rc);

for i=1:9
    data{i}=data{i}(1:ul);
end

% Total number of points in the radius
tot3=length(data{1});

%% seperate data type 
[ix,ix] = sort(data{4});

for i=1:9
    data{i} = data{i}(ix,:);
end

nCG = nnz(data{4} == 0);


for i=1:9
    NLDNg{i}=data{i}(1:nCG);
    NLDNc{i}=data{i}(nCG+1:end);
end


clear data;

%% finally convert all to array and sort by detecting time
NLDNc = cell2mat(NLDNc);
NLDNg = cell2mat(NLDNg);

if isempty(NLDNc)==0
    NLDNc=sortrows(NLDNc,1);
else
    NLDNc=[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
end

if isempty(NLDNg)==0
    NLDNg=sortrows(NLDNg,1);
else
    NLDNg =[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
end

%% Data output columns
% 1 - detection time in seconds from midnight
% 2 - x
% 3 - y
% 4 - z
% 5 - Lat
% 6 - Lon
% 7 - Ip (in kA)
% 8 - arrival time
% 9 - horizontal distance

% figure
% hold all
% plot(NLDNg(:,8),NLDNg(:,4),'ro','markerfacecolor','r')
% plot(NLDNc(:,8),NLDNc(:,4),'bo','markerfacecolor','b')
% 
% %plot(0,0,'ko')
% axis square

%% Printing Infomation
if output
    fprintf('\nTOTAL NUMBER OF NLDN POINTS :')
    fprintf('\n\tIn the file \t\t\t= %d\n',tot1)
    fprintf('\tIn the time range \t\t= %d (%.1f%%)\n',tot2,tot2*100/tot1)
    fprintf('\tIn the the radius %d km \t= %d (%.1f%%)\n',rc/1000,tot3,tot3*100/tot1)
end



