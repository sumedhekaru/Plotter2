function data=pbfaExtract(fn,t1,t2,rc,x0,y0,z0,output)
%Points by fast anttenna extract
% Command: data=pbfaExtract(fn,t1,t2,rc,x0,y0,z0,output);
% Where data contains
%   data{1} - distance to the pulse from the given sensor 
%   data{2} - Occuring time
%   data{3} - x
%   data{4} - y
%   data{5} - z
%   data{6} - Detection time at the sensor
%   data{7} - I
%   data{8} - 
%   data{9} - dt
%   data{10} - dx
%   data{11} - dy
%   data{12} - dz
%   data{13} - ki_sqrd


% fn='PBFA_20100817_1600.txt';
% 
% t1=0;
% t2=80000;
% 
% rc=500000;
% x0=0;
% 
% y0=0;
% z0=0;

try
    temp = output;
catch
    output = 1;
end

% Speed of light
c=299792458.0;

% Open and read all data in the file
fid = fopen(fn);
data=textscan(fid,'%f %f %f %f %f %f %f %s %f %f %f %f %f','HeaderLines',2);
fclose(fid);

% data{1} - flash ID
% data{2} - Occuring time
% data{3} - x
% data{4} - y
% data{5} - z
% data{6} - N  If 10<N<20 --> Calculated by hand
% data{7} - I
% data{8} - comma sepeated string containing the sensors used
% data{9} - dt
% data{10} - dx
% data{11} - dy
% data{12} - dz
% data{13} - ki_sqrd

data{8} = data{6};

% find distance to each pulse - store at 1nd column
data{1}=sqrt((data{3}-x0).^2+(data{4}-y0).^2+(data{5}-z0).^2);

% Detection time at the sensor - store at 6th column
if x0 == 0 && y0 == 0 && z0 == 0
    data{6}=data{2};
else
    data{6}=data{2}+data{1}/c;
end

try
    data=cell2mat(data);
catch
    % Old PBFA files didn't have 13 columns. If above command failing, it 
    % could be this. Let's have 13 column.
    data{13} = nan(size(data{1}));
    data=cell2mat(data);
end

% After this point
%   data{1} - distance to the pulse from the given sensor 
%   data{2} - Occuring time
%   data{3} - x
%   data{4} - y
%   data{5} - z
%   data{6} - Detection time at the sensor

% Total number of points
tot1=length(data(:,1));

%% Choose only data between t1 and t2
data=sortrows(data,2);


lol=nnz(data(:,6)<t1)+1  ;     % Index of the lower matrix element
ul=nnz(data(:,6)<t2)  ;                        % Index of the upper matrix element  

data=data(lol:ul,:);

% totol number of points in the time range
tot2=length(data(:,1));

%% Choose only data inside the radius rc
% sort data according to the distacne
data=sortrows(data,1);

ul=nnz(data(:,1)<=rc);
data=data(1:ul,:);

% Total number of points in the radius
tot3=length(data(:,1));

%xx = data(:,3)


%% Choose the data inside the box
sen_set = open('sensor_setting.mat');
tot4 = NaN;
if sen_set.ldar_box(1)
    lb = sen_set.ldar_box;
    
    % Choose points within x range
    inds1 = find(data(:,3) > lb(2)*1000);
    inds2 = find(data(:,3) < (lb(2) + lb(4))*1000);
    xinds = intersect(inds1,inds2);
    
    % Choose points within y range
    inds1 = find(data(:,4) > lb(3)*1000);
    inds2 = find(data(:,4) < (lb(3) + lb(5))*1000);
    yinds = intersect(inds1,inds2);
    
    inds = intersect(xinds,yinds);
    
    data = data(inds,:);
    
    tot4=length(inds);
end

% finally convert all to array and sort by detecting time
if isempty(data)==0
    data=sortrows(data,6);
else
    data=[NaN,NaN,NaN,NaN,NaN,NaN];
end

%% Printing Infomation
if output
    fprintf('\nTOTAL NUMBER OF PBFA POINTS :')
    fprintf('\n\tIn the file \t\t\t= %d\n',tot1)
    fprintf('\tIn the time range \t\t= %d (%.1f%%)\n',tot2,tot2*100/tot1)
    fprintf('\tIn the the radius %d km \t= %d (%.1f%%)\n',rc/1000,tot3,tot3*100/tot1)
    if ~isnan(tot4)
        fprintf('\tIn the box (%0.1f,%0.1f,%0.1f,%0.1f)km\t= %d (%.1f%%)\n',lb(2:5),tot4,tot4*100/tot1)
    end
end

