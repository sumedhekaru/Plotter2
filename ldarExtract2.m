function [CGnn,CALnn,DLSnn]=ldarExtract2(fn,t1,t2,rc,x0,y0,z0,graphon,output)

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
c=299792458.0;


% Open and read all data in the file
fid = fopen(fn);
data=textscan(fid,'%f %f:%f:%f:%f %f %f %f %s %s','HeaderLines',12);
fclose(fid);

% Total number of pints
tot1=length(data{1});

% find occuring time in seconds - store at 1st column
data{1}=data{2}*3600+data{3}*60+data{4}+data{5}*(10^(-6));

% % find distance to each pulse - store at 10th column
temp =sqrt((data{6}-x0).^2+(data{7}-y0).^2+(data{8}-z0).^2);

% find horizontal distance
data{11} = sqrt((data{6}).^2+(data{7}).^2);

% Detection time at the sensor - store at 10th column
if x0 == 0 && y0 == 0 && z0 == 0
    data{10}=data{1};
else
    data{10}=data{1}+temp/c;
end

%% Choose only data between t1 and t2
lol=length(data{10})- nnz(data{10}>t1)+1 ;      % Index of the lower matrix element
ul=nnz(data{10}<t2)  ;                         % Index of the upper matrix element  

for i=1:11
    data{i}=data{i}(lol:ul);
end

% totol number of points in the time range
tot2=length(data{1});

%% Choose only data inside the radius rc
% sort data according to the distacne

[ix,ix] = sort(data{11});

for i=1:11
    data{i} = data{i}(ix,:);
end

ul=nnz(data{11}<=rc);

for i=1:11
    data{i}=data{i}(1:ul);
end

% Total number of points in the radius
tot3=length(data{1});


%% Choose the data inside the box
sen_set = open('sensor_setting.mat');
tot4 = NaN;
if sen_set.ldar_box(1)
    lb = sen_set.ldar_box;
    
    % Choose points within x range
    inds1 = find(data{6} > lb(2)*1000);
    inds2 = find(data{6} < (lb(2) + lb(4))*1000);
    xinds = intersect(inds1,inds2);
    
    % Choose points within y range
    inds1 = find(data{7} > lb(3)*1000);
    inds2 = find(data{7} < (lb(3) + lb(5))*1000);
    yinds = intersect(inds1,inds2);
    
    inds = intersect(xinds,yinds);
    
    for i=1:11
        data{i}=data{i}(inds);
    end
    
    tot4=length(data{1});
end


%% Choose data by type

[ix,ix] = sort(data{9});

for i=1:11
    data{i} = data{i}(ix,:);
end

nIC=sum(strcmp('4DLSS',data{9}));
nCAL=sum(strcmp('CAL',data{9}));
nCG=sum(strcmp('CGLSS',data{9}));

for i=1:11
    DLSnn{i}=data{i}(1:nIC);
    CALnn{i}=data{i}(nIC+1:nIC+nCAL);
    CGnn{i}=data{i}(nIC+nCAL+1:end);
end

clear data

% Occuring time should be at column 9 but not column 1
DLSnn{9}=DLSnn{1};
CALnn{9}=CALnn{1};
CGnn{9}=CGnn{1};

% finally convert all to array and sort by detecting time
DLSnn=cell2mat(DLSnn);
CALnn=cell2mat(CALnn);
CGnn=cell2mat(CGnn);


if isempty(DLSnn)==0
    DLSnn=sortrows(DLSnn,10);
else
    DLSnn=[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
end

if isempty(CALnn)==0
    CALnn=sortrows(CALnn,10);
else
    CALnn=[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
end

if isempty(CGnn)==0
    CGnn=sortrows(CGnn,10);
    % Impose all z coordinates to zero
    CGnn(:,8)=CGnn(:,8)-CGnn(:,8);

else
    %CGnn=[0,0,0,0,0,0,0,0,0,0,0];
    CGnn=[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
end

%% Printing Infomation
if output
    fprintf('\nTOTAL NUMBER OF LDAR POINTS :')
    fprintf('\n\tIn the file \t\t\t= %d\n',tot1)
    fprintf('\tIn the time range \t\t= %d (%.1f%%)\n',tot2,tot2*100/tot1)
    fprintf('\tIn the radius %d km \t= %d (%.1f%%)\n',rc/1000,tot3,tot3*100/tot1)
    
    if ~isnan(tot4)
        fprintf('\tIn the box (%0.1f,%0.1f,%0.1f,%0.1f)km\t= %d (%.1f%%)\n',lb(2:5),tot4,tot4*100/tot1)
    end
    
    
    fprintf('\n In side the radius %d km',rc/1000)
    fprintf('\n\t Number of CG points \t\t= %d',nCG)
    fprintf('\n\t Number of IC points \t\t= %d',nIC)
    fprintf('\n\t Number of CAL points \t\t= %d\n',nCAL)
end




