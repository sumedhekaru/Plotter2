function read_Ip_Ep_r_file
clc

fn = 'C:\Users\sumedhe\Desktop\Ip_vs_Epxr\Ip_Ep_data_20110814_1800-2400.txt';

% Open and read all data in the file
fid = fopen(fn);
data=textscan(fid,'%f %f %s %f %f %f');
% Data contains the following
% 1 - index
% 2 - seconds from midnight
% 3 - sensor name and frequency
% 4 - Peak electric field
% 5 - Distance to the sensor
% 6 - peak current from CGLSS

fclose(fid);

%% Get only data between R1 and R2
R1 = 30;
R2 = 80;
% sort data according to distance
[ix,ix] = sort(data{5});

for i=1:6
    data{i} = data{i}(ix,:);
end

lol=nnz(data{5}<R1)+1 ;      % Index of the lower matrix element
ul= nnz(data{5}<R2)  ;       % Index of the upper matrix element  

for i=1:6
    data{i}=data{i}(lol:ul);
end


%% Get only the current only between I1 and I2
I1= -500;
I2 = 500;

% sort data according to current
[ix,ix] = sort(data{6});

for i=1:6
    data{i} = data{i}(ix,:);
end

lol=nnz(data{6}<I1)+1 ;      % Index of the lower matrix element
ul= nnz(data{6}<I2)  ;       % Index of the upper matrix element  

for i=1:6
    data{i}=data{i}(lol:ul);
end


%% Sort data in to sensors
L = length(data{1});

k02 = nan(L,3);
c=0;

sen = 'STC:ch3:S-0.50MHz';


fac=1;
% fac = 1350.25/9.94; %K02
% fac = 1214.5/10.18;	%K14
% fac = 1456/11.01;	%K24
% fac = 1241.573/18.85; %BCC
% fac = 1456/11.01;  %K17
% fac = 1804/8.76; % EDW
 fac = 551.57/8.65; %STC
%fac = 444.931/8.33; % FLT
%fac = 559.9/8.81; % OVD


for i = 1:L
    str = data{3}(i);
    
    if strcmp(str,sen) && abs(data{4}(i)) > 4
        c = c + 1;
        k02(c,1)=data{6}(i);
        k02(c,2)=data{4}(i)/fac;
        k02(c,3)=data{5}(i);
        
    end
    
    k02(c+1:end,:)=[];
end

%N = [1 1.13];
N=[1.13];
%str = sprintf('R = %0.0f-%0.0fkm    I = %0.0f - %0.0fkA   %s',R1,R2,I1,I2,sen);
str = sprintf('R = %0.0f-%0.0fkm   %s',R1,R2,sen);

peakCurrBestFit(abs(k02(:,1)),abs(k02(:,2)),k02(:,3),N,str)



Ip = abs(data{6});
Ep = abs(data{4});
R  = data{5};
%N = [1 1.13];
N=[];
str = sprintf('R = %0.0f-%0.0fkm     I = %0.0f - %0.0fkA',R1,R2,I1,I2);

%peakCurrBestFit(Ip,Ep,R,N,str)