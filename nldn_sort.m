function nldn_sort

% NLDN data file from G. Lu came with all data. This program will devide
% data in to each day.

%% User Inputs
startDate = '2011-04-07';
endDate   = '2011-10-06';
fileName  = 'C:\Users\Sumedhe\Desktop\report7114557stroke Copy.txt';


%% Read the file
fID = fopen(fileName);
data=textscan(fID,'%s %f:%f:%f %f %f %f %s %s','HeaderLines',0);
fclose(fID);

d1 = datenum(startDate);
d2 = datenum(endDate);

d.dates = data{1};
d.hh = data{2};
d.mm = data{3};
d.ss = data{4};
d.lat = data{5};
d.lon = data{6};
d.Ip = data{7};
d.Ip_unit = data{8};
d.type = data{9};

ul = 0;
L = 1;

for i = d1:d2
    if L > 0
        lol = ul+1;
    end
    
    n = i-d1+1;
    dstr = {datestr(i,'yyyy-mm-dd')};
    
    L = sum(strcmp(d.dates,dstr));
    
    ul = lol + L - 1;    
    fprintf('%s\t%i\t%i\t%i\n',dstr{:},lol,ul,L);
    
    write_file(i,lol,ul,d);
end

function write_file(i,lol,ul,d)
 
% create the file name

dn = 'C:\Users\Sumedhe\Desktop\NLDN2\';


dateStr = datestr(i,'yyyymmdd');

year = dateStr(1:4);
month = dateStr(5:6);
date = dateStr(7:8);

dn = [dn year '\' month '\'];

if ~exist(dn,'dir')
    mkdir(dn)
end

fn = [dn 'NLDN2_' year month date '.txt'];

fID = fopen(fn,'a+');

for i = lol:ul;
    fprintf(fID,'%s %2.2i:%2.2i:%012.9f %08.4f %08.4f %0.1f %s %s\n',d.dates{i,:},d.hh(i),d.mm(i),d.ss(i), ...
        d.lat(i),d.lon(i),d.Ip(i),d.Ip_unit{i,:},d.type{i,:});
end

fclose(fID)