function download_radar_data

% Firt you should request the data from http://www.ncdc.noaa.gov/oa/radar/radardata.html 
% and request level 2 data and get the link.
% and add that link to the following link variable.
% select the correct Download folder dbf
%run nextRAD2netCDF.m to convert the file into netcdf

%% User input

%link = 'http://www1.ncdc.noaa.gov/pub/has/HAS010655463/0004/';
% link = 'http://www1.ncdc.noaa.gov/pub/has/HAS010656086/0004/';
 %link = 'http://www1.ncdc.noaa.gov/pub/has/HAS010674871/0004/';
%link = 'http://www1.ncdc.noaa.gov/pub/has/HAS010674818/';
link = 'http://www1.ncdc.noaa.gov/pub/has/HAS011013285/0002/';
% Download folder
 dbf = 'C:\Users\daqop\Desktop\radar_data_2016\RADAR\MEM\';
%dbf = 'C:\Users\daqop\Desktop\radar_data_2016\';
%dbf = 'C:\data\2010-07-01 -- 2010-08-19\data\RADAR\MEL\';


% Version of data (2 or 3)
version = 2;

% Station
%station = 'KJAX'; % Tampa
%station = 'KAMX'; % Miami
%station = 'KMLB';  % Melborn

station = 'KNQA';  % Memphis
%station = 'KGWX';  % colombus
%% Program
data = urlread(link);


indx = strfind(data,station);

L = length(indx);

wbh = waitbar(0,'Downloading RADAR data','name','RADAR');

if version == 3    
    for i = 1:L
       
        fn = data(indx(i):indx(i)+30);
        
        try
            waitbar(i/L,wbh)
        catch
            disp('Data download interuppted.')
            return
        end
        
        % Create the correct folder for data
        year = fn(20:23);
        month = fn(24:25);
        date = fn(26:27);
        hour = floor(str2double(fn(28:29)));
        
        bf = sprintf('%s/%s/%s/%s/%2.2i/',dbf,year,month,date,hour);
        
        if ~exist(bf,'dir')
            mkdir(bf)
        end

        urlwrite([link fn],[bf fn]);
        
    end
end

if version == 2 
    
    for i = 1:L
        tic
        
        fn = data(indx(i):indx(i)+25);
        
        try
            waitbar(i/L,wbh,['Downloading ' fn '...'])
        catch
            disp('Data download interuppted.')
            return
        end
        
        % Create the correct folder for data
        year = fn(5:8);
        month = fn(9:10);
        date = fn(11:12);
        hour = floor(str2double(fn(14:15)));
        
        bf = sprintf('%s/%s/%s/%s/%2.2i/',dbf,year,month,date,hour);
        
        if ~exist(bf,'dir')
            mkdir(bf)
        end
        
        % download data
        urlwrite([link fn],[bf fn]);

        % unzip and delete original zip file
        gunzip([bf fn],bf)
        delete([bf fn])        
        
    end
end

fprintf('Done downloading RADAR data')

delete(wbh)
