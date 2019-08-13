function download_radar_data

% Firts you should request the data from www.ncdc.noaa.gov and get the link
% and add that link to the following link variable.

%% User input
link = 'http://ftp3.ncdc.noaa.gov/pub/has/HAS010340788/';

% Download folder
dbf = 'C:\Users\Sumedhe\Desktop\Radar\';

% Version of data (2 or 3)
version = 2;

%% Program
data = urlread(link);

indx = strfind(data,'KMLB');

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

delete(wbh)
