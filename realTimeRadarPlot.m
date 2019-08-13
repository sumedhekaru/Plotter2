function realTimeRadarPlot(stationId)

if nargin < 1
    % Memphis KNQA
    stationId = 'KNQA';
end

% Base radar downloading site (obtained from GRAnalysit)
base_site = 'http://mesonet-nexrad.agron.iastate.edu/level2/raw/';

% Base folder for saving data
rbf = 'C:\Users\sumedhe\Desktop\RealTimeRadar\';


%% Start the program
str = urlread([base_site stationId]);

indx = strfind(str,stationId);

% Just before last file (last file might be still writing).
lastFile = str(indx(end-2):indx(end-2)+17);

year = lastFile(6:9);
month = lastFile(10:11);
date = lastFile(12:13);
hour = lastFile(15:16);
mm = lastFile(17:18);

% Download this file

bf1 = sprintf('%s/%s/Compressed/%s/%s/%s/%s/',rbf,stationId,year,month,date,hour);

if ~exist(bf1,'dir')
    mkdir(bf1)
end

% If not already downloaded, let's download it
if ~exist([bf1 lastFile],'file')
    urlwrite([base_site stationId '/' lastFile],[bf1 lastFile]);
end

bf2 = sprintf('%s/%s/%s/%s/%s/%s/',rbf,stationId,year,month,date,hour);

if ~exist(bf2,'dir')
    mkdir(bf2)
end


% If have not Extracted it, let's do it
if ~exist([bf2 lastFile '.netcdf'],'file')
    cmd = sprintf('"C:/Program Files/Java/jre6/bin/java" -classpath "C:/Users/sumedhe/Desktop/Plotter2/toolsUI-4.3.jar" ucar.nc2.FileWriter -in "%s%s" -out "%s%s.netcdf"',...
        bf1,lastFile,bf2,lastFile);
    system(cmd)
end


% All done, let's plot it
sen_set = open('sensor_setting.mat');
sen_set.radarFn = [bf2 lastFile '.netcdf'];
sen_set.radEleAngInd = 1;
radar_plot(sen_set)





