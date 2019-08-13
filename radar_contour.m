function radar_contour(fileName,eleAngNum)
% This function is use to draw a contour line arround a user clicked
% positon.
%   Usage:
%       (1) Draw a radar plot
%       (2) Go to Radar in the menu and select "Radar Contour" or run
%       'radar_contour' in matlab command window or run this file
%       (3) If the current azimuth angle info available in figure handle,
%       you are done! If not it will ask you wich angle you want to plot.
% 
%   Modification History
%       2014-07-09 Created by Sumedhe Karunarathne
%

 a = guidata(gcf);

if nargin < 1
    % File name not provided. Try to obtain if it is available in the
    % figure handle.
    if isempty(a) || ~isfield(a,'fileName')
        msgbox('No file name found, please provide file name manually')
        return
    end
    
    % If it come this far, we have a file name. Let's see it is in the disk
    if ~exist(a.fileName,'file')
        msgbox(a.fileName,'File not found!')
        return
    end
    
    fileName = a.fileName;
end

   
% If it come this far, we have a file. Now we need a radar angle. 
if nargin < 2
    if isempty(a) || ~isfield(a,'eleAngNum')
        answer = inputdlg('Enter Radar Elevation Angle Number (Ex. 5)');        
        eleAngNum = str2double(answer);
        
        if isnan(eleAngNum) || isempty(eleAngNum)
            msgbox('Error occure while getting Radar Elevation Angle.')
            return
        end
    else
        eleAngNum = a.eleAngNum;
    end
end

%eleAngNum = 11;


% Ok, now we have file name and elevation angle. Let's get the (x,y)
% coordinates (in km)
[x_click,y_click] = ginput(1);
%x_click = 12.2;
%y_click = -12.1;



% Open the radar file
ncid1 = netcdf.open(fileName,'nowrite');

% Station latitude, longitude
stLat = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),'StationLatitude','double');
stLon = netcdf.getAtt(ncid1,netcdf.getConstant('NC_GLOBAL'),'StationLongitude','double');

% In horizontal coordinates
[x0, y0] = latlon2xy(stLat,-stLon);
x0 = x0/1000;
y0 = y0/1000;

% Determine which number (horizontal scan number) is to load
if eleAngNum == 1
    eleAngVar = 'elevationR_HI';
    k = 1;
elseif eleAngNum == 2
    eleAngVar = 'elevationR_HI';
    k = 3;
else
    eleAngVar = 'elevationR';
    k = eleAngNum-2;
end

varid = netcdf.inqVarID(ncid1,eleAngVar);
eleAngs = netcdf.getVar(ncid1,varid);

% Elevation angles
try
    eleAng = mean(eleAngs(:,k));
catch
    msgbox('This angle does not exist')
    return
end


netcdf.close(ncid1);
hold all
R = sqrt((x0 - x_click)^2+(y0 - y_click)^2);

% Draw a cirle with R
cA = 0:0.01:2*pi;
cX = R*cos(cA);
cY = R*sin(cA);
plot(cX+x0,cY+y0,'k--','linewidth',1)

% Calculate the altitude of the circle
Z = R*tan(eleAng*pi/180);

% add the text to the cirlce
text(x_click,y_click',sprintf('%0.2f km',Z))


set(gcf,'renderer','painters')


