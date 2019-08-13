function [x,y] = latlon2xy(lat, lon, disp, lat0,lon0);

if nargin < 1    
    lat = 44.02;
    lon = 103.29;
end
    

if nargin < 4
    %lat0 = 44.033642;
    %lon0 = 103.285475;
    lat0 = 28.538486111;
    lon0 = 80.642633333;
end

if nargin < 3
    disp = 0;
end


% Calculating x and y and z
x = 111319.491*cos((lat0).*2.*pi./360).*(lon0-(lon));
y = 111319.491*(-lat0+lat);

if nargout < 1 || disp == 1;
    fprintf('\tx = %0.1f\n\ty = %0.1f\n',x,y)
end


