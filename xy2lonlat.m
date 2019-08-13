function [lon,lat] = xy2lonlat(x,y,disp,lon0,lat0)

if nargin < 1    
    x = -80620;
    y = 10956;
end
    

if nargin < 4
    lat0 = 28.538486111;
    lon0 = 80.642633333;
end

if nargin < 3
    disp = 0;
end

lat = y./111319.491+lat0;
lon = lon0-x./111319.491./cos(lat0.*2.*pi./360);

if nargout < 1 || disp == 1;
    fprintf('\tlon = %0.6f\n\laty = %0.6f\n',lon,lat)
end