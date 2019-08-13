function plot_scan(nc_file)
% plot_scan -- Plot the variables in a Chilbolton scanning NetCDF file
%
% plot_scan('nc_file') plots the variables in the file 'nc_file' in a
%   new A4 figure window.
%
%   To use this function you need the NetCDF Toolbox for
%   Matlab-5. This can be downloaded free from:
%     http://crusty.er.usgs.gov/~cdenham/MexCDF/nc4ml5.html

if nargin < 1, help(mfilename), return, end

f = netcdf(nc_file, 'nowrite');
if isempty(f)
  error(['Could not open ' nc_file 'as a NetCDF file']);
  return
end

% Colour axes for known 2D variables
known_var = {'Zh', 'Zdr', 'Ldr', 'v', 'v_unfolded', 'phi_dp', ...
	     'Z', 'v1', 'v2'};
known_min = [-20 -1 -40 -15 -40 -10 -20 -7.5 -15];
known_max = [40 4 -10 15 40 40 40 7.5 15];
default_min = -15;
default_max = 15;

% Get a few key variables
range = f{'range'};
elev = f{'elev'};
azim = f{'azim'};
scantype = f.scantype(:);
location = f.location(:);

% If no scantype attribute, assume RHI
if isempty('scantype')
  scantype = 'RHI';
end
if isempty('location')
  location = 'Chilbolton'
end

% Count the number of 2D fields
names_2d = {};
n_2d = 0;
names = ncnames(var(f));
for ii = 1:length(names)
  the_size = size(f{names{ii}}(:));
  if length(find(the_size > 1)) > 1
    n_2d = n_2d + 1;
    names_2d{n_2d} = names{ii};
  end
end

% If less than 3 2D fields then have blank space at the foot of the page
n_subplots = n_2d;
if n_subplots < 3
  n_subplots = 3;
end

% Make an A4 figure
handle = figure('Units','inches','paperposition',[0.5 0.5 7.5 10.5], ...
		'position',[0.5 0.8 7.5 10.5],'Resize','off', ...
		'papertype','a4','paperorientation','portrait');

set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultTextFontSize',12);

% Generate the coordinate arrays
r_earth = 6371; % radius of the earth (km)
r = cos(elev.*(pi/180)) * range';
z = sin(elev.*(pi/180)) * range' + sqrt(r.^2 + r_earth.^2) - r_earth;
x = (sin(azim.*(pi/180)) .* cos(elev.*(pi/180))) * range';
y = (cos(azim.*(pi/180)) .* cos(elev.*(pi/180))) * range';

% Plot each 2D field in turn
for ii = 1:n_2d
  % Create a subplot
  subplot(n_subplots,1,ii);
  % Remove missing data
  missing_value = f{names_2d{ii}}.missing_value(:);
  if ~isempty(missing_value)
    var = f{names_2d{ii}}(:,:);
    var(find(var == missing_value)) = NaN;
  end
  if strcmp(scantype, 'PPI')
    pcolor(x,y,var)
    xlabel(['Distance east of ' location ' (km)']);
    ylabel(['Distance north of ' location ' (km)']);
    daspect([1 1 1]);
    aux_string = sprintf('Elevation %4.1f\\circ', mean(f{'elev'}(:)));
  elseif strcmp(scantype, 'Fixed');
    pcolor((f{'time'}(:)-f{'time'}(1))*60,range',var');
    xlabel('Time (minutes)');
    ylabel('Height (km)');
    ylim([0 12])
    aux_string = sprintf('Elevation %4.1f\\circ Azimuth %5.1f\\circ', ...
			 mean(f{'azim'}(:)), mean(f{'elev'}(:)));
  else % assume RHI
    pcolor(r,z,var)
    xlabel(['Distance from ' location ' (km)']);
    ylabel('Height (km)');
    ylim([0 12])
    aux_string = sprintf('Azimuth %5.1f\\circ', mean(f{'azim'}(:)));
  end
  shading flat
  % If the variable is recognised then set the colour axis
  caxis_min = default_min;
  caxis_max = default_max;
  for jj = 1:length(known_var)
    if strcmp(known_var{jj}, names_2d{ii})
      caxis_min = known_min(jj);
      caxis_max = known_max(jj);
    end
  end
  caxis([caxis_min caxis_max]);
  title(['\bf' f{names_2d{ii}}.long_name(:)]);
  set(gca,'Layer','top');
  
  h = colorbar;
  axes(h);
  set(title([names_2d{ii} ' (' f{names_2d{ii}}.units(:) ')']), ...
      'interpreter', 'none');
end

% Title of figure

% First the date string
date_num = datenum(f.year(:), f.month(:), f.day(:));
date_string = datestr(date_num, 'dd-mmm-yyyy');
date_string(find(date_string == '-')) = ' ';

% Then the time string
start_time = f{'time'}(1);
start_hour = floor(start_time);
start_min = floor(60.0 * (start_time-start_hour));
start_sec = round(3600.0 * (start_time-start_hour-start_min/60));
time_string = sprintf('%02d:%02d:%02d', start_hour, start_min, start_sec);

title_string = f.title(:);
if isempty(title_string)
  title_string = f.system(:);
end

axes('position', [0.04 0.96 0.9 0.025]);
set(gca,'Visible','off');
axis([0 1 0 1])
text(0, 0.5, ['{\bf' title_string '}  ' 10 ...
	      'Tape ' num2str(f.file(:),'%04d') ...
	      ' Raster ' num2str(f.raster(:)) ...
	      ' Scan ' num2str(f.scan(:)) ...
	      '  ' date_string ' ' time_string ' UTC' ...
	      '  ' aux_string]);

close(f);

