function var_list = load_nc(nc_file, prefix, destination)
% load_nc -- Load NetCDF variables and global attributes.
%
% load_nc('nc_file') loads all variables and global
%   attributes of 'nc_file' into the Matlab workspace of the 
%   caller of this routine. The names of the loaded variables
%   and attributes are listed on the screen, and returned or
%   assigned to "ans". If any variable has a "missing_value"
%   attribute then those values in the corresponding variable
%   are set to NaN.
%
% load_nc('nc_file', 'prefix') loads all variables and global
%   attributes, but prefixes the name of each with "prefix".
%   This is useful if you load radar and lidar data at the same
%   time as both have variables called time, range etc.
%
% load_nc('nc_file', 'prefix', destination) operates as above, but
%   sends the variables to destination, which can be 'caller' (the
%   default), or 'base'.
%
%   To use this function you need the NetCDF Toolbox for
%   Matlab-5. This can be downloaded free from:
%     http://crusty.er.usgs.gov/~cdenham/MexCDF/nc4ml5.html
 
if nargin < 1, help(mfilename), return, end

if nargin < 2, prefix=''; end

if nargin < 3, destination='caller'; end

result = [];
if nargout > 0, var_list = result; end

f = netcdf(nc_file, 'nowrite');
if isempty(f), return, end
disp([name(f) ':']);

% Write variables to caller's workspace
names = ncnames(var(f));
for ii = 1:length(names)
   newname = names{ii};
   newname(find(newname == '-')) = '_';
   missing_value = f{names{ii}}.missing_value(:);
   if ~isempty(missing_value)
     if ischar(missing_value)
       missing_value = str2num(missing_value);
     end
      tmp_var = f{names{ii}}(:);
      tmp_var(find(tmp_var == missing_value)) = NaN;
      assignin(destination, [prefix newname], tmp_var);
      
      clear tmp_var      
   else
      assignin(destination, [prefix newname], f{names{ii}}(:))
   end
   the_size = size(f{names{ii}}(:));
   long_name = clean_up_string(f{names{ii}}.long_name(:));
   if ~isempty(long_name)
     long_name = clean_up_string(long_name);
     long_name = [' (' long_name ')'];
   end
   units = clean_up_string(f{names{ii}}.units(:));

   names{ii} = newname;

   namefill = blanks(max(0,14-length(newname)));

   if the_size(1) == 1 & the_size(2) == 1
     disp([namefill prefix newname ':' 9 num2str(f{names{ii}}(1)) ' ' units 9 9 ...
	   long_name]);
   elseif length(the_size) > 2
      disp([namefill prefix newname ':' 9 '[' num2str(the_size(1)) 'x' ...
	    num2str(the_size(2)) 'x' num2str(the_size(3)) '] ' units ...
	    9 long_name]); 
   else
      disp([namefill prefix newname ':' 9 '[' num2str(the_size(1)) 'x' ...
	    num2str(the_size(2)) '] ' units 9 long_name]); 
   end
end

% Write global attributes to caller's workspace
attnames = ncnames(att(f));
for ii = 1:length(attnames)
   eval(['attribute = f.' attnames{ii} '(:);']);
   assignin(destination, [prefix attnames{ii}], attribute);
   namefill = blanks(max(0,14-length(attnames{ii})));
   disp([namefill prefix attnames{ii} ': ' clean_up_string(attribute)]);
end

result = [names attnames];

close(f)

if nargout > 0
   var_list = result;
else
   ncans(result)
end

function newstr = clean_up_string(oldstr)
newstr = num2str(oldstr);
if length(newstr) > 1
  if newstr(end-1) == '\' & newstr(end) == '0'
    newstr = deblank(newstr(1:end-2));
  end
end
