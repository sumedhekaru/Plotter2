%% Open netcdf file

filename = 'C:\Users\Sumedhe\Desktop\RADAR\KMLB_N0R_20110814_210100.nc';

ncid = netcdf.open(filename,'NC_NOWRITE');



%% Explore the Contents

[numdims,nvars,natts] = netcdf.inq(ncid);



%% Get Global attributes Information

for ii = 0:natts-1

    fieldname = netcdf.inqAttName(ncid, netcdf.getConstant('NC_GLOBAL'), ii);

    fileinfo.(fieldname) = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'), fieldname );

end

 

% allocate structure

dimension = repmat(struct('name', '', 'length', 0), numdims, 1);

 

for ii = 1:numdims

 [dimension(ii).name, dimension(ii).length] = netcdf.inqDim(ncid,ii-1);



 % padding name for table layout

    padlength   = min(0, length(dimension(ii).name));

    name_padded = [dimension(ii).name repmat(' ', padlength+1)];

 

    fprintf('%s\t\t%d\n', name_padded, dimension(ii).length);

end

 

%% Get the Data

for ii = 1:nvars

    [name, ~, ~, natts] = netcdf.inqVar(ncid,ii-1);

   

    % Get Variable Attributes

    tmpstruct = struct();

    for jj = 1:natts

        fieldname = netcdf.inqAttName(ncid, ii-1, jj-1);

        tmpstruct.(fieldname) = netcdf.getAtt(ncid, ii-1, fieldname(2:end) );

    end

 

    % Get raw data

    data = netcdf.getVar(ncid,ii-1);

   

    % Replace Missing Numbers (if necessary

    if (isfield(tmpstruct, 'missing_value') )

        data( data == tmpstruct.missing_value ) = NaN;

    end

   

    % Scale data (if necessary)

    if( isfield(tmpstruct, 'scale_factor') )

        data = double(data) * tmpstruct.scale_factor;

    end

   

    % Apply offset (if necessary)

    if( isfield(tmpstruct, 'add_offset') )

        data = data + tmpstruct.add_offset;

    end

   

    % Transpose data from column major to row major

    if( isnumeric(data) && ndims(data) > 2 )

        data = permute(data, [2 1 3:ndims(data)]);

    elseif ( isnumeric(data) && ndims(data) == 2 )

        data = data';

    end

% store attribute and data with appropriate name

    varinfoname = [name '_info'];

    assignin('caller', varinfoname, tmpstruct);

    assignin('caller', name, data);

end

 

%% Close File

netcdf.close(ncid);

 

%% Animate and Create AVI

%% get days since the year 0000.

%% The maximum time dimensions have to be found before the script is executed

%% This could be done by typing time_info in the MATLAB screen.

%% And the decision to add date_offset should be dependent on the time dimension.

%% The time is here divided by 24 since the file used had hours since jan 1 1.

date_offset = datenum('01-JAN-1');

 

animate_netcdf(lat, lon, data(:, :,1:365), datestr((time(1:365)/24)+date_offset), 0.5, 'radiation_1990.avi');

 

%% Clean Up temporary variables

clear ndims nvars natts ii jj tmpstruct idx ncid filename

clear fieldname name vartype dimids varinfoname data date_off