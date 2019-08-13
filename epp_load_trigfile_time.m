function [ sEpp ] = epp_load_trigfile_time( filename, daysec )
% SYNOPSIS: Loads the trigger containing the specified time (if any)
%
% PURPOSE: Attempts to load a trigger at the specified time in seconds
%          since midnight (UT) within the specified file.  If no trigger
%          is found, t_i and y_i will have no elements.  If no time is 
%          specified, then all events are returned.  This program
%          should be compatible with all triggered data files acquired by
%          the "epp2" data program (LASA2/CMCN/etc).
%
% USAGE:  sEpp = epp_load_trigfile_time( filename, daysec );
%
%   INPUTS:  filename      The path and name of the epp2 data file
%            daysec        The target time in day-seconds [sec]
%
%   OUTPUT:  sEpp          EPP structure with elements:
%               .t_i       Times of data in day-seconds [sec]
%               .y_i       Data voltages [V]
%               .sHdr      Header structure from file (see
%                           "epp_read_header" for more information)
%
% LIMITATIONS:
%  1) Can not deal with merged files only.  Requires a header file in
%     addition to the data file.  (OK, since only LANL may be using
%     merged files when this function was created)
%
% MODIFICATION HISTORY:
% '09May26     Created by Mark Stanley


%%%%  Initialize  %%%%

sEpp.t_i = [];
sEpp.y_i = [];
yRaw_i   = [];         % Raw digital values

% Exit if user is stupid
if (nargin < 1)            
    fprintf('ERROR (epp_load_trigfile_time): No file specified!');
    return
end


%%%%  Process file header  %%%%

[sHdr, ok] = epp_read_header( filename );
if (ok ~= 1)
    fprintf('ERROR [epp_load_trigfile_time]: Unable to process %s\n', ...
            filename);
    return
end;
sEpp.sHdr = sHdr;

% Verify that this is a triggered data file (not low-pass continuous)
if (sHdr.fileType ~= sHdr.TRIG)
    fprintf(['ERROR [epp_load_trigfile_time]: File %s is not a ' ...
            'triggered data file!\n'], filename);
    return
end;


%%%%  Determine key data parameters  %%%%

% Sample interval:
dT = 1. / sEpp.sHdr.sTrig.freq;

% Trigger length in time:
trigLen = sEpp.sHdr.sTrig.pre + sEpp.sHdr.sTrig.post;
tLength = trigLen * dT;


%%%%  Load header file  %%%%

[filePath fileText] = fileparts( filename );

hdrFilename = [filePath filesep fileText '.h' ...
               sprintf('%i', sEpp.sHdr.sChan.number)];
fHdr = fopen( hdrFilename, 'r');
if (fHdr == -1) 
    fprintf(['ERROR [epp_load_trigfile_time]: Unable to open header '...
            'file %s\n'], hdrFilename);
    return;
end
tH0_i = fread( fHdr, inf, 'double');
tH1_i = tH0_i + tLength - dT;
fclose( fHdr );


%%%%  Load trigger(s) from data file  %%%%

% Open data file (note: if we got this far, the file must exist)
%  and skip past main header
fData = fopen( filename, 'r' );
fread( fData, sEpp.sHdr.masterSize, 'int8' );

% Load all data if no 2nd argument:
if (nargin == 1)
    for i=1:numel(tH0_i)
        % Discard (redundant) header if this is a merged file
        if (sEpp.sHdr.sTrig.merged == 1)
            fread( fData, 1, 'double');
        end
        
        % Read digital data from trigger and append
        yTrigRaw_i = fread( fData, trigLen, 'uint16');
        yRaw_i = [yRaw_i, yTrigRaw_i'];      %#ok<AGROW>
        
        % Determine times of data points & store
        tTrig_i  = tH0_i(i):dT:tH1_i(i);
        sEpp.t_i = [sEpp.t_i, tTrig_i];
    end
else
    % Locate trigger of interest (if any)
    iTrig = find( (tH0_i <= daysec) & (tH1_i >= daysec) );
    
    if ( numel(iTrig) > 0 )
        % Size of each trigger in bytes?
        trigBytes = trigLen * 2;
        if (sEpp.sHdr.sTrig.merged == 1)
            % Skip past time-header in merged file
            trigBytes = trigBytes + 8;
        end

        % Fast forward to trigger
        fseek( fData, (iTrig(1) - 1) * trigBytes, 'cof' );
        
        % Discard (redundant) header if this is a merged file
        if (sEpp.sHdr.sTrig.merged == 1)
            fread( fData, 1, 'double');
        end
        
        % Read digital data from trigger and append
        yRaw_i = fread( fData, trigLen, 'uint16');
        yRaw_i = yRaw_i';
        
        % Determine times of data points & store
        sEpp.t_i  = tH0_i(iTrig(1)):dT:tH1_i(iTrig(1));    
    end
end

fclose( fData );


%%%%  Convert raw data to volts  %%%%

voltRange = 2. * sEpp.sHdr.sChan.range;
maxVal = (sEpp.sHdr.sChan.sum * 2^sEpp.sHdr.sDaq.bitDepth) - 1;
midVal = maxVal / 2.;
sEpp.y_i = (voltRange / maxVal) * (yRaw_i - midVal); 


end

