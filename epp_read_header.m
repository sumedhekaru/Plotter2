function [sHdr, ok] = epp_read_header(filename)
% epp_read_header: returns the master header of an EPP data file
%  
%  COMMENTS: This returns the contents of the master header in the form of 
%            a structure (of structures).  This routine works on both
%            triggered data files (ch* or m* extensions) and continuous 
%            data files (lp* extension).
%
%  USAGE:  [sHdr, ok] = epp_read_header(filename)
%
%    INPUT:   filename       The name (and path) of the data file
%
%    OUTPUT:  sHdr           Master header structure:
%                .masterSize    Size of master header in bytes
%                .generation    EPP Program generation (2)
%                .version       Version of the EPP program
%                .fileType      File type:  1: triggered  2: continuous
%                .sDaq        DAQ board structure:
%                   .manuf      Manufacturer (2 char code)
%                   .model      Board model number
%                   .count      Board count (usually 0 if 1st or only)
%                   .bitDepth   Bit-depth of A/D
%                .ppsSep      Raw sampling frequency when file opened [Hz]
%                               (before summation, if any)
%                .sGps        GPS structure:
%                   .abbrev     Station abbreviation (up to 6 char)
%                   .id         ID of GPS (see documentation)
%                   .lat        Latitude of station
%                   .lon        Longitude...
%                   .z          Altitude... [m]
%                   .numPts     Number of GPS points used for lat/lon/z
%                .sDate       Date structure:
%                   .year       Year (4-digit)
%                   .month      Month (1-12)
%                   .day        Day of month (1-31)
%                .sTime       Time structure:
%                   .hour       Hour (0-23)
%                   .min        Minute (0-59) at start of file
%                   .sec        Second (0-59) at start of file
%                   .daySec     Seconds since midnight (UT) for file start
%                .sChan       Channel structure:
%                   .number     Channel number
%                   .range      Voltage range (half of total -or- max V)
%                   .major      Major data code (see documentation)
%                   .minor      Minor...
%                   .sum        No. points summed into each data (1-16)
%                .sTrig       Triggered structure (valid for trig files):
%                   .pre        Pretrigger length in samples
%                   .post       Postrigger...
%                   .mode       Mode of triggering (bitwise and of:
%                                 1: fixed  2: low-freq float  4: slave
%                   .vNeg       Fixed-level negative voltage threshold
%                   .vPos       ...positive...
%                   .dVneg      Floating-level relative voltage threshold
%                   .dVpos      ...positive...
%                   .master     Master channel for slave triggering
%                   .merged     Are headers merged with data in file?
%                   .freq       Actual (final) sampling frequency [Hz]
%               .sCont        Continuous structure (valid for cont files):
%                   .mode       Mode of continuous sampling (0: average,
%                                 1: 2nd-order Butterworth filter)
%                   .fCutoff    Cutoff frequency (if Butterworth filter)
%                   .numAvg     Number of points in average
%                   .fTarget    Target sampling frequency
%                   .freq       Actual...
%                   .t0         Time of first data point (day-sec)
%               .TRIG         Constant (1) used to denote trigger data file
%               .CONT         ...(2)...continuous data file
%
%              ok        True(1) if read was OK, False(0) otherwise
%
%
%  MODIFICATION HISTORY:
%  '09Apr02    Created by Mark Stanley
%     May26    Fixed bug where 'sHdr.sTrig.merged' was not set for
%               older data versions


%%%%  Constants  %%%%

TRIG = 1;      % Triggered data file-type
CONT = 2;      % Continuous...

sHdr.TRIG = TRIG;
sHdr.CONT = CONT;


%%%% Open input file (exit if unable to) %%%%

fId = fopen(filename, 'r');
if (fId == -1) 
    fprintf(1, 'ERROR (epp_read_header): Unable to open file %s\n', ...
            filename);
    sHdr.masterSize = -1;
    ok = 0;
    return;
end;


%%%%  Critical information at start of header  %%%%

sHdr.masterSize = fread(fId, 1, 'uint16');
sHdr.generation = fread(fId, 1, 'uint8');
sHdr.version    = fread(fId, 1, 'float');

if (sHdr.version >= 1.38)
    sHdr.fileType = fread(fId, 1, 'uint8');
else
    % Only triggered data before v1.38
    sHdr.fileType = TRIG;
end;

% Bug patch (rare to encounter since was during development phase)
if ((sHdr.fileType == CONT) && (sHdr.version <= 1.385))
    sHdr.masterSize = sHdr.masterSize + 2;
end;


%%%%  DAQ board stuff  %%%%

if (sHdr.version >= 1.5)
    sHdr.sDaq.manuf = fread(fId, 4, 'uint8=>char')';
end;

% Board model/count/bit-depth important for 1.3 and beyond
if (sHdr.version > 1.30)
    if (sHdr.version < 1.31)
        sHdr.sDaq.model = fread(fId, 1, 'int32');
    else
        sHdr.sDaq.model = fread(fId, 1, 'int16');
    end;
    
    if (sHdr.version >= 1.35)
        sHdr.sDaq.count = fread(fId, 1, 'uint8');
    else
        sHdr.sDaq.count = 0;
    end;
    
    if (sHdr.version >= 1.31)
        sHdr.sDaq.bitDepth = fread(fId, 1, 'uint8');
    else
        sHdr.sDaq.bitDepth = 16;
    end;
else
    % epp2 only compatible with MC4020 before v1.30
    sHdr.sDaq.model    = 4020;
    sHdr.sDaq.bitDepth = 12;
end;

% Before Comedi (v1.5-), only 2 possibilities for manuf...
if (sHdr.version < 1.5)
    if (sHdr.sDaq.model == 4020)
        sHdr.sDaq.manuf = 'MC';
    else
        sHdr.sDaq.manuf = 'NI';
    end;
end;
    
% The raw sample frequency (PPS spacing) when file written:
if ((sHdr.version >= 1.5) && (sHdr.version <= 1.504) && ...
        (sHdr.fileType == CONT))
    % Bug where 64-bit long briefly used (64-bit machine)
    sHdr.ppsSep = fread(fId, 1, 'uint64');
else
    sHdr.ppsSep = fread(fId, 1, 'uint32');
end;

% This is only stored for triggered data:
if (sHdr.fileType == TRIG)
    sHdr.sTrig.ppsSepTarget = fread(fId, 1, 'uint32');
end;


%%%%  Station information  %%%%

sHdr.sGps.abbrev = fread(fId, 6, 'uint8=>char')';
sHdr.sGps.id     = fread(fId, 1, 'uint8');
sHdr.sGps.lat    = fread(fId, 1, 'double');
sHdr.sGps.lon    = fread(fId, 1, 'double');
sHdr.sGps.z      = fread(fId, 1, 'double');
sHdr.sGps.numPts = fread(fId, 1, 'uint32');


%%%%  Date/time of file  %%%%

sHdr.sDate.year  = fread(fId, 1, 'uint16');
sHdr.sDate.month = fread(fId, 1, 'uint8');
sHdr.sDate.day   = fread(fId, 1, 'uint8');

sHdr.sTime.hour   = fread(fId, 1, 'uint8');
sHdr.sTime.min    = fread(fId, 1, 'uint8');
sHdr.sTime.sec    = fread(fId, 1, 'uint8');
sHdr.sTime.daySec = fread(fId, 1, 'int32');

% Continuous data: time of first data point
if (sHdr.fileType == CONT)
    dT0 = fread(fId, 1, 'float');
    sHdr.sCont.t0 = sHdr.sTime.daySec + dT0;
end;


%%%%  Data channel setup  %%%%

sHdr.sChan.number = fread(fId, 1, 'uint8');

% Voltage range (max V) became a float when NI driver added:
if (sHdr.version >= 1.3)
    sHdr.sChan.range = fread(fId, 1, 'float');
else
    sHdr.sChan.range = fread(fId, 1, 'uint8');
end;

% Major and minor codes for data type
sHdr.sChan.major = fread(fId, 1, 'uint16');
sHdr.sChan.minor = fread(fId, 1, 'uint16');


%%%%  Triggered data setup (mostly)  %%%%

if (sHdr.fileType == TRIG)
    sHdr.sTrig.pre  = fread(fId, 1, 'uint32');
    sHdr.sTrig.post = fread(fId, 1, 'uint32');
    sHdr.sChan.sum  = fread(fId, 1, 'uint16');   % Channel setup
    sHdr.sTrig.mode = fread(fId, 1, 'uint16');
    
    % Thresholds in terms of voltages starting with v1.37
    if (sHdr.version >= 1.37)
        sHdr.sTrig.vNeg  = fread(fId, 1, 'float');
        sHdr.sTrig.vPos  = fread(fId, 1, 'float');
        sHdr.sTrig.dVneg = fread(fId, 1, 'float');
        sHdr.sTrig.dVpos = fread(fId, 1, 'float');
    else
        % Earlier versions used a single pair of digital threshold
        %  values regardless of mode
        digNeg = fread(fId, 1, 'uint16');
        digPos = fread(fId, 1, 'uint16');
        
        % Conversion factors:
        trueDepth = (2^sHdr.sDaq.bitDepth) * sHdr.sChan.sum;
        digMid    = trueDepth / 2.;
        dig2volt  = 2. * sHdr.sChan.range / trueDepth;
        
        if (sHdr.sTrig.mode && 1)
            % Fixed threshold mode:
            sHdr.sTrig.vNeg = dig2volt * (digNeg - digMid);
            sHdr.sTrig.vPos = dig2volt * (digPos - digMid);
        end;
        if (sHdr.sTrig.mode && 2)
            % Floating-DC threshold mode (always symmetric)
            dDig = (digPos - digNeg) / 2.;
            sHdr.sTrig.dVneg = -dDig * dig2volt;
            sHdr.sTrig.dVpos =  dDig * dig2volt;
        end;
    end;
        
    sHdr.sTrig.master = fread(fId, 1, 'uint8');
    if (sHdr.version >= 1.393)
        sHdr.sTrig.merged = fread(fId, 1, 'uint16');
        fread(fId, 1, 'uint16');
    else
        sHdr.sTrig.merged = 0;
    end;

    % Sampling frequency:
    sum = sHdr.sChan.sum;
    if (sum == 0)
        sum = 1;     % Sum should never be zero, but just in case...
    end;
    sHdr.sTrig.freq = sHdr.ppsSep / sum;
end;


%%%%  Continuous data setup (mostly)  %%%%

if (sHdr.fileType == CONT)
    sHdr.sChan.sum     = fread(fId, 1, 'uint16');    % Channel setup
    sHdr.sCont.mode    = fread(fId, 1, 'uint8');
    sHdr.sCont.fCutoff = fread(fId, 1, 'float');
    sHdr.sCont.numAvg  = fread(fId, 1, 'uint32');
    sHdr.sCont.fTarget = fread(fId, 1, 'float');
    sHdr.sCont.freq    = fread(fId, 1, 'float');

    % Development version bug fix:
    if (sHdr.version <= 1.385)
        dT0 = dT0 / sHdr.sCont.numAvg;
        sHdr.sCont.t0 = sHdr.sTime.daySec + dT0;
    end;
end;
    

%%%%  Cleanup  %%%%

% If we got this far everything should be OK    
ok = 1;

fclose(fId);