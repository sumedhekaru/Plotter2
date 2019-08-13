function [tn,vn]=SA_Extract1(fn,t1,t2,tshift)

[data_i t_i] = read_lp_file(fn,tshift);

% Do we need to load data from Next or Previous 5 miniutes due to time
% shift?
if tshift ~= 0
    hh = str2double(fn(end-9:end-8));
    mm = str2double(fn(end-7:end-6));
    ss = str2double(fn(end-5:end-4));
    
    YYYY = str2double(fn(end-18:end-15));
    MM   = str2double(fn(end-14:end-13));
    DD   = str2double(fn(end-12:end-11));
    
    tmin = hh*3600+mm*60;
    tmax = tmin+300;
    
    
    if (t1 + tshift) > tmin && t1 < (tmin+tshift) 
        % need to load previous 5 mins of data
        % Previous file name
   
        now = datenum(YYYY,MM,DD,hh,mm,ss);
        new_t = datevec(now - 300/86400);
        
        % modify file name
        fn = sprintf('%s%4.4i/%2.2i/%2.2i/%s%4.4i%2.2i%2.2i_%2.2i%2.2i00%s', ...
            fn(1:end-34),new_t(1),new_t(2),new_t(3),fn(end-23:end-19), ...
            new_t(1),new_t(2),new_t(3),new_t(4),new_t(5),fn(end-3:end));
         
        if exist(fn,'file')
        % Load previous 5 mins data
            [temp_y temp_t] = read_lp_file(fn,tshift);
        
            data_i = [temp_y ; data_i];
            t_i    = [temp_t  t_i];
        else
            fprintf('Prev 5min file not found\n\t%s\n',fn)
        end
        
        %disp('Previous file is loaded')
    end

        
        
    if (t2 + tshift) < tmax && t2 > tmax+tshift
        
        % Need to load next 5 mins data
                % Previous file name
   
        now = datenum(YYYY,MM,DD,hh,mm,ss);
        new_t = datevec(now + 300/86400);
        
        % modify file name
        fn = sprintf('%s%4.4i/%2.2i/%2.2i/%s%4.4i%2.2i%2.2i_%2.2i%2.2i00%s', ...
            fn(1:end-34),new_t(1),new_t(2),new_t(3),fn(end-23:end-19), ...
            new_t(1),new_t(2),new_t(3),new_t(4),new_t(5),fn(end-3:end));
                 
        % Load previous 5 mins data
        if exist(fn,'file')
            [temp_y temp_t] = read_lp_file(fn,tshift);
            
            data_i = [data_i ; temp_y];
            t_i    = [t_i  temp_t];
        else
            fprintf('Prev 5min file not found\n\t%s\n',fn)
        end
        
        %disp('Next file is loaded')
    end
end

lol = nnz(t_i<t1)+1;           % Index of the lower matrix element
ul  = nnz(t_i<=t2);            % Index of the upper matrix element  
vn=data_i(lol:ul);             % data for given time interval
tn=t_i(lol:ul);


function [data_i t_i] = read_lp_file(fn,tshift)

% Read the header:
sHdr = epp_read_header( fn );

% What is the data point spacing (in sec)?
dT = 1. / (sHdr.sCont.freq);

% Read past header & read all data:
fId = fopen(fn, 'r');
discard = fread(fId, sHdr.masterSize, 'uint8');
data_i  = fread(fId, inf, 'float');
fclose( fId );

% Generate corresponding times relative to file start
%  (for absolute times, add 'sHdr.sCont.t0')
t_i = (0:dT:dT*(length(data_i)-1))+sHdr.sCont.t0+tshift;

