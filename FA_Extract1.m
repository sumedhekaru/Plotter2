function [t,y,ch_freq_str]=FA_Extract1(fn, hfn ,t1,t2,tshift,settings,sn)

[t,y,ch_freq_str]= read_ch_file(fn, hfn ,t1,t2,tshift,settings,sn);


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
        fn = sprintf('%s%4.4i/%2.2i/%2.2i%s%4.4i%2.2i%2.2i_%2.2i%2.2i00%s', ...
            fn(1:end-34),new_t(1),new_t(2),new_t(3),fn(end-23:end-19), ...
            new_t(1),new_t(2),new_t(3),new_t(4),new_t(5),fn(end-3:end));
        
        hfn = sprintf('%s%4.4i/%2.2i/%2.2i%s%4.4i%2.2i%2.2i_%2.2i%2.2i00%s', ...
            fn(1:end-34),new_t(1),new_t(2),new_t(3),fn(end-23:end-19), ...
            new_t(1),new_t(2),new_t(3),new_t(4),new_t(5),hfn(end-2:end));
        
        
        % Load previous 5 mins data if files are exsist
        if exist(fn,'file') && exist(hfn,'file')
            [temp_t temp_y ch_freq_str ] = read_ch_file(fn, hfn ,t1,t2,tshift,settings,sn);
        
            y = [temp_y NaN y];
            t = [temp_t NaN t];
        else
            fprintf('Prev 5min files not found\n\t%s\n\t%s\n',fn,hfn)
        end
        
        %disp('Previous file is loaded')
    end

        
        
    if (t2 + tshift) < tmax && t2 > tmax+tshift
        
        % Need to load next 5 mins data
                % Previous file name
   
        now = datenum(YYYY,MM,DD,hh,mm,ss);
        new_t = datevec(now + 300/86400);
        
        % modify file name
        fn = sprintf('%s%4.4i/%2.2i/%2.2i%s%4.4i%2.2i%2.2i_%2.2i%2.2i00%s', ...
            fn(1:end-34),new_t(1),new_t(2),new_t(3),fn(end-23:end-19), ...
            new_t(1),new_t(2),new_t(3),new_t(4),new_t(5),fn(end-3:end));
                 
       hfn = sprintf('%s%4.4i/%2.2i/%2.2i%s%4.4i%2.2i%2.2i_%2.2i%2.2i00%s', ...
            fn(1:end-34),new_t(1),new_t(2),new_t(3),fn(end-23:end-19), ...
            new_t(1),new_t(2),new_t(3),new_t(4),new_t(5),hfn(end-2:end));
                 
        % Load previous 5 mins data
        if exist(fn,'file') && exist(hfn,'file')        
            [temp_t temp_y ch_freq_str ] = read_ch_file(fn, hfn ,t1,t2,tshift,settings,sn);
        
             y = [y NaN temp_y];
             t = [t NaN temp_t];
        else
            fprintf('Next 5min files not found\n\t%s\n\t%s\n',fn,hfn)
        end
        
        %disp('Next file is loaded')
    end
end



% Finding data in given time range
lol=nnz(t<t1)+1;       % Index of the lower matrix element
ul=nnz(t<=t2) ;                    % Index of the upper matrix element
y=y(lol:ul);
t=t(lol:ul);



function [t,y,ch_freq_str]= read_ch_file(fn, hfn ,t1,t2,tshift,settings,sn)

ch_freq_str = '';


t=[];   % Variable for time
y=[];   % Variable for voltage
% Plotting up to 60 fast antenna files

% Load corresponding header file times
fId = fopen(hfn, 'r');
trigs  = fread(fId, inf, 'double') ;
fclose( fId );

% Introduce time shift before filtering time range
trigs = trigs + tshift;
% Finding the triggers between the given time range
lol=length(trigs)- nnz(trigs>t1)+1;       % Index of the lower matrix element
ul=nnz(trigs<t2);                      % Index of the upper matrix element
% Let's find two more triggers from both ends

if lol>1 ;    lol=lol-1; end

if ul < length(trigs) ;   ul=ul+1; end

% Triggers between given time range and remove time shift
trigs=trigs(lol:ul)-tshift;

% Loading fast antenna data
for j=1:length(trigs)
    sFa = epp_load_trigfile_time(fn,trigs(j));
    %%%%%%%%%%%%%%%% Reducing Sampling Rate %%%%%%%%%%%%%%%%%%%%%
    %find the current frequency
    
    if ~isempty(sFa.t_i)
        
        f=round(1/(sFa.t_i(2)-sFa.t_i(1)));
              
        
        switch settings.chSumMethod
            case 2 % User needs downsampling
                % Number of intervals
                ni = round(f/settings.chfreq(sn));
                
                if ni > 1
                    sFa.t_i = downsample(sFa.t_i,ni);
                    sFa.y_i = downsample(sFa.y_i,ni);
                    ch_freq_str = sprintf(':D-%.2fMHz',f/1e6/ni);
                else
                    ch_freq_str = sprintf(':O-%.2fMHz',f/1e6);
                end
                
                
                
            case 1 % User needs summing
                
                % Number of intervals
                ni = round(f/settings.chfreq(sn));
                
                if ni > 1
                    len = floor(length(sFa.t_i)/ni);
                    
                    temp1 = nan(ni,len);
                    temp2 = temp1;
                    
                    for indn = 1:ni
                        tempt= downsample(sFa.t_i,ni,indn-1);
                        tempy= downsample(sFa.y_i,ni,indn-1);
                        temp1(indn,:)=tempt(1:len);
                        temp2(indn,:)=tempy(1:len);
                        clear tempt tempy
                        
                    end
                    
                    sFa.t_i=mean(temp1);
                    sFa.y_i=mean(temp2);
                    ch_freq_str = sprintf(':S-%.2fMHz',f/1e6/ni);
                else
                    ch_freq_str = sprintf(':O-%.2fMHz',f/1e6);
                end
                
                
                
            otherwise
                % User needs nothing
                ch_freq_str = sprintf(':O-%.2fMHz',f/1e6);
                
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y=[y,NaN,sFa.y_i];
    t=[t,NaN,sFa.t_i];   

end

t = t+tshift;