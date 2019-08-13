function PBFA_auto

% This function intend to find PBFA locations automatically

try; h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

fn = 'C:\Users\sumedhe\Desktop\PBFA_auto_data_test.txt';

% Write the header if not exit
if exist(fn,'file')
    fID = fopen(fn,'a+');
else
    fID = fopen(fn,'a+');
    fprintf(fID,'Index\ttpbfa\txpbfa\typbfa\tIpbfa\tNpbfa\n');
end

fprintf(fID,'\n');

settings = h.sen_set;
g = h.g;

tshift=settings.t_shift;

factor = g.factor;

ch_legend=cell(1,60);
for i=1:20;
    ch_legend{i*3-2}=[settings.sen_IDs{i} ':ch1'];
    ch_legend{i*3-1}=[settings.sen_IDs{i} ':ch2'];
    ch_legend{i*3}=[settings.sen_IDs{i} ':ch3'];
end


%% Load CGLSS data
cglss_fn = sprintf('%s/cglss/%s/%s/KSCCGLSS%s%s%s.dat',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},g.MM{:},g.DD{:});
x0 = 0;
y0 = 0;
z0 = 0;

% Load CGLSS data
cglss = CGLSS_extract(cglss_fn,g.t1,g.t2,str2double(settings.ldar_r),x0,y0,z0);

% Total number of CGLSS sources
tot_CG = length(cglss(:,1));

%% Load fast antenna data

TW = 1e-3;              % total time window to load
tw = 100e-6;            % Peak finding window size

[ch_fn h_fn] = generate_ch_fn(settings,g);

wbh = waitbar(0,sprintf('Loading Fast Antenna data. (%.2f%%)',0));
figure
hold all
tools2fig

ts = [];
sns = [];

nOfPks  = [];
peak_t  = [];
peak_v  = [];
delta_t = [];

for i=1:60
    
    try
        val = i/60;
        waitbar(val,wbh,sprintf('Loading Fast Antenna data. (%.2f%%)',val*100))
    catch
        return
    end
    
    
    if strcmp(ch_fn{i},'')==0
        
        [t y ch_freq_str] = FA_Extract1(ch_fn{i},h_fn{i},g.t1,g.t2,tshift(i),settings,i);
        %plot(t,y)
        
        % distance
        sn = ceil(i/3);
        
        if ~isempty(t)
            
            % Replace NaNs with zero so that filtfilt will work fine
            yn = y;
            yn(isnan(yn))=0;
            
            % filter out low frequencies
            Fs = 1/(t(2)-t(1));
            
            if isnan(Fs)
                Fs = 1/(t(5) - t(4));
            end
            
            
            [z,p] = butter(5,5000/(Fs/2),'high'); % Create a High-pass butterworth filter;
            
            yn = filtfilt(z,p,yn);    % filter the data.
            
            
            % Hilbert transform of voltages
            yn_hil = abs(hilbert(yn));            
            
            
            
            %save 1 data set for comparison
            if ~exist('y_comp','var')
                
                y_comp = yn;
                t_comp = t;
                dt = 0;
                
            else
                
                % finding convolution
                dt = t_comp(2) - t_comp(1);
                
                cy = xcorr(yn,y_comp);
                
                [mm mmn ] = max(cy);
                
                lag = median([1:length(cy)]) - mmn;
                
                dt = lag*dt;
                
            end
                    
            plot(t+dt,yn_hil)

            
            
            % Find peaks
            [ pksLocs,  pks] = peakfinder(yn_hil,0.03,0.03);
%           
            sns     = [sns sn];
            nOfPks  = [nOfPks length(pksLocs)];
            peak_t  = [peak_t t(pksLocs) + dt NaN];
            peak_v  = [peak_v pks NaN];
            delta_t = [delta_t dt];
            
            plot(t(pksLocs)+dt,pks,'o')
            
                       
            % Down sample to 0.5MHz to get peak e-field
            f=round(1/(t(2)-t(1)));
            % Number of intervals
            ni = round(f/0.5e6);
            
            if ni > 1
                t = downsample(t,ni);
                y = downsample(y,ni);
            end
            
            
            % Derivative of voltage (or the difference of the voltage)
            dy = y(2:end)-y(1:end-1);
            
            %plot(t,y*factor(i))
            %plot(t(1:end-1)+dt/2,dy*factor(i))
            
%             [maxdE ind] = max(abs(dy));
%             maxdE = maxdE*factor(i);
%             
%             plot(t(ind),yn(ind),'ro','MarkerFaceColor','r')
            maxdE_name= [ch_legend{i} ch_freq_str];
            
        end
    end
end

find_pbfa_positions(sns,nOfPks,peak_t,peak_v,delta_t,settings)


delete(wbh)

function [ch_fn h_fn] = generate_ch_fn(settings,g)
% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

% checking witch graphs are on
ch_on=g.chgraphs;

for i=1:60
    if ch_on(i)==1
        % finding the file extention number
        ext=mod(i,3);
        if ext==0
            ext=3;
        end
        
        % Finding the stattion ID
        sid=settings.sen_IDs{ceil(i/3)};
        
        % finding the file name
        filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        
        hfilename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        
        
        % If file is not exist don't store the file name
        if exist(filename,'file')==0 || exist(hfilename,'file')==0
            ch_fn{i}='';
            h_fn{i}='';
        else
            ch_fn{i}=filename;
            h_fn{i}=hfilename;
        end
    else
        ch_fn{i}='';
        h_fn{i}='';
    end
end

function find_pbfa_positions(sen_num,nOfPks,peak_t,peak_v,delta_t,settings)
%% find PBFA positions

% Let's find the longest record of peaks
[L maxInd] = max(nOfPks);
m = length(sen_num);

pk_t = NaN(L,m);
pk_v = pk_t;

start = 1;

x0 = settings.x0;
y0 = settings.y0;
z0 = settings.z0;

% Make line array to squire matrix for  easy search
for j = 1:m
    pk_t(1:nOfPks(j),j) = peak_t(start:start+nOfPks(j)-1);
    pk_v(1:nOfPks(j),j) = peak_v(start:start+nOfPks(j)-1);
    start = start+nOfPks(j)+1;
end



% time window +-tw
tw = 5.0e-6;

counter_start = 1;

arg.method = 2;
arg.outer = [];      % Outer sensors 
arg.t_out = [];

% [AX,H1,H2]=plotyy(nan,nan,nan,nan);
% hold(AX(2), 'on')
 

%figure
for i = 1:L
    
    temp = pk_t(i,maxInd);
    
    arg.inner = [];      % Inner sensors      
    arg.t_in  = [];
        
    for j=1:m
        [mini ind] = min(abs(pk_t(:,j)-temp));
        
        % is it within the time window
        if mini < tw
            arg.inner = [arg.inner sen_num(j)];
            arg.t_in  = [arg.t_in  (pk_t(ind,j) - delta_t(j))];
            pk_t(ind,j) =inf; 
        end        
    end
    
%     arg.inner'
%     arg.t_in'
%     
%     return
    
    n = length(arg.t_in);
    
    if n > 4
        [x1,y1,z1,t1]=pbfa_finder(arg);
        if isreal(z1) && z1 > 3000    && z1 < 15000 ...
                      && abs(x1) < 100000 && abs(y1) < 100000
            fprintf('%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',counter_start,t1,x1,y1,z1,n)
            delt = sqrt((x1-x0)^2+(y1-y0)^2+(z1-z0)^2)/3e8;
            plot(t1+delt-delta_t(2),0,'ro','MarkerFaceColor','r','MarkerSize',4)
            counter_start = counter_start + 1;
        end
    end      
end
%     
% function find_pbfa_positions(sen_num,nOfPks,peak_t,peak_v,delta_t,settings)
% %% find PBFA positions
% 
% % Let's find the longest record of peaks
% [L maxInd] = max(nOfPks);
% m = length(sen_num);
% 
% pk_t = NaN(L,m);
% pk_v = pk_t;
% 
% start = 1;
% 
% x0 = settings.x0;
% y0 = settings.y0;
% z0 = settings.z0;
% 
% % Make line array to squire matrix for  easy search
% for j = 1:m
%     pk_t(1:nOfPks(j),j) = peak_t(start:start+nOfPks(j)-1);
%     pk_v(1:nOfPks(j),j) = peak_v(start:start+nOfPks(j)-1);
%     start = start+nOfPks(j)+1;
% end
% 
% 
% 
% % time window +-tw
% tw = 5.0e-6;
% 
% counter_start = 1;
% 
% arg.method = 1;
% arg.outer = [];      % Outer sensors 
% arg.t_out = [];
% 
% % [AX,H1,H2]=plotyy(nan,nan,nan,nan);
% % hold(AX(2), 'on')
% 
% 
% % figure
% % plot(peak_t,peak_t,'ro')
% 
% peak_t = sort(peak_t);
% 
% time = peak_t(1);
% index = 0;
% 
% for i=2:length(peak_t);
%     if ~isnan(peak_t(i))
%         if (peak_t(i)-time(1)) < tw
%             time = [time peak_t(i)];
%         else
%             if length(time) > 4
%                 index = index + 1
%                 time'
%                 time = peak_t(i);
%             else
%                 time = peak_t(i);
%             end
%         end
%     end    
% end
% 
% 
% % %figure
% % for i = 1:L
% %     
% %     
% %     
% %     
% %     temp = pk_t(i,maxInd);
% %     
% %     arg.inner = [];      % Inner sensors      
% %     arg.t_in  = [];
% %         
% %     for j=1:m
% %         [mini ind] = min(abs(pk_t(:,j)-temp));
% %         
% %         % is it within the time window
% %         if mini < tw
% %             arg.inner = [arg.inner sen_num(j)];
% %             arg.t_in  = [arg.t_in  (pk_t(ind,j) - delta_t(j))];
% %             pk_t(ind,j) =inf; 
% %         end        
% %     end
% %     
% % %     arg.inner'
% % %     arg.t_in'
% % %     
% % %     return
% %     
% %     n = length(arg.t_in);
% %     
% %     if n > 4
% %         [x1,y1,z1,t1]=pbfa_finder(arg);
% %         if isreal(z1) && z1 > 3000    && z1 < 15000 ...
% %                       && abs(x1) < 100000 && abs(y1) < 100000
% %             fprintf('%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',counter_start,t1,x1,y1,z1,n)
% %             delt = sqrt((x1-x0)^2+(y1-y0)^2+(z1-z0)^2)/3e8;
% %             plot(t1+delt-delta_t(2),0,'ro','MarkerFaceColor','r','MarkerSize',4)
% %             counter_start = counter_start + 1;
% %         end
% %     end      
% end