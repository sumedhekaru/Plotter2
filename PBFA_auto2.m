function PBFA_auto2
%clc

try; h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

settings = h.sen_set;

t1 = 77636.333; % start time
t2 = 77637; % end

g = h.g;

TW = 5e-3;              % total time window to load

wbh1 = waitbar(0,sprintf('Loading data and calculating PBFA (%.2f%%)',0));


fn = 'C:\Users\sumedhe\Desktop\PBFA_auto_data_test.txt';
% Write the header if not exit
if exist(fn,'file')
    fID = fopen(fn,'a+');
else
    fID = fopen(fn,'a+');
    fprintf(fID,'Index\ttpbfa\txpbfa\typbfa\tIpbfa\tNpbfa\n');
end

fprintf(fID,'\n');

n = floor((t2-t1)/TW);

index = 1;

for i=1:n
    g.t1 = t1+(i-1)*TW;
    g.t2 = g.t1+TW;
    
    index = PBFA_auto1(settings,g,index,fID);

    try
        waitbar(i/n,wbh1,sprintf('Calculating PBFA %s UT   (%.2f%%)',...
            sec2hhmmss(g.t1),i/n*100));
    catch
        return
    end
    %return
end

delete(wbh1)
fclose(fID);

function index=PBFA_auto1(settings,g,index,fID)

% This function intend to find PBFA locations automatically

tshift=settings.t_shift;

factor = g.factor;

ch_legend=cell(1,60);
for i=1:20;
    ch_legend{i*3-2}=[settings.sen_IDs{i} ':ch1'];
    ch_legend{i*3-1}=[settings.sen_IDs{i} ':ch2'];
    ch_legend{i*3}=[settings.sen_IDs{i} ':ch3'];
end


% %% Load CGLSS data
% cglss_fn = sprintf('%s/cglss/%s/%s/KSCCGLSS%s%s%s.dat',...
%     settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},g.MM{:},g.DD{:});
% x0 = 0;
% y0 = 0;
% z0 = 0;
% 
% % Load CGLSS data
% cglss = CGLSS_extract(cglss_fn,g.t1,g.t2,str2double(settings.ldar_r),x0,y0,z0);
% 
% % Total number of CGLSS sources
% tot_CG = length(cglss(:,1));

%% Load fast antenna data

[ch_fn h_fn] = generate_ch_fn(settings,g);

% wbh = waitbar(0,sprintf('Loading Fast Antenna data. (%.2f%%)',0));


% figure
% hold all
% tools2fig

sns = [];

nOfPks  = [];
peak_t  = [];
peak_v  = [];
delta_t = [];

for i=1:60
    
%     try
%         val = i/60;
%         waitbar(val,wbh,sprintf('Loading Fast Antenna data. (%.2f%%)',val*100))
%     catch
%         return
%     end
    
    
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
            %yn_hil = abs(hilbert(yn));            
            
            
            % find the current frequency
            f=round(1/(t(2)-t(1)));
            % Number of intervals to sum 1MHz            
            ni = round(f/0.5e6);
            
                     
            if ni > 1
                len = floor(length(t)/ni);
                
                temp1 = nan(ni,len);
                temp2 = temp1;
                
                for indn = 1:ni
                    tempt= downsample(t,ni,indn-1);
                    tempy= downsample(yn,ni,indn-1);
                    temp1(indn,:)=tempt(1:len);
                    temp2(indn,:)=tempy(1:len);
                    clear tempt tempy
                    
                end
                
                tn1=mean(temp1);
                yn1=mean(temp2);
            end
            
            
            %save 1 data set for comparison
            if ~exist('y_comp','var') && exist('yn1','var')
                
                y_comp = yn1;
                t_comp = tn1;
                dt = 0;
                
                d.t{1} = tn1;
                d.y{1} = yn1*factor(i);
                data_cnt = 1;
                d.sn{1} = sn;
                
                
            elseif exist('yn1','var')
                
                % finding convolution
                dt = t_comp(2) - t_comp(1);
                
                cy = xcorr(yn1,y_comp);
                
                [mm mmn ] = max(cy);
                
                lag = median([1:length(cy)]) - mmn;
                
                dt = lag*dt;
                
                data_cnt = data_cnt + 1;
                d.t{data_cnt} = tn1;
                d.y{data_cnt}  = yn1*factor(i);
                d.sn{data_cnt} = sn;
                
                
            else
                
                dt = 1; % 
                
            end
            
                    
            %plot(t+dt,yn_hil)
            %plot(t+dt,yn)

            
            
            % Find peaks
            if dt < 240e-6 % time to travel 70km
                [ pksLocs,  pks] = peakfinder(yn,0.03,0.03,-1);
                
                sns     = [sns zeros(1,length(pksLocs))+sn];
                nOfPks  = [nOfPks length(pksLocs)];
                peak_t  = [peak_t t(pksLocs) + dt];
                peak_v  = [peak_v pks NaN];
                delta_t = [delta_t zeros(1,length(pksLocs))+dt];
                
                %plot(t(pksLocs)+dt,pks,'o')
            end
        end
    end
end

if nOfPks > 4
    index = find_pbfa_positions(sns,peak_t,delta_t,settings,index,d,fID);
end


% g.pbfa = 1;
% g.chgraphs = zeros(1,60);
% g.chgraphs(4) = 1;
% plot_all3(g)

%delete(wbh)


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

    
function index=find_pbfa_positions(sen_num,peak_t,delta_t,settings,index,d,fID)
%% find PBFA positions

% time window +-tw
tw = 2.5e-6;


% [AX,H1,H2]=plotyy(nan,nan,nan,nan);
% hold(AX(2), 'on')


% figure
% plot(peak_t,peak_t,'ro')

peaks = NaN(length(peak_t),3);

peaks(:,1) = peak_t';
peaks(:,2) = sen_num';
peaks(:,3) = delta_t';

peaks = sortrows(peaks,1);

peak_t = peaks(:,1);
sen_num = peaks(:,2);
delta_t = peaks(:,3);


time = peak_t(1);
ts = peak_t(1) - delta_t(1);
sns  = sen_num(1);
%index = 0;

for i=2:length(peak_t);
    if ~isnan(peak_t(i))
        if (peak_t(i)-time(1) ) < tw
            time = [time peak_t(i)];
            ts   = [ts (peak_t(i)-delta_t(i))];
            sns  = [sns sen_num(i)];
        else
            n = length(ts);
            
            if n > 4
                [ts sns] = clean_t(time,ts,sns);
            end
            
            n = length(ts);
            
            if n > 4
               
                                        
                index = calculate_PBFA(ts,sns,index,settings,d,fID);
                
                %return
                
                time = peak_t(i);
                sns  = sen_num(i);
                ts = peak_t(i) - delta_t(i);
            else
                time = peak_t(i);
                sns  = sen_num(i);
                ts = peak_t(i) - delta_t(i);
            end
        end
    end    
end


function index = calculate_PBFA(ts,sns,index,settings,d,fID)
 

arg.method = 2;
arg.outer = [];      % Outer sensors 
arg.t_out = [];

n = length(ts);

% xs = NaN(1,n);
% ys = NaN(1,n);
% zs = NaN(1,n);

% for i=1:n
%     t_temp = ts;
%     sns_temp= sns;
%     t_temp(i)= [];
%     sns_temp(i)=[];
    
   
    
    [xs,ys,zs,t1]=find_XYZ(ts,sns,settings);
    
    if isreal(zs) && zs > 3000    && zs < 18000 ...
            && abs(xs) < 100000 && abs(ys) < 100000
        
        Ip = Ip_cal(xs,ys,zs,t1,d,settings); 
        %plot(t1,0,'ro','MarkerFaceColor','r')
        fprintf('%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
            index,t1,xs(1),ys(1),zs(1),n,Ip)
        
         fprintf(fID,'%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
            index,t1,xs(1),ys(1),zs(1),n,Ip);
        
        index = index+1;
        
    end
% end


function [x1,y1,z1,t1]=find_XYZ(t,sns,sen_set)
%clc
% Speed of light in air
v = 299792458/1.0003;

% t'

% Creating K and g arrays
sns_tmp = sns;
t_tmp = t;

tmin    = min(t);
t = t - tmin;

x       = sen_set.x(sns);
y       = sen_set.y(sns);
z       = sen_set.z(sns);

t(6:end) = [];
x(6:end) = [];
y(6:end) = [];
z(6:end) = [];


% x'
% y'
% z'


for i=1:length(t)-1;
    K(i,1) = x(i+1) - x(i);
    K(i,2) = y(i+1) - y(i);
    K(i,3) = z(i+1) - z(i);
    K(i,4) = ( t(i+1) - t(i) );
    
    g(i,1) = 0.5*(x(i+1)^2 + y(i+1)^2 + z(i+1)^2 ....
                - x(i)^2 - y(i)^2 - z(i)^2 ...
                - v^2*( t(i+1)^2-t(i)^2 ));
end


%f=K\g;
% f = pinv(K)*g;
[f,flag,relres]= lsqr(K,g,eps*2,1000);

%flag

%relres

x1 = f(1);
y1 = f(2);
z1 = f(3);
t1 = -f(4)/v^2; 

% Find distance to each sensor
r = sqrt((x -x1).^2+ (y-y1).^2);

% Calculate z value
indx = find(r < 20000);

% remove BCC and FFI as they are 1MHz data (not good to calculate z)
indx(indx==5) = [];
indx(indx==11)= [];


% Closest sensor number
%csn = sns(ind);

z11 = sqrt(v^2*(t(indx)-t1).^2 - (x(indx)-x1).^2 - (y(indx) - y1).^2) - z(indx);
% z11'
% r'

t1 = t1 + tmin;
z1 = mean(z11);


function [ts sns] = clean_t(time,ts,sns)
% this function check for repeated sensors and remove those

for i=min(sns):max(sns)
    
    % Repititions
    rep = find(sns == i);
    
    if sum(rep) > 1
               
        dts = abs(time(rep) - mean(time));
        [m ind]     = min(dts);
        rep(ind)    =[];
        time(rep)   =[];
        ts(rep)     =[];
        sns(rep)    =[];

    end
end
    

function Ip = Ip_cal(xs,ys,zs,t1,d,settings)
% figure
% hold all

% +- time window
tw = 15e-6;
L = length(d.sn);
Ip = NaN(1,L);

for i=1:L
    
    sn = d.sn{i};
    r  = (((settings.x(sn)-xs)^2+(settings.y(sn)-ys)^2+(settings.z(sn)-zs)^2)^0.5);
    
    if r > 50000
        % Propagation time shift
        dtt = r/299792458;
        
        % choose data for current pulse
        lol = sum(d.t{i} < (t1 + dtt - tw));
        ul  = sum(d.t{i} < (t1 + dtt + tw));
        
        t = d.t{i}(lol:ul);
        y = d.y{i}(lol:ul);
        
        % Derivative of voltage (or the difference of the voltage)
        dy = abs(y(2:end)-y(1:end-1));
        % dt = t(2)-t(1);
        %plot(t,y*factor(i))
        %plot(t(1:end-1)+dt/2+dtt,dy)
        
        [maxdE ind] = max(dy);
        
        if maxdE > 0.5           
            Ip(i) = 0.0347*maxdE*(r/1000)^1.13+ 0.9568;
        end

%         
%         maxdE_name= [ch_legend{i} ch_freq_str];
        
        
       %plot(t+dtt,y)
    end
    
end

Ip = nanmean(Ip);



