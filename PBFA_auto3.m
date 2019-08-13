function PBFA_auto3
%clc

try; h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

settings = h.sen_set;
g = h.g;
tshift=settings.t_shift;

t1 = 77636.333; % start time
t2 = 77637.333; % end
%t1 = 77636.333; % start time
%t2 = 77636.339; % end



tw = 50e-3;              % total time window to load
TW = 1;               % data load window




wbh1 = waitbar(0,sprintf('Loading data and calculating PBFA (%.2f%%)',0));


fn = 'C:\Users\sumedhe\Desktop\test.txt';
% Write the header if not exit
if exist(fn,'file')
    fID = fopen(fn,'a+');
else
    fID = fopen(fn,'a+');
    fprintf(fID,'Index\ttpbfa\txpbfa\typbfa\tIpbfa\tNpbfa\n');
end

index = 1;

fprintf(fID,'\n');

fID2 = fopen([fn(1:end-4) '_peaks.txt'],'a+');


n = floor(TW/tw);
n2 = floor((t2-t1)/TW);

g.t1 = t1-tw;
g.t2 = t1;



for k=1:n2
    
    [ch_fn h_fn] = generate_ch_fn(settings,g);
    
    tL = t1 + (k-1)*TW;
    tR = tL + TW;
    
    for i=1:60
        
        if strcmp(ch_fn{i},'')==0            
            [dd.t{i} dd.y{i} ~] = FA_Extract1(ch_fn{i},h_fn{i},tL,tR,tshift(i),settings,i);
        else
            dd.t{i}=[];
            dd.y{i}=[];            
        end
    end
     
       
    for m=1:n
        g.t1 = g.t1+tw;
        g.t2 = g.t2+tw;
        %fprintf('%.4f\n',g.t1)     
        index = PBFA_auto1(settings,g,index,fID,fID2,dd);
        
        try
            val = ((k-1)*n+m)/(n2*n);
            waitbar(val,wbh1,sprintf('Calculating PBFA %s UT   (%.2f%%)',...
                sec2hhmmss(g.t1),val*100));
        catch
            return
        end
    end
end

delete(wbh1)
fclose(fID);
fclose(fID2);

function index=PBFA_auto1(settings,g,index,fID,fID2,dd)

factor = g.factor;

% wbh = waitbar(0,sprintf('Loading Fast Antenna data. (%.2f%%)',0));


% figure
% hold all
% tools2fig

sns = [];
nOfPks  = [];
peak_t  = [];
peak_v  = [];
delta_t = [];

sns1 = [];
nOfPks1  = [];
peak_t1  = [];
peak_v1  = [];
delta_t1 = [];

for i=1:33
    
    
    if ~isempty(dd.t{i})
        
        lol = sum(dd.t{i} < g.t1) + 1;
        ul  = sum(dd.t{i} < g.t2);
        
        t = dd.t{i}(lol:ul);
        y = dd.y{i}(lol:ul);
        
%         plot(t,y)
        
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
            
            
            [z,p] = butter(5,2000/(Fs/2),'high'); % Create a High-pass butterworth filter;
            
            try
                yn = filtfilt(z,p,yn);    % filter the data.
            catch
                yn = [];
            end
            
            if ~isempty(yn)
                
                
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
                    try
                        dt = t_comp(2) - t_comp(1);
                    catch
                        dt = 0.2e-6;
                    end
                    
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
                    [ pksLocs,  pks] = peakfinder(yn,0.01,0.01,+1);
                    
                    sns     = [sns zeros(1,length(pksLocs))+sn];
                    nOfPks  = [nOfPks length(pksLocs)];
                    peak_t  = [peak_t t(pksLocs) + dt];
                    peak_v  = [peak_v pks NaN];
                    delta_t = [delta_t zeros(1,length(pksLocs))+dt];
                    
                    [ pksLocs,  pks] = peakfinder(yn,0.01,0.01,+1);
                    
                    sns1     = [sns1 zeros(1,length(pksLocs))+sn];
                    nOfPks1  = [nOfPks1 length(pksLocs)];
                    peak_t1  = [peak_t1 t(pksLocs) + dt];
                    peak_v1  = [peak_v1 pks NaN];
                    delta_t1 = [delta_t1 zeros(1,length(pksLocs))+dt];
                    
                    %plot(t(pksLocs)+dt,pks,'o')
                end
            end
        end
    end
end

if length(nOfPks) > 4
    index = find_pbfa_positions(sns,peak_t,delta_t,settings,index,d,fID,fID2);
end

if length(nOfPks1) > 4
    index = find_pbfa_positions(sns1,peak_t1,delta_t1,settings,index,d,fID,fID2);
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


hms = sec2hhmmss(g.t1);
hh = str2double(hms(1:2));
mm = floor(str2double(hms(4:5))/5)*5;

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
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},hh,mm,0,ext);
        
        hfilename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},hh,mm,0,ext);
        
        
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

    
function index=find_pbfa_positions(sen_num,peak_t,delta_t,settings,index,d,fID,fID2)
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

try    
    time = peak_t(1);
    ts = peak_t(1) - delta_t(1);
    sns  = sen_num(1);
catch
    return
end
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
               
                                        
                index = calculate_PBFA(ts,sns,index,settings,d,fID,fID2);
                
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


function index = calculate_PBFA(ts,sns,index,settings,d,fID,fID2)
 
v = 299792458/1.0003;
arg.method = 5;
arg.outer = [];      % Outer sensors 
arg.t_out = [];
arg.inner = sns;
arg.t_in = ts;



n = length(ts);



[xs,ys,zs,t1]=pbfa_finder(arg);

if isreal(zs)
    
    ki_sqrd = cal_ki_sqrd(t1,xs,ys,zs,arg,settings);
    
    if n > 5 && ki_sqrd < 5
        
        Ip = Ip_cal(xs,ys,zs,t1,d,settings);
        
        snsStr = sprintf('%i,',sort(sns));    
        
        fprintf('%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\tNaN\tNaN\tNaN\tNaN\t%0.1f\n',...
            index,t1,xs,ys,zs,n+50,Ip,snsStr(1:end-1),ki_sqrd)
        
        fprintf(fID,'%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\tNaN\tNaN\tNaN\tNaN\t%0.1f\n',...
            index,t1,xs,ys,zs,n+50,Ip,snsStr(1:end-1),ki_sqrd);
                
        store_peaks_info(index,[ts',sns'],fID2)
        
        index = index+1;
       
    end
end


% 
% x = settings.x;
% y = settings.y;
% z = settings.z;

%[xs,ys,zs,t1]=find_XYZ(ts,sns,settings);

% Let's check the point is valid

% if isreal(zs) && zs > 3000    && zs < 18000 ...
%         && abs(xs) < 100000 && abs(ys) < 100000
%     r =( (x(sns)-xs).^2+(y(sns)-ys).^2+(z(sns)-zs).^2).^.5;
%     [dt ind] = max(abs(ts-t1-r/v)*1e6);
%     
%     % if dt is more than a micro second, let's drop the closest sensor
% %     if dt > .50
% %         %[m ind] = min(r);
% %         ts(ind) = [];
% %         sns(ind) = [];
% %         
% %         if length(ts) > 4
% %             fprintf('\n')       
% %             index = calculate_PBFA(ts,sns,index,settings,d,fID);
% %             return
% %         end
% %     else
%                 
%         Ip = Ip_cal(xs,ys,zs,t1,d,settings);
%         
%         %plot(t1,0,'ro','MarkerFaceColor','r')
%         fprintf('%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\t\t%0.2f\n',...
%             index,t1,xs(1),ys(1),zs(1),n,Ip,dt)
%         
%         fprintf(fID,'%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
%             index,t1,xs(1),ys(1),zs(1),n,Ip);
%         
%         index = index+1;
% %     end
%     
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

% t(6:end) = [];
% x(6:end) = [];
% y(6:end) = [];
% z(6:end) = [];


% x'
% y'
% z'
L = length(t);
K = NaN(L-1,4);
g = NaN(L-1,1);

for i=1:L-1;
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

try
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
end
Ip = nanmean(Ip);

function ki_sqrd = cal_ki_sqrd(t,x,y,z,arg,sen_set)

        % Ki -sqrd
        tt1 = arg.t_in';
        
        deg_free = length(arg.t_in);
        
        tt2 = t+sqrt((x - sen_set.x(arg.inner)).^2 + ...
                     (y - sen_set.y(arg.inner)).^2 + ...
                     (z - sen_set.z(arg.inner)).^2)'/3e8;
                          
        ki_sqrd = 1/deg_free*sum(((tt1 - tt2)/0.2e-6).^2);

function store_peaks_info(index,peaks,fID2)
        
    peaks = sortrows(peaks,2);
    
    fprintf(fID2,'\n%i',index);    

    j = 1;
    
    for i = 1:20
        try
            if i == peaks(j,2)
                fprintf(fID2,'\t%0.8f',peaks(j,1));
                j = j + 1;
            else
                fprintf(fID2,'\tNaN');
            end
        catch
            fprintf(fID2,'\tNaN');
        end
    end
    