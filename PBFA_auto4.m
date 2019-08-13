function PBFA_auto4(h)
% 2013-09-25 Changed from PBFA_auto3 to get better PBFA results. Here, the aim is not
% find every single pulse, insetead find few big pulses.

if nargin < 1
    try; h=guidata(findall(0,'Tag','plotter2'));
    catch; disp('Run plotter2 first!'); return;
    end
    
    t1 = 77636.333; % start time
    t2 = 77637.339; % end
else
    t1 = h.g.t1;
    t2 = h.g.t2;
end

settings = h.sen_set;
g = h.g;
tshift=settings.t_shift;



if (t2 - t1) < 2
    TW = t2-t1;               % data load window
else
    TW = 2;
end

tw = 2e-3;              % total time window to load

if tw > TW
    tw = TW;
end

wbh1 = waitbar(0,sprintf('Loading data and calculating PBFA (%.2f%%)',0));

%fn = 'C:\Users\sumedhe\Desktop\PBFA_auto_TOAM5\PBFA_auto_test_20110814_new.txt';

% Generate file name based on g.t1

if ispc

    if g.mm < 30
        ext=0;
    else
        ext=30;
    end

    fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
    disp('Data will be saved in the following file')
    disp(['    ' fn])
    
else

    PBFA_files = ls([tempdir 'PBFA_*']); 
    [s1 s2] = size(PBFA_files);
    fn = sprintf('%sPBFA_%5.5i.txt',tempdir,s1/2+1);
    
    disp('**********************************************************')
    disp('MAC computers cannot write into Windows hard drives remotely.')
    disp('Therefore info will write to the following file')
    disp(['     ' fn])
    disp('Please contact Sumedhe with the above file to add data to PBFA directory')
    disp('**********************************************************\n')
           
end

% Write the header if not exit
if exist(fn,'file')
    fID = fopen(fn,'a+');
else
    
    % Make the folder it not available
    [pathstr,~,~] = fileparts(fn);
    
    if ~exist(pathstr,'dir')
        mkdir(pathstr)
    end
    
    fID = fopen(fn,'a+');
    fprintf(fID,'Index\ttpbfa\txpbfa\typbfa\tIpbfa\tNpbfa\n');
end

index = 1;

fprintf(fID,'\n');

fID2 = fopen([fn(1:end-4) '_peaks.txt'],'a+');

n = floor(TW/tw);

tLs = t1:TW:t2;
tRs = tLs + TW;
tRs(end) = t2;

n2 = length(tLs);

g.t1 = t1-tw;
g.t2 = t1;

g.chgraphs = [0 0 1 0 1 1 0 1 1 0 0 0 0 1 1 ...
              0 0 1 0 1 1 0 1 1 0 1 1 0 1 1];
          
settings.pk_th = 0.01; % Peak Threshold
settings.selection = 200; % Peak windows in number of samples

%figure(100)
%tools2fig
%hold all


for k=1:n2
    
    [ch_fn h_fn] = generate_ch_fn(settings,g);
    
   
    tL = tLs(k);
    tR = tRs(k);
    
    if tL >= tR
        break
    end
    
    for i=[2 3 5 6 8 9 17 18 20 21 23 24 26 27 29 30] % 3:3:30
        
        if strcmp(ch_fn{i},'')==0            
            [dd.t{i} dd.y{i} ~] = FA_Extract1(ch_fn{i},h_fn{i},tL,tR,tshift(i),settings,i);
        else
            dd.t{i}=[];
            dd.y{i}=[];            
        end
        
        try
            waitbar(i/30,wbh1,sprintf('Loading 5 MHz data (%.0f%%)',i/30*100));
        catch
            return
        end
        
    end
    
    % Number of elements in a window (assuming 5 Mhz)
    nelm = round(tw/0.2e-6);
    g.lol = 1-nelm;
    g.ul = 0; 
       
    for m=1:n
        g.t1 = g.t1+tw;
        g.t2 = g.t2+tw;
        g.lol = g.lol + nelm;
        g.ul = g.ul + nelm;
        %fprintf('%.4f\n',g.t;1)     
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


% wbh = waitbar(0,sprintf('Loading Fast Antenna data. (%.2f%%)',0));


% figure
% hold all
% tools2fig

sns = [];
nOfPks  = [];
peak_t  = [];
peak_v  = [];
peak_vH = [];
delta_t = [];

sns1 = [];
nOfPks1  = [];
peak_t1  = [];
peak_v1  = [];
peak_vH1 = [];
delta_t1 = [];


%figure(100)
%hold all
for i=[3 6 9 18 21 24 27 30]
    
    L = length(dd.t{i});
    saturated = 0;
        
    if ~isempty(dd.t{i})
        
        lol = g.lol;
        
        if g.ul > L
            ul = L;
        else
            ul = g.ul;
        end

        t = dd.t{i}(lol:ul);
        y = dd.y{i}(lol:ul);
        
        % find out ch3 is saturated or not        
        if range(y) < 0.01 
           saturated = 1;
           try
            y = dd.y{i-1}(lol:ul);            
           end
        end
        
        % sensor
        sn = ceil(i/3);

        if ~isempty(t)
  
            [yn yH] = hill_tra(t,y,10000);
            %plot(t,yn)

            if ~isempty(yn)
                
                yH = abs(yH);
                yn1 = yH;
                tn1 = t;
                
                %save 1 data set for comparison
                if ~exist('y_comp','var') && exist('yn1','var')
                    
                    y_comp = yn1;
                    t_comp = tn1;
                    dt = 0;
                    
                    d.t{1} = tn1;
                    d.y{1} = yn1;
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
                    
                    lag = median(1:length(cy)) - mmn;
                    
                    dt1 = lag*dt;
                    
                    [mxv mxi1] = max(y_comp);
                    [mxv mxi2] = max(yn1);
                    
                    dt2 = t_comp(mxi1) - tn1(mxi2);
                    
                    if dt2 > dt1 - 0.001 && dt2 < dt1 + 0.001
                        dt = dt1;
                        % disp('Hooray!')
                    else
                        dt = dt2;
                        % disp('nope!')
                    end
                    
                    data_cnt = data_cnt + 1;
                    d.t{data_cnt} = tn1;
                    d.y{data_cnt}  = yn1;
                    d.sn{data_cnt} = sn;
                    
                    
                else
                    
                    dt = 1; %
                    
                end
                
                
                %plot(t+dt,yH)
                %plot(t+dt,yn)
                
                
                
                % Find peaks
                if dt < 240e-6 % time to travel 70km
                    
                    % CG type peaks
                    [ pksLocs,  pks] = peakfinderSum(yn,settings.selection,settings.pk_th,-1);                    
                                      
                    sns     = [sns zeros(1,length(pksLocs))+sn];
                    nOfPks  = [nOfPks length(pksLocs)];
                    peak_t  = [peak_t t(pksLocs) + dt];
                    peak_v  = [peak_v pks NaN];                    
                    delta_t = [delta_t zeros(1,length(pksLocs))+dt];
                    
                    % if ch3 saturated, no point of storing hilbert peak
                    % info to calculate Ip
                    if saturated
                        peak_vH = [peak_vH NaN(size(pksLocs))];
                    else
                        peak_vH = [peak_vH yH(pksLocs)];
                    end
                    
                   %plot(t(pksLocs)+dt,yn(pksLocs),'o')
                    
                    % IC type peaks
                    [ pksLocs,  pks] = peakfinderSum(yn,settings.selection,settings.pk_th,+1);
                    
                    sns1     = [sns1 zeros(1,length(pksLocs))+sn];
                    nOfPks1  = [nOfPks1 length(pksLocs)];
                    peak_t1  = [peak_t1 t(pksLocs) + dt];
                    peak_v1  = [peak_v1 pks NaN];
                    delta_t1 = [delta_t1 zeros(1,length(pksLocs))+dt];
                    
                    
                    % if ch3 saturated, no point of storing hilbert peak
                    % info to calculate Ip
                    if saturated
                        peak_vH1 = [peak_vH1 NaN(size(pksLocs))];
                    else
                        peak_vH1 = [peak_vH1 yH(pksLocs)];
                    end
                    
                    
                    %plot(t(pksLocs)+dt,pks,'o')
                end
            end
        end
    end
end

% Find CG types locations
if length(nOfPks) > 4
    pulse_type = -1;
    index = find_pbfa_positions(sns,peak_t,delta_t,settings,index,d,fID,fID2,peak_vH,pulse_type);
end

% Find IC types locations
if length(nOfPks1) > 4
    pulse_type = 1;
    index = find_pbfa_positions(sns1,peak_t1,delta_t1,settings,index,d,fID,fID2,peak_vH1,pulse_type);
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

for i=1:30
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

    
function index=find_pbfa_positions(sen_num,peak_t,delta_t,settings,index,d,fID,fID2,peak_vH,pulse_type)
%% find PBFA positions

% time window +-tw
tw = 3.5e-6;


% [AX,H1,H2]=plotyy(nan,nan,nan,nan);
% hold(AX(2), 'on')


% figure
% plot(peak_t,peak_t,'ro')

peaks = NaN(length(peak_t),3);

peaks(:,1) = peak_t';
peaks(:,2) = sen_num';
peaks(:,3) = delta_t';
peaks(:,4) = peak_vH';

peaks = sortrows(peaks,1);

peak_t = peaks(:,1);
sen_num = peaks(:,2);
delta_t = peaks(:,3);
peaks_vH = peaks(:,4);

try    
    time = peak_t(1);
    ts = peak_t(1) - delta_t(1);
    sns  = sen_num(1);
    vHs  = peaks_vH(1);
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
            vHs  = [vHs peaks_vH(i)];
        else
            n = length(ts);
            
            if n > 4
                [ts sns vHs] = clean_t(time,ts,sns,vHs);
            end
            
            n = length(ts);
            
            if n > 4               
                                        
                index = calculate_PBFA(ts,sns,index,settings,d,fID,fID2,vHs,pulse_type);
                
                %return
                
                time = peak_t(i);
                sns  = sen_num(i);
                ts = peak_t(i) - delta_t(i);
                vHs = peaks_vH(i);
            else
                time = peak_t(i);
                sns  = sen_num(i);
                ts = peak_t(i) - delta_t(i);
                vHs = peaks_vH(i);
            end
        end
    end    
end


function index = calculate_PBFA(ts,sns,index,settings,d,fID,fID2,vHs,pulse_type)
 

arg.outer = [];      % Outer sensors 
arg.t_out = [];
arg.inner = sns;
arg.t_in = ts;
arg.sen_set = settings;

n = length(ts);


% Method 5 answer
arg.method = 5;
[xs,ys,zs,t1]=pbfa_finder(arg);
ki_sqrd = cal_ki_sqrd(t1,xs,ys,zs,arg,settings);

% Method 3 answer
arg.method = 3;
[xs2,ys2,zs2,t2]=pbfa_finder(arg);
ki_sqrd2 = cal_ki_sqrd(t2,xs2,ys2,zs2,arg,settings);

if ki_sqrd < ki_sqrd2
    arg.method = 5;
else
    arg.method = 3;
    xs = xs2; ys = ys2; zs = zs2; t1 = t2;
    ki_sqrd = ki_sqrd2;
end

% Determine we are looking for a RS?
if ~isreal(zs) || zs < 3000
    arg.method = 6;
    [xs,ys,zs,t1]=pbfa_finder(arg);
    pulse_type = 0;
end



if n >= 5 && ki_sqrd < 5
    
    Ip = Ip_cal2(xs,ys,zs,settings,sns,vHs);
    
    snsStr = sprintf('%i,',sort(sns));
    
    if pulse_type < 0
        pulseID = -1000-n-10*arg.method;
    else
        pulseID = pulse_type*1000+n+10*arg.method;
    end
    
    fprintf('%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%4.4i\t%0.1f\t%s\tNaN\tNaN\tNaN\tNaN\t%0.1f\n',...
        index,t1,xs,ys,zs,pulseID,Ip,snsStr(1:end-1),ki_sqrd)
    
    fprintf(fID,'%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\tNaN\tNaN\tNaN\tNaN\t%0.1f\n',...
        index,t1,xs,ys,zs,pulseID,Ip,snsStr(1:end-1),ki_sqrd);
    
    %plot(t1,zs/5000,'k*')
    
    store_peaks_info(index,[ts',sns'],fID2)
    
    index = index+1;
    
end





function [ts sns vHs] = clean_t(time,ts,sns,vHs)
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
        vHs(rep)    =[];
    end
end

function Ip = Ip_cal2(xs,ys,zs,settings,sns,vHs)


m =  [472.4 387.1 418.7 NaN NaN 411.9 676.4  204.2 165.6 209.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
r = sqrt((settings.x - xs).^2 + (settings.y - ys).^2 + (settings.z - zs).^2)/1000;
rs = r(sns);
ms = m(sns);
Ips = ms.*vHs.*((rs/100).^1.13);
indx = find(rs > 30);

Ip = nanmean(Ips(indx));


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
    
function    [ pksLocs,  pks] = peakfinderSum(yn,selection,pk_th,pk_type);



% Number of elements
L = length(yn);

% Number of slices
N = L/selection;

% Derivative of voltage
dydtn = nan(1,L);
dydtn(2:L) = yn(2:end) - yn(1:end-1);

% 2nd derivativve of the voltage
dydtn2 = nan(1,L);
dydtn2(2:L) = dydtn(2:end) - dydtn(1:end-1);

lol = 1:selection:L;
ul  = selection-1:selection:L;
ul(end) = L;

%figure
%hold all
%tools2fig

pksLocs = NaN(1,200);
pks = pksLocs;
k = 0;

for i = 1:N
    y = yn(lol(i):ul(i));
    %t = tn(lol(i):ul(i));
    
    if pk_type < 0
        [mm ind] = min(y);
        indx = (i-1)*selection + ind;
        meanv = sum(y)/selection;
        
        if (mm-meanv) < -pk_th  && indx > 1 && indx < L
            %&& dydtn2(indx-1) < 0 && dydtn2(indx+1) > 0
            
            k = k+1;
            pksLocs(k) = indx;
            pks(k) = mm;
            %plot(tn(indx),yn(indx),'ro')
        end
    else
        [mm ind] = max(y);
        indx = (i-1)*selection + ind;
        meanv = sum(y)/selection;
       
        if (mm-meanv) > pk_th  && indx > 1 && indx < L
            %&& dydtn2(indx-1) < 0 && dydtn2(indx+1) > 0
            
            k = k+1;
            pksLocs(k) = indx;
            pks(k) = mm;
            %plot(tn(indx),yn(indx),'ro')
        end
    end
end

% Cleanup
pksLocs = pksLocs(1:k);
pks = pks(1:k);

% figure(200)
% hold all
% plot(tn,yn)
% plot(tn(pksLocs),pks,'o')