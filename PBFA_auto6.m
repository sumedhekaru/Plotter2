function PBFA_auto6(h)
% in this function aim is to find location of positive NBP

if nargin < 1
    try; h=guidata(findall(0,'Tag','plotter2'));
    catch; disp('Run plotter2 first!'); return;
    end
    
    %t1 = 77636.333; % start time
    %t2 = 77637.339; % end
    t1 = h.g.t1;
    t2 = h.g.t2;
else
    t1 = h.g.t1;
    t2 = h.g.t2;
end

settings = h.sen_set;
g = h.g;
tshift=settings.t_shift;





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

    %fn = 'C:\Users\sumedhe\Desktop\pbfa_test.txt';
    
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
d.fID = fID;
index = 1;

fprintf(fID,'\n');

fID2 = fopen([fn(1:end-4) '_peaks.txt'],'a+');

g.chgraphs = [1 0 0 1 0 0 1 0 0 0 0 0 1 0 0  ...
    1 0 0 1 0 0 1 0 0 1 0 0 1 0 0];

% g.chgraphs = [0 0 1 0 0 1 0 0 1 0 0 0 0 0 1 0 0  ...
%     1 0 0 1 0 0 1 0 0 1 0 0 1 0 0];


%figure(100)
%tools2fig
%hold all

% Trial time shifts 
t_shift_tr = zeros(1,60);

d.sen_set = settings;
d.maxKiSqrd = 3;

warning('off','signal:findpeaks:largeMinPeakHeight')

d.indx = 0;


[ch_fn, h_fn] = generate_ch_fn(settings,g);

wbh1 = waitbar(0,'Calculating PBFA');

% Initialize dd structure
d.dd(1,30).t = [];
d.dd(1,30).y = [];
d.sn_ind = [];
d.sns = [];

for i= 1:30
    
    if strcmp(ch_fn{i},'')==0
        % Get data
        [d.dd(i).t, d.dd(i).y] = FA_Extract1(ch_fn{i},h_fn{i},g.t1,g.t2,tshift(i),settings,i);
        % record some statistics
        d.dd(i).max = nanmax(d.dd(i).y);
        d.dd(i).min = nanmin(d.dd(i).y);
        d.dd(i).mean = nanmean(d.dd(i).y);
        d.dd(i).std = nanstd(d.dd(i).y);
        d.dd(i).range = d.dd(i).max - d.dd(i).min;
        d.sn_ind = [d.sn_ind i];
        d.sns = [d.sns ceil(i/3)];
        
        % Get the hilbert transform and filtered data
        [d.dd(i).yF, d.dd(i).yH] = hill_tra(d.dd(i).t,d.dd(i).y,1000);
        
        % Normalize data
        d.dd(i).yN = d.dd(i).yF/range(d.dd(i).yF);
        
    end
    
    try
        waitbar(i/30,wbh1,sprintf('Loading 5 MHz data (%.0f%%)',i/30*100));
    catch
        return
    end
    
end

d = get_general_cross_corr(d);


% Time window times
tWinL = g.t1;
tWinR = g.t2;

% Let's go through time windows


% Find lol and ul for a single data strip
%figure; hold all
for j = 1:30
    if ~isempty(d.dd(j).t)
        lol = 1;
        ul = length(d.dd(j).t);
        %plot(d.dd(j).t(lol:ul) - d.dd(j).timeDiff1,d.dd(j).yN(lol:ul))
        d.dd(j).lol = lol;
        d.dd(j).ul = ul;
    end
end

% Upsample data to 10 MHz
d = up_sample_data(d);

% Let's find new cross correlation for this time window
%d = get_window_correlation(d);

% Let's find peaks
d = find_peaks_for_strip(d);

% find PBFA points
d = find_PBFAs(d);






delete(wbh1)
fclose(fID);
fclose(fID2);

function d = find_PBFAs(d)

% sensor i with max number of peaks
i0 = d.maxNofPeaksi;

L = d.maxNofPeaks;

% dt thershold
dt_tr = 1e-6;

arg.outer = [];      % Outer sensors
arg.t_out = [];
arg.sen_set = d.sen_set;
arg.method = 5;
%[ax p1 p2] = plotyy(nan,nan,nan,nan);
%hold(ax(2), 'on')

%figure(100); hold all
% cluster peaks and find pbfa
for i = 1:L
    %t0  = d.dd(i0).pkts(i) - d.dd(i0).timeDiff2 - d.dd(i0).timeDiff1;
    t0  = d.dd(i0).pkts(i) - d.dd(i0).timeDiff1;
    arg.inner = [];
    arg.t_in = [];
    
    for j = 1:30
        if ~isempty(d.dd(j).t)
            %td2 = d.dd(j).timeDiff2;
            td2 = 0;
            td1 = d.dd(j).timeDiff1;
            
            ts = d.dd(j).pkts - td1 - td2;
            dts = abs(ts-t0);
            
            % minimum dt
            [tmin, mind] = min(dts);
            
            if tmin < dt_tr
                % This could be a similar pulse
                arg.inner = [arg.inner, ceil(j/3)];
                arg.t_in = [arg.t_in, ts(mind) + td1 + td2];
            end
        end
    end
    
    if length(arg.inner) >= 5
        [xs,ys,zs,t1]=pbfa_finder(arg);
        ki_sqrd = cal_ki_sqrd(t1,xs,ys,zs,arg,d.sen_set); 
        
        if ki_sqrd < d.maxKiSqrd
            d.indx = d.indx + 1;
            fprintf('%4.4i\t%12.7f\t%10.1f\t%10.1f%10.1f%10.1f\n',d.indx,t1,xs/1000,ys/1000,zs/1000,ki_sqrd)
            %plot(t1,zs/100,'ro')
            
            snsStr = sprintf('%i,',sort(arg.inner));
            pulseID = -inf;
            Ip = 0;
            fprintf('%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%4.4i\t%0.1f\t%s\tNaN\tNaN\tNaN\tNaN\t%0.1f\n',...
                d.indx,t1,xs,ys,zs,pulseID,Ip,snsStr(1:end-1),ki_sqrd)
            
            fprintf(d.fID,'%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\tNaN\tNaN\tNaN\tNaN\t%0.1f\n',...
                d.indx,t1,xs,ys,zs,pulseID,Ip,snsStr(1:end-1),ki_sqrd);
            %plot3(xs,ys,zs,'ro','markerfacecolor','r','markersize',2)
        end
    end
end

%daspect([1 1 1])
%view(3)
%linkaxes(ax,'x')
%xlim([min(d.dd(1).t), max(d.dd(1).t)])

function d = up_sample_data(d)
% We will upsample the filtered data

%figure; hold all;
for i = 1:30
    lol = d.dd(i).lol;
    ul = d.dd(i).ul;
    if ~isempty(d.dd(i).t(lol:ul))
        y = d.dd(i).yF(lol:ul);
        t = d.dd(i).t(lol:ul);
        
        % Upsample data to 10 Mhz
        y = resample(y,10,1);
        t = t(1):1e-7:t(end);
        
               
        L1 = length(t);
        L2 = length(y);
        L = min([L1 L2]);
        
        d.dd(i).yU = y(1:L);
        d.dd(i).tU = t(1:L);
        
        %plot(t(1:L) - d.dd(i).timeDiff1,y(1:L))
    end
end




function d = find_peaks_for_strip(d)

% Threshold
thr = 0.0001;

% Minimum Peak Height
mph = 0.004;

% Minimum distance to close by peak (number of samples)
mpd = 20;


%figure; hold all;

L = -inf;

for i = 1:30
    lol = d.dd(i).lol;
    ul = d.dd(i).ul;
    if ~isempty(d.dd(i).t(lol:ul))
        y = d.dd(i).yU;
        t = d.dd(i).tU;

        [pks, locs] = findpeaks(y,'MinPeakHeight',mph,'MinPeakDistance',mpd);
        %[pks, locs] = max(y);
        
        d.dd(i).pkts = t(locs);
        
        % Number of peaks
        L2 = length(locs);
        
        if L < L2
            L = L2;
            d.maxNofPeaks = L2;
            d.maxNofPeaksi = i;
        end
        
        
        
        %td2 = d.dd(i).timeDiff2;
        td2 =0 ;
       
        td1 = d.dd(i).timeDiff1;
        
        %plot(t-td1-td2,y)        
        %plot(t(locs)-td1-td2,-pks,'o')
        
        
    end
end




function d = get_window_correlation(d)
% The good sensor to correlate (found from generel cross correlation
% function)
 ref_i = d.ref_i;


% Let's find cross correlations between reference sensor and others




%figure; hold all;
for i = 1:30
    if ~isempty(d.dd(i).t)

        
        [acor,lag] = xcorr(d.dd(i).yU,d.dd(ref_i).yU);
        [~,I] = max(abs(acor));
        %plot(acor)
        lagDiff = lag(I);
        td2 = lagDiff/10e6;
        d.dd(i).timeDiff2 = td2;
        
        %td1 = d.dd(i).timeDiff1;
        %plot(d.dd(i).tU-td1-td2,d.dd(i).yU)
        
    end
end


function d = get_general_cross_corr(d)


% Let's find none saturated sensor with good signal to correlate others
% I think one with the medium data range should be not too close and not too
% far


ranges = [d.dd.range];
ranges = sort(ranges);

medianV = ranges(round(length(ranges)/2));

for i = 1:30
    if d.dd(i).range == medianV
        % We found a sensor to work with
        ref_i = i;
        break
    end
end

d.ref_i = ref_i;
d.ref_sns = ceil(ref_i/3);


% Let's find cross correlations between reference sensor and others

%figure; hold all;
for i = 1:30
    if ~isempty(d.dd(i).t)
        
        
        [acor,lag] = xcorr(d.dd(i).yN,d.dd(ref_i).yN,250);
        [~,I] = max(abs(acor));
        %plot(acor)
        lagDiff = lag(I);
        timeDiff = lagDiff/1e6;
        d.dd(i).timeDiff1 = timeDiff;
        %plot(d.dd(i).t-timeDiff,d.dd(i).yN)
        
    end
end




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
for i=[3 6 9 15 18 21 24 27 30]-2
    
    
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
                    
                    % Increase sampling rate
                    
                    tTemp = t(1):1e-7:t(end);
                    yTemp = spline(t,yn,t(1):1e-7:t(end));
                    
                    % CG type peaks
                    %[ pksLocs,  pks] = peakfinderSum(yTemp,settings.selection,settings.pk_th,-1);
                    [pks,pksLocs] = findpeaks(yTemp,'MinPeakHeight',0.001);
                    
                    sns     = [sns zeros(1,length(pksLocs))+sn];
                    nOfPks  = [nOfPks length(pksLocs)];
                    peak_t  = [peak_t tTemp(pksLocs) + dt];
                    peak_v  = [peak_v pks NaN];
                    delta_t = [delta_t zeros(1,length(pksLocs))+dt];
                    
                    % if ch3 saturated, no point of storing hilbert peak
                    % info to calculate Ip
                    if saturated
                        peak_vH = [peak_vH NaN(size(pksLocs))];
                    else
                        peak_vH = [peak_vH yTemp(pksLocs)];
                    end
                    
                    %plot(t(pksLocs)+dt,yn(pksLocs),'o')
                    
                    % IC type peaks
                    [ pksLocs,  pks] = peakfinderSum(yTemp,settings.selection,settings.pk_th,+1);
                    
                    sns1     = [sns1 zeros(1,length(pksLocs))+sn];
                    nOfPks1  = [nOfPks1 length(pksLocs)];
                    peak_t1  = [peak_t1 tTemp(pksLocs) + dt];
                    peak_v1  = [peak_v1 pks NaN];
                    delta_t1 = [delta_t1 zeros(1,length(pksLocs))+dt];
                    
                    
                    % if ch3 saturated, no point of storing hilbert peak
                    % info to calculate Ip
                    if saturated
                        peak_vH1 = [peak_vH1 NaN(size(pksLocs))];
                    else
                        peak_vH1 = [peak_vH1 yTemp(pksLocs)];
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