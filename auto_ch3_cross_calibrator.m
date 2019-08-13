function auto_ch3_cross_calibrator
% This function was written to auto calibrate ch2 aginst ch2 for
% for outer sensors. For this perpose, this function will find CGLSS
% points that occured about the same distance from two sensors.
% The variable dist_tol control the percent tollerence of the distance.

%% Get plotter2 data
try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

settings = h.sen_set;
g = h.g;


%% User inputs
sn      = 11;                % Sensor number
t1      = 61200;            % Begining time
t2      = 86400;            % End time
dt      = 0.02e-3;          % +/- plot range
dist_tol = 5;               % Distance tolerance
channel  = 1;               % Calibrate channel 2 or 3

sid     = settings.sen_IDs{sn};
% directory name to save fig files
dn      = sprintf('C:/Users/Sumedhe/Desktop/Sumedhe/Calibration/Auto/');

%% Creating directories and file names to save

% dir names
d.d{1} = sprintf('%s/%s%s%s/%s/ch3_cross/with_K02/',dn,g.YYYY{:},g.MM{:},g.DD{:},sid);
d.d{2} = sprintf('%s/%s%s%s/%s/ch3_cross/with_K14/',dn,g.YYYY{:},g.MM{:},g.DD{:},sid);
d.d{3} = sprintf('%s/%s%s%s/%s/ch3_cross/with_K24/',dn,g.YYYY{:},g.MM{:},g.DD{:},sid);
d.d{6} = sprintf('%s/%s%s%s/%s/ch3_cross/with_K17/',dn,g.YYYY{:},g.MM{:},g.DD{:},sid);

% crating directories
if ~exist(d.d{1},'dir'); mkdir(d.d{1}); end
if ~exist(d.d{2},'dir'); mkdir(d.d{2}); end
if ~exist(d.d{3},'dir'); mkdir(d.d{3}); end
if ~exist(d.d{6},'dir'); mkdir(d.d{6}); end


% file name to save calibration data
d.fID(1) = fopen(sprintf('%s%s_ch3_cross_calibration_with_K02.txt',d.d{1},sid),'a+');
d.fID(2) = fopen(sprintf('%s%s_ch3_cross_calibration_with_K14.txt',d.d{2},sid),'a+');
d.fID(3) = fopen(sprintf('%s%s_ch3_cross_calibration_with_K14.txt',d.d{3},sid),'a+');
d.fID(6) = fopen(sprintf('%s%s_ch3_cross_calibration_with_K14.txt',d.d{6},sid),'a+');


%% Load CGLSS data
cglss_fn = sprintf('%s/cglss/%s/%s/KSCCGLSS%s%s%s.dat',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},g.MM{:},g.DD{:});

x0 = 0; y0 = 0; z0 = 0;

CG = CGLSS_extract(cglss_fn,t1,t2,500000,x0,y0,z0);
[L ~] = size(CG);

k02_cnt = 0;
k14_cnt = 0;
k17_cnt = 0;
k24_cnt = 0;

%% Load ch data
wbh = waitbar(0,' ','name',sprintf('%s Ch2 Cross Calibration',sid));

for i = 1:L
    
    msg = sprintf('Calibrating Ch2 channel... (%0.2f%%)',i/L*100);
    
    try
        waitbar(i/L,wbh,msg)
    catch
        disp('User stopped Ch1 auto calibration')
        fclose(d.fID(1));   fclose(d.fID(3));
        fclose(d.fID(2));   fclose(d.fID(6));
        return
    end
    
    tcg  = CG(i,1);
    xcg  = CG(i,2);
    ycg  = CG(i,3);
    
    % Distance to sensors
    D = sqrt((settings.x - xcg).^2 + (settings.y - ycg).^2);
   
    % Distance to inrersting sensor
    D0 = D(sn);
    
    % Is this point is same distance from any of the sensor
    dD = abs(D - D0);
    
    if dD(1)/D0*100 < dist_tol
        calibrate_with(1,tcg,dt,settings,g,sid,sn,D,d,channel)
        k02_cnt = k02_cnt + 1;
    end
    
    if dD(2)/D0*100 < dist_tol
        calibrate_with(2,tcg,dt,settings,g,sid,sn,D,d,channel)
        k14_cnt = k14_cnt + 1;
    end
    
    
    if dD(3)/D0*100 < dist_tol
        calibrate_with(3,tcg,dt,settings,g,sid,sn,D,d,channel)
        k24_cnt = k24_cnt + 1;
    end
    
    if dD(6)/D0*100 < dist_tol
        calibrate_with(6,tcg,dt,settings,g,sid,sn,D,d,channel)
        k17_cnt = k17_cnt + 1;
    end
end

try delete(wbh); end;
try fclose(d.fID(1)); end;
try fclose(d.fID(2)); end;
try fclose(d.fID(3)); end;
try fclose(d.fID(6)); end;





function calibrate_with(csn,tcg,dt,settings,g,sid,sn,D,d,channel)

data = load_data(tcg,dt,settings,g,sn,csn,D,channel);

if ~isempty(data.t1) && ~isempty(data.t2)
    
    r1 = range(data.v1);
    r2 = range(data.v2);
    
    if r1 > 0.005 && r2 > 0.005
                
        % Normalize data to 100km
        data.v1 = data.v1*D(csn)/100000;
        data.v2 = data.v2*D(sn)/100000;
        
        % Time shift CSN so that both peaks are the same
        [mm ind1] = min(data.v1);
        [mm ind2] = min(data.v2);
        
        mini_tshift = data.t1(ind1) - data.t2(ind2);
        
        data.t1 = data.t1-mini_tshift;        
        
        r1 = range(data.v1);
        r2 = range(data.v2);
        gain = r2/r1;
        
        
        m1 = max(data.v1);
        m2 = max(data.v2/gain);
        
        t0 = min(data.t1);
        
        fg = figure;
        hold all
        plot((data.t1-t0)*1e6,data.v1)
        plot((data.t2-t0)*1e6,data.v2/gain+m1-m2)
        tit = sprintf('%s-%s-%s Auto Calibration %s ch%1.1i\nGain = %0.3f    ch2_range = %0.4f    ch1_range = %0.4f',...
            g.YYYY{:},g.MM{:},g.DD{:},sid,channel,gain,r1,r2);
        title(tit)
        legend({sprintf('%s:ch%1.1i',settings.sen_IDs{csn},channel) ...
                sprintf('%s:ch%1.1i',settings.sen_IDs{sn},channel)})
        box on
        %xlim([0 2*dt]*1e6)
        ylabel('dV (Volts)')
        xlabel(['\mus after ' sprintf('%0.6f',t0) 's UT'])
        
        secnd = floor(tcg);
        milsec = round((tcg-secnd)*1000);
        figfn = sprintf('%s_ch%1.1i_cross_%s%s%s_%5.5i_%3.3i.fig',...
            sid,channel,g.YYYY{:},g.MM{:},g.DD{:},secnd,milsec);
        
        saveas(fg,[d.d{csn} figfn])
        delete(fg)
         
        %Write data in to the file
        fprintf(d.fID(csn),'%s\t%0.3f\t%0.4f\t%0.4f\n',...
            figfn,gain,r1,r2);
    end
    
end


function data = load_data(tcg,dt,settings,g,sn,csn,D,channel)
    % Let's make summing for ch2 data
    settings.chSumMethod = 1;
   
    sid1 = settings.sen_IDs{csn};
    sid2 = settings.sen_IDs{sn};
    
    % If sensor is BCC then it should be 1 MHz
    if sn == 5
        settings.chfreq(sn*3-(3-channel)) = 1e6;
        settings.chfreq(csn*3-(3-channel)) = 1e6;
    else
        % Make ch frequency to 5MHz
        settings.chfreq(sn*3 - (3-channel)) = 5e6;
        settings.chfreq(csn*3 -(3-channel)) = 5e6;
    end
   
    % Generate ch file names according to the current CGLSS poing
    hhmmss = sec2hhmmss(tcg);
    g.hh = str2double(hhmmss(1:2));
    g.mm = floor(str2double(hhmmss(4:5))/5)*5;

    % Base folder
    bfolder=sprintf('%s/%s/%s/%s/', ...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

    % finding ch2 file name
    ch2fn1 = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
        bfolder,sid1,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,channel);

    ch2hf1 = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
        bfolder,sid1,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,channel);
    
    % finding ch1 file name
    ch2fn2 = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
        bfolder,sid2,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,channel);

    ch2hf2 = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
        bfolder,sid2,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,channel);
    
    t1 = tcg + D(csn)/3e8 - dt;
    t2 = t1 + 2*dt;   
    tshift=settings.t_shift(csn*3);
    
    if exist(ch2fn1,'file') && exist(ch2hf1,'file')
        [data.t1 data.v1 ~] = FA_Extract1(ch2fn1,ch2hf1,t1,t2,tshift,settings,csn*3-(3-channel));
    else
        fprintf('Missing file: %s\n',ch2fn1)
        data.t1 = [];
        data.v1 = [];
    end
    
    t1 = tcg + D(sn)/3e8 - dt;
    t2 = t1 + 2*dt;
    tshift=settings.t_shift(sn*3);
    
    if exist(ch2fn2,'file') && exist(ch2hf2,'file')
        [data.t2 data.v2 ~] = FA_Extract1(ch2fn2,ch2hf2,t1,t2,tshift,settings,sn*3-(3-channel));
    else
        fprintf('Missing file: %s\n',ch2fn2)
        data.t2 = [];
        data.v2 = [];
    end

  
        
    
