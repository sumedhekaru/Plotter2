function auto_ch3_calibrator
% This function was written to auto calibrate ch3 aginst ch1

%% Get plotter2 data
try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

settings = h.sen_set;
g = h.g;


%% User inputs

sn      = 9;                % Sensor number
t1      = 71907;            % Begining time
t2      = 86400;            % End time
dt      = 0.005e-3;             % +/- plot range
maxG    = 12;               % Maximum gain
minG    = 7;                % Minimum gain
sid     = settings.sen_IDs{sn};

%%
% directory name to save fig files
dn      = 'C:\Users\daqop\Desktop\Sumedhe\Calibration\Auto\';

dn      = sprintf('%s/%s%s%s/%s/ch3_auto/',dn,g.YYYY{:},g.MM{:},g.DD{:},sid);

% Create the directory if not exist
if ~exist(dn,'dir')
    mkdir(dn)
end

% file name to save calibration data
fn  = sprintf('%s_ch3_auto_calibration_data.txt',sid);
fID = fopen([dn fn],'a+');

%% Load CGLSS data
cglss_fn = sprintf('%s/cglss/%s/%s/KSCCGLSS%s%s%s.dat',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},g.MM{:},g.DD{:});

x0 = 0; y0 = 0; z0 = 0;

CG = CGLSS_extract(cglss_fn,t1,t2,500000,x0,y0,z0);
L    = length(CG);

%% Turn on the correct graphs
g.lpgraphs = zeros(1,60);
g.chgraphs = zeros(1,60);
g.chgraphs((sn-1)*3+2) = 1; % ch2
g.chgraphs((sn-1)*3+3) = 1; % ch3

tH = -1;

wbh = waitbar(0,' ','name','Ch3 Auto Calibration');

for i = 1:L
    
    msg = sprintf('Calibrating Ch3 channel... (%0.2f%%)',i/L*100);
    
    try
        waitbar(i/L,wbh,msg)
    catch
        disp('User stopped Ch3 auto calibration')
        fclose(fID);
        return
    end
    
    tcg  = CG(i,1);
    xcg  = CG(i,2);
    ycg  = CG(i,3);
    
    if tcg > tH
        tH = ceil(tcg);
        data = load_data(tcg,settings,g,sid,sn);
    end
    
    % If the data loaded lets find
    if ~isempty(data.ch2t) || ~isempty(data.ch3t)
        % time shift due to propagation
        t_shift = sqrt((xcg - settings.x(sn)).^2 + (ycg - settings.y(sn)).^2)/3e8;
        
        % Time range to plot
        tr1 = tcg + t_shift - dt;
        tr2 = tcg + t_shift + dt;
        
        % Clip ch2 data
        lol = sum(data.ch2t < tr1)+1;
        ul  = sum(data.ch2t < tr2);
        ch2t = data.ch2t(lol:ul);
        ch2v = data.ch2v(lol:ul);
        
        % Clip ch3 data
        lol = sum(data.ch3t < tr1)+1;
        ul  = sum(data.ch3t < tr2);
        ch3t = data.ch3t(lol:ul);
        ch3v = data.ch3v(lol:ul);
        
        if ~isempty(ch2t) && ~isempty(ch3t)
            
            r1 = range(ch2v);
            r2 = range(ch3v);
            gain = r2/r1;
            
            
            if r1 > 0.005 && r2 > 0.005  && gain < maxG && gain > minG
                
                m1 = max(ch2v);
                m2 = max(ch3v/gain);
                
                t0 = min(ch2t);
                
                fg = figure;
                hold all
                plot((ch2t-t0)*1e6,ch2v)
                plot((ch3t-t0)*1e6,ch3v/gain+m1-m2)
                tit = sprintf('%s-%s-%s Auto Calibration %s ch3\nGain = %0.3f    ch2_range = %0.4f    ch3_range = %0.4f',...
                    g.YYYY{:},g.MM{:},g.DD{:},sid,gain,r1,r2);
                title(tit)
                legend({'ch2' 'ch3'})
                box on
                xlim([0 2*dt]*1e6)
                ylabel('dV (Volts)')
                xlabel(['\mus after ' sprintf('%0.6f',t0) 's UT'])
                
                secnd = floor(tcg);
                milsec = round((tcg-secnd)*1000);
                figfn = sprintf('%s_ch3_auto_%s%s%s_%5.5i_%3.3i.fig',...
                    sid,g.YYYY{:},g.MM{:},g.DD{:},secnd,milsec);
                
                saveas(fg,[dn figfn])
                delete(fg)
                
                % Write data in to the file
                fprintf(fID,'%s\t%0.3f\t%0.4f\t%0.4f\n',...
                    figfn,gain,r1,r2);
            end
        end
    end
    
end

try delete(wbh); end;
try fclose(fID); end;



function data = load_data(tcg,settings,g,sid,sn)

    t1 = floor(tcg);
    t2 = t1 + 1;
    
    tshift=settings.t_shift(sn*3);
    
    % Generate ch file names according to the current CGLSS poing
    hhmmss = sec2hhmmss(tcg);
    g.hh = str2double(hhmmss(1:2));
    g.mm = floor(str2double(hhmmss(4:5))/5)*5;

    % Base folder
    bfolder=sprintf('%s/%s/%s/%s/', ...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

    % finding ch2 file name
    ch2fn = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch2', ...
        bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss);

    ch2hf = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h2', ...
        bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss);
    
    % finding ch3 file name
    ch3fn = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch3', ...
        bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss);

    ch3hf = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h3', ...
        bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss);
    
    if exist(ch2fn,'file') && exist(ch2hf,'file')
        [data.ch2t data.ch2v ~] = FA_Extract1(ch2fn,ch2hf,t1,t2,tshift,settings,sn*3-1);
    else
        fprintf('Missing file: %s\n',ch2fn)
        data.ch2t = [];
        data.ch2v = [];
    end
    
    if exist(ch3fn,'file') && exist(ch3hf,'file')
        [data.ch3t data.ch3v ~] = FA_Extract1(ch3fn,ch3hf,t1,t2,tshift,settings,sn*3);
    else
        fprintf('Missing file: %s\n',ch3fn)
        data.ch3t = [];
        data.ch3v = [];
    end

  
        
    
