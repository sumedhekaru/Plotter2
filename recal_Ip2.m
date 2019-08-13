function recal_Ip2
clc
% Peak current calculations of PBFA using Ch3 data has messed up. It is
% always showing a lower value than CGLSS. As a solution Ch1 data has used
% to calculate peak currents. This program will find the new peak currents
% for already calculated PBFA points for comparison.

% The differences between recal_Ip vs recal_Ip2 are
%   - Calculate peakcurrents for both PBFA and PBFA-Auto
%   - Instead of only doing for RSs, it will do for all PBFA points


%% User inputs

%% Start the program
% Get Plotter2 data
try h=guidata(findall(0,'Tag','plotter2'));
catch
    disp('Run Plotter2 first')
    return
end

g = h.g;
settings = h.sen_set;

x0 = 0;
y0 = 0;
z0 = 0;


% Load PBFA data
if g.mm < 30
    ext=0;
else
    ext=30;
end

pbfa_fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    settings.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

if ~exist(pbfa_fn,'file')
    % Check whether the file is exists
    fprintf('File not exist\n%s\n',pbfa_fn); return;
end

PBFA=pbfaExtract(pbfa_fn,g.t1,g.t2,str2double(settings.ldar_r),...
    x0,y0,z0,0);

if ~isnan(PBFA(1,1))    
    L = length(PBFA(:,1));
    
    % Turn on only ch1 data
    g.chgraphs = [1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 ...
        1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 ...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    
    [ch_fn h_fn] = generate_ch_fn(settings,g);
    tshift=settings.t_shift;
    factor = g.factor;
    
    
    fprintf('\nPBFA-auto data\nTime\t\t\tOld Ip\tNew_Ip\n')
    
    for k = 1:L
        t0 = PBFA(k,2);
        
        % Let's see there is a CGLSS point to go with this
        t1 = t0 - 0.05e-3;
        t2 = t0 + 0.05e-3;
        
        x = PBFA(k,3);
        y = PBFA(k,4);
        z = PBFA(k,5);
        
        Ip = PBFA(k,7);
        
        R = sqrt((settings.x-x).^2 + (settings.y-y).^2+ (settings.z-z).^2);
        
        Eps = nan(1,10);
        
        for i=1:3:30
            
            if strcmp(ch_fn{i},'')==0
                
                % distance
                sn = ceil(i/3);
                r = R(sn);
                
                % Travel time to sensor
                travelT = r/3e8;
                
                
                
                [t y ch_freq_str] = FA_Extract1(ch_fn{i},h_fn{i},t1+travelT,t2+travelT,tshift(i),settings,i);
                
                
                
                if ~isempty(t)
                    
                    [yF yH] = hill_tra(t,y,2000);
                    yH = abs(yH);
                    
                    yPeakH = max(yH);
                    yPeakF = min(yF);
                    
                    if abs(yPeakF*factor(i)) > 1
                        Eps(sn) = yPeakH;
                    end
                    
                end
            end
        end
        % Calculate the peak current
        
        Ipn = cal_Ip_using_ch1(Eps,R);
        
        fprintf('%12.6f\t%5.1f\t%5.1f\n',t0,Ip,Ipn)
        
        
    end
end

pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    settings.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

if ~exist(pbfa_fn,'file')
    % Check whether the file is exists
    fprintf('File not exist\n%s\n',pbfa_fn); return;
end

PBFA=pbfaExtract(pbfa_fn,g.t1,g.t2,str2double(settings.ldar_r),...
    x0,y0,z0,0);

if ~isnan(PBFA(1,1))
    
    L = length(PBFA(:,1));
    
    % Turn on only ch1 data
    g.chgraphs = [1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 ...
        1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 ...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    
    [ch_fn h_fn] = generate_ch_fn(settings,g);
    tshift=settings.t_shift;
    factor = g.factor;
    
    
    fprintf('\nPBFA data\nTime\t\t\tOld Ip\tNew_Ip\n')
    
    for k = 1:L
        t0 = PBFA(k,2);
        
        % Let's see there is a CGLSS point to go with this
        t1 = t0 - 0.05e-3;
        t2 = t0 + 0.05e-3;
        
        x = PBFA(k,3);
        y = PBFA(k,4);
        z = PBFA(k,5);
        
        Ip = PBFA(k,7);
        
        R = sqrt((settings.x-x).^2 + (settings.y-y).^2+ (settings.z-z).^2);
        
        Eps = nan(1,10);
        
        for i=1:3:30
            
            if strcmp(ch_fn{i},'')==0
                
                % distance
                sn = ceil(i/3);
                r = R(sn);
                
                % Travel time to sensor
                travelT = r/3e8;
                
                
                
                [t y ch_freq_str] = FA_Extract1(ch_fn{i},h_fn{i},t1+travelT,t2+travelT,tshift(i),settings,i);
                
                
                
                if ~isempty(t)
                    
                    [yF yH] = hill_tra(t,y,2000);
                    yH = abs(yH);
                    
                    yPeakH = max(yH);
                    yPeakF = min(yF);
                    
                    if abs(yPeakF*factor(i)) > 1
                        Eps(sn) = yPeakH;
                    end
                    
                end
            end
        end
        % Calculate the peak current
        
        Ipn = cal_Ip_using_ch1(Eps,R);
        
        fprintf('%12.6f\t%5.1f\t%5.1f\n',t0,Ip,Ipn)
        
        
    end
end






function Ip = cal_Ip_using_ch1(Eps,R)

ind = find(R(1:10) >= 30000);

% factor for y = mx calibration
% xFac = [264.0 217.1 276.7 NaN NaN 245.2 345.2 105.8 84.9 199.9];
% Ip = Eps(ind).*((R(ind)/100e3).^1.13).*xFac(ind)

% % factor and C for y = mx + c calibration
% xFac = [274.6 225.2 289.2 NaN NaN 255.2 362.2 112.0 88.8  212.7];
% c    = [ -1.4  -1.3  -1.6 NaN NaN  -1.5  -1.8  -1.9 -1.5   -2.1];
% Ip = Eps(ind).*((R(ind)/100e3).^1.13).*xFac(ind)+c(ind);

% y = mx+c calibration using collocated PBFA ch1
xFac = [265.0 219.7 277.0 NaN NaN 247.5 336.5 105.1 84.3  205.1];
c    = [ -0.5  -0.5  -.7  NaN NaN  -0.5  -0.3  -0.8 -0.5   -0.9];
Ip = Eps(ind).*((R(ind)/100e3).^1.13).*xFac(ind)+c(ind);


Ip = nanmean(Ip);



Ip = nanmean(Ip);




function [ch_fn h_fn] = generate_ch_fn(settings,g)
% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

% checking witch graphs are on
ch_on=g.chgraphs;

for i=1:30
    if ch_on(i)
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
        if ~exist(filename,'file')
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

