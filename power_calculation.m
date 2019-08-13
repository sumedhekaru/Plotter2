function power_calculation(sns,ch,g,sen_set,mean_power_curve,x0,y0,z0,t0)

%% User inputs

% sensor number(s)
if nargin < 1
    sns = [1 2 3 6 7 8 9 10 ];
end

% Channel number 
if nargin < 2
    ch = 1;
end

% get plotter2 data
if nargin < 3
    h=guidata(findall(0,'Tag','plotter2'));
    g = h.g;
    sen_set = h.sen_set;
end

% Turn on mean power curve as default
if nargin < 5
    mean_power_curve = 1;
end
%Get the location(s) of this pulse
if nargin < 6
    [x0 y0 z0] = get_location(g,sen_set);
end

g.hh = floor(t0/3600);
g.mm = floor(((t0 - 3600*g.hh)/60)/5)*5;
    
    
%% Start the program
R = sqrt((x0 - sen_set.x).^2+(y0 - sen_set.y).^2+(z0 - sen_set.z).^2);

fg = figure;
set(fg,'Visible','off')
hold all
lg = {};

% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

data = [];
sensors = [];
Ptots = [];

wbh = waitbar(0,'Power Spectral Density','name','PSD');

Lsns = length(sns);

R = sqrt((sen_set.x-x0).^2 + (sen_set.y - y0).^2 + (sen_set.z - z0).^2);

for sn = sns
    %% Create file name
    % sensor ID
    
    waitbar(sn/Lsns,wbh)
    
    sid = sen_set.sen_IDs{sn};
   
    
    % create file name
    fn = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
        bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ch);
    
    hfn = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
        bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ch);
    
    % Chanel ID (Ex K14 channel 1 is 4)
    if ch == 3;
        chID = sn*3;
    elseif ch == 2;
        chID = sn*3 - 1;
    else
        chID = sn*3 - 2;
    end
    
       
    %% Get the data
    if exist(fn,'file')
        
        g.t1 = t0 + R(sn)/3e8 - 500e-6;
        g.t2 = t0 + R(sn)/3e8 + 500e-6;
        
        [t, y ,ch_freq_str] = FA_Extract1(fn,hfn,g.t1,g.t2,sen_set.t_shift(chID),sen_set,chID);
        
        %calibrated data
        vshift=sen_set.vshift; % Manual offset
        factor = g.factor;      % Calibration
        
        y = y*factor(chID)+vshift(chID);
        
        %filter data
        [t,y] = ch_high_pass(t,y,2000);
        
        % remove DC
        %offset = nanmean(y(1:50));
        %y = y-offset;
        
        if range(y) < 1 % May be no useful data here. Let's ignore it.
            t = nan(size(t));
        end
    else
        t = nan;
        y = nan;
    end
    
    if length(t) ~= sum(isnan(t))       
        
        
        %% Power calculations
        % figure
        % plot(t,y)
        % tools2fig
        
        
        %% FFT
        % Sample time
        for i = 1:1000
            T = t(i+1)-t(i);
            if ~isnan(T)
                break
            end
        end
        
        Fs = 1/T ;                   % Sampling frequency
        L = length(t);                     % Length of signal
        NFFT = 2^nextpow2(L); % Next power of 2 from length of y
        %NFFT = 4096
        
        %% FFT
        % Y = fft(y,NFFT)/L;
        % f = Fs/2*linspace(0,1,NFFT/2+1);
        % % Plot single-sided amplitude spectrum.
        % %plot(log10(f(1:end)),2*abs(Y(1:NFFT/2+1)),'rp-')
        % title('Power Density')
        % xlabel('Log(f) Hz')
        % ylabel('|Y(f)|')
        %
        % Pyy = 2*Y(1:NFFT/2+1).*conj(Y(1:NFFT/2+1));
        %
        % figure
        % box on
        % hold all
        % plot(log10(f),Pyy,'r.')
        % x = log10(f(2:end));
        % % new x
        % xn = 0:0.01:x(end);
        % yn = csaps(x,Pyy(2:end));
        % sn = fnxtr(yn);
        %
        % % Le'ts poly fit
        %
        % %plot(xn,yn,'r')
        % %fnplt(sn,[0 x(end)])
        % vn = fnval(sn,xn);
        % ind = find(vn <0);
        % vn(ind) = [];
        % xn(ind) = [];
        % plot(xn,vn,'r')
        
        
        %ylim([0 5])
        
        
        
        
        %% Pwelch analysis
        
        noverlap = 0;
        %figure

        %[Pxx,fx] = pwelch(y,[],noverlap,NFFT,Fs);
        %[Pxx,fx] = pwelch(y,[],[],NFFT,Fs);
        [Pxx, fx] = pwelch(y,[], NFFT/2, NFFT, Fs);
        
        
        
        % Here Pxx = E^2. Total source power can be obtained using
        %       Pt = 4*pi/3*R^2*epsilon0*c*Ep^2 See Krider 1992 JGR
        
        epsilon0 = 8.854e-12;
        c        = 299792458/1.0003;
        
        Pt = 4*pi/3*R(sn)^2*epsilon0*c*Pxx;
        
        
        plot(log10(fx),Pt);
        Ptot = trapz(fx,Pt);
        lg = [lg sprintf('%i-%s%s-%0.2e W %0.1f km',sn,sen_set.sen_IDs{sn},ch_freq_str,Ptot,R(sn)/1000)];
        
        data = [data Pt];
        
        % available sensors
        sensors = [sensors, sn];
        
        % total powers
        Ptots = [Ptots, Ptot];
        
    end
    
end

delete(wbh);

meanP = mean(data,2);
plot(log10(fx'),meanP,'linewidth',2);
Pm = trapz(fx,meanP);


box on
ylabel('Power/freqency (W/Hz)')
xlabel('log_{10}[Frequency (Hz)]')
title(sprintf('Power spectrum %s-%s-%s     %0.6f - %0.6f\n<P_{tot}> = %0.2e W',...
        g.YYYY{:},g.MM{:},g.DD{:},g.t1,g.t2,Pm))
legend([lg,'Average'],'Location','NorthWest')
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8])
tools2fig

% Save data in figure handle for later use
figD.dataP = data;
figD.dataMeanP = meanP;
figD.sns = sensors;
figD.x0 = x0;
figD.y0 = y0;
figD.z0 = z0;
figD.totalPs = Ptots;
figD.totalMeanP = Pm;
figD.fx = fx;
guidata(gcf,figD)
set(fg,'Visible','on')

if mean_power_curve
    figure
    plot(log10(fx'),meanP,'linewidth',2);
    box on
    ylabel('Power/freqency (W/Hz)')
    xlabel('log_{10}[Frequency (Hz)]')
    title(sprintf('Power spectrum %s-%s-%s     %0.6f - %0.6f\n<P_{tot}> = %0.2e W',...
        g.YYYY{:},g.MM{:},g.DD{:},g.t1,g.t2,Pm))
    legend('Mean Power','Location','NorthWest')
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8])
    tools2fig
end




% % Create spectrum object and call msspectrum method.
% h = spectrum.welch;
% hmss = msspectrum(h,y,'Fs',Fs); % retuns a data object
% 
% psdd = psd(h,y,'Fs',Fs);
% 
% Ptot = trapz(psdd.Frequencies,psdd.Data)
% 
% assignin('base', 'hmss', psdd)
% 
% %sqrt(max(hmss.data))
% 
% figure
% % View the power in dB = 10*log10(hmss.data).
% %plot(hmss)
% %plot(psdd.Frequencies,psdd.Data,'r.--')
%  plot(psdd)

function [x, y, z, type] = get_location(g,sen_set)


x = [];
y = [];
z = [];
type = [];

if g.mm < 30
    ext=0;
else
    ext=30;
end

% Is there a PBFA point?
fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
data = pbfaExtract(fn,g.t1,g.t2,str2double(sen_set.ldar_r),0,0,0,0);

if ~isnan(data(1))
    x = nanmean(data(3));
    y = nanmean(data(4));
    z = nanmean(data(5));
    type = 1;
    if z > 3000
     return
    end
end

% Is there a PBFA-Auto point?
fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
data = pbfaExtract(fn,g.t1,g.t2,str2double(sen_set.ldar_r),0,0,0,0);

if ~isnan(data(1))
    x = nanmean(data(3));
    y = nanmean(data(4));
    z = nanmean(data(5));
    type = 2;
    return
end
    

