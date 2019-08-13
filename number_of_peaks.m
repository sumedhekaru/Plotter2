function  number_of_peaks
% calculate number of begining pulses in time window by counting number of
% negative peaks in the range. Specify the channel and the time range from
% the plotter2 main window.

%% Get plotter2 data
try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end
settings = h.sen_set;
g = h.g;

%% User inputs
% senosr number
sn = 1;
ch = 3;
f = 1e6;                        % Frequency
peak_type = -1;                  % +1 for IC and -1 for CG
thr = 0.7;                      % threshold

%%
sid = settings.sen_IDs(sn);
% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

% finding ch2 file name
fn = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
    bfolder,sid{:},g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ch);

hf = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
        bfolder,sid{:},g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ch);

% Make ch frequency to 1MHz
settings.chfreq(sn*3-2) = f;
settings.chfreq(sn*3-1) = f;
settings.chfreq(sn*3)   = f;

tshift = settings.t_shift((sn-1)*3+ch);

[t v ~] = FA_Extract1(fn,hf,g.t1,g.t2,tshift,settings,(sn-1)*3+ch);

% Is calibration should be included
if settings.plot_calibrated
    factor = g.factor;
else
    factor = zeros(1,60)+1;
end
factor = factor((sn-1)*3+ch);

v = v*factor;

% High pass filter
% Replace NaNs with zero so that filtfilt will work fine
v(isnan(v))=0;

% filter out low frequencies
Fs = 1/(t(2)-t(1));
[z,p] = butter(5,5000/(Fs/2),'high'); % Create a low-pass butterworth filter;
v = filtfilt(z,p,v);    % filter the data.

pL = peakfinder(v,thr,thr,peak_type);

%% Getting LDAR points

if g.mm < 30
    ext=0;
else
    ext=30;
end

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);

if ~exist(ldar_fn,'file')
    % If file is not exist don't store the file name
   ldar_fn='';
end

% Loading LDAR2 data
[CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        0,0,0,0);
    
[L1 L2] = size(DLS);
%     assignin('base', 'lx', [DLS(:,6);CG(:,6)])
%     assignin('base', 'ly', [DLS(:,7);CG(:,7)])

R = mean((DLS(:,6).^2+DLS(:,7).^2).^0.5)/1000;





tit = sprintf('%s-%s-%s      %0.7f     %s:Ch%i \n%i pulses within %.1fms     %i (%0.2f) LDAR2 detection    R = %0.1f km', ...
    g.YYYY{:},g.MM{:},g.DD{:},g.t1,settings.sen_IDs{sn},ch, length(pL),(g.t2-g.t1)*1000,L1,L1/length(pL)*100, R);
figure

plot(t,v)
hold all
tools2fig
xlim([g.t1 g.t2])
plot(t(pL),v(pL),'mo','markerfacecolor','m','markersize',2)
xlabel('Time (s)')
ylabel('V/m')
title(tit)

