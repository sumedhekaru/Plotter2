function plot_hsv_intensity
% This function will plot HSV intensity data with slow/fast antenna data.
% HSV intensity data should availbel under HSVI folder. Don't forget to
% enter apropriate values for user inputs.

%% Get plotter2 data
try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

settings = h.sen_set;
g = h.g;

%% User inputs
% First choose date time and time range from plotter2

% HSVI file name (txt file for HSV Intensity, could be generated using "extract_data_from_figure")
hsvifn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110724\FI_Jul24-Flash01part_Intensity_figure_cropped_180-0_80_120.txt';
%alignT = 59734.3122298;
%alignfN = -44674;
alignT = 64131.754013;
alignfN = -46470;
sn = 6;  % To plot data
ch = 3;   % channel to plot
f  = 5e6; % Data frequency
normFrame = -46440; % When this is zero, it will normalize the entire range. If not sepecify the exact image number to normalize.
minFrame = -46478; % When this is zero, program auto calculate the minimum frame. Otherwise this will consider as the minimum brightness frame.

%% Get the HSV intensity file
fID = fopen(hsvifn,'r');
data = textscan(fID,'%f %f');
N = data{1};
I = data{2};
I0 = I(1);


%% Plot data


% Load FFI channel chanle "ch" data data

sid = settings.sen_IDs{sn};
tshift=settings.t_shift(sn*3);
ind = (sn-1)*3+ch;
settings.chfreq(ind) = f;

if settings.plot_calibrated == 1
    factor = g.factor(ind);
    ylab = 'dE (V/m)';
else
    factor = 1;
    ylab = 'dV (Volts)';
end

bfolder=sprintf('%s/%s/%s/%s/', ...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

fn = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
    bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ch);

hf = sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
    bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ch);

% Load data
if exist(fn,'file') && exist(hf,'file')
    [t v freq_str] = FA_Extract1(fn,hf,g.t1,g.t2,tshift,settings,ind);
    v = v*factor;
else
    fprintf('Missing file: %s\n',fn)
    t = [];
    v = [];
end

if isempty(t)
    return
end

I_orig = I;
N_orig = N;


T = alignT + (N-alignfN)*20e-6;
% Trim data
lol = sum(T < g.t1)+1;
ul  = sum(T < g.t2);
T = T(lol:ul);
I = I(lol:ul);
N = N(lol:ul);

if minFrame == 0
    [minIntensity, mInd] = min(I);
    fprintf('Minimum intensity frame number %i.\n',N(mInd))
else
    ind = find(N_orig==minFrame);
    minIntensity = I_orig(ind);
end

 
if normFrame == 0   
    [maxIntensity, mInd] = max(I);  
    fprintf('Intensty normalized to frame number %i.\n',N(mInd))
else
    ind = find(N_orig==normFrame);
    maxIntensity = I_orig(ind);
end

I = I - minIntensity;
I = I/(maxIntensity - minIntensity); 
  
  
figure
tools2fig
%[AX,H1,H2] = plotyy(t,v,T,I);
[AX,H1,H2] = plotyy(t,v,T-20e-6,I,'stairs');


set(H2,'Color','r','LineStyle','-','marker','d','markerfacecolor','r','markersize',3); 

%mmm = max(I)*1.2;

ylimits = get(AX(1),'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(AX(1),'YTick',ylimits(1):yinc:ylimits(2))

%set(AX(2),'yLim',[0 mmm],'YTick',-mmm:mmm/5:mmm);

ylimits = get(AX(2),'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(AX(2),'YTick',ylimits(1):yinc:ylimits(2))


set(AX(2),'xtick',[])
set(AX, 'YColor', [0 0 0])
set(get(AX(2),'Ylabel'),'String','Normalized Light Intensity');
set(AX,'XAxisLocation','bottom')
set(AX(2),'YAxisLocation','right')


legend({[sid freq_str] ,'HSV-Intensity'})
linkaxes(AX,'x')
xlim([g.t1 g.t2])

xlabel('Time (s)')
ylabel(ylab)

box on
grid on
[~,tit] = fileparts(hsvifn);
title(tit)




%% Do some work to plot with Hilbert Transform

%plot(t,v)
%hold all
    
%% Apply a hiss pass filter

% Replace NaNs with zero so that filtfilt will work fine
v(isnan(v))=0;

% filter out low frequencies
Fs = 1/(t(2)-t(1));

if isnan(Fs)
    Fs = 1/(t(5) - t(4));
end


[z,p] = butter(5,5000/(Fs/2),'high'); % Create a High-pass butterworth filter;

try
    v = filtfilt(z,p,v);    % filter the data.
catch
    disp('Couldn''t filter data')
    return
end

v = abs(hilbert(v));

figure
[AX,H1,H2] = plotyy(t,v,T,I);


ylimits = get(AX(1),'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(AX(1),'YTick',ylimits(1):yinc:ylimits(2))


ylimits = get(AX(2),'YLim');
yinc = (ylimits(2)-ylimits(1))/5;
set(AX(2),'YTick',ylimits(1):yinc:ylimits(2))


set(AX(2),'xtick',[])
set(AX, 'YColor', [0 0 0])
set(get(AX(2),'Ylabel'),'String','Normalized Light Intensity');



legend({[sid freq_str] ,'HSV-Intensity'})
linkaxes(AX,'x')
xlim([g.t1 g.t2])

xlabel('Time (s)')
ylabel('magnitude of HT')

box on
grid on
[~,tit] = fileparts(hsvifn);
title(tit)
tools2fig





    






