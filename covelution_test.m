function covelution_test
% clc
% Let's load the nevest plotter values
try
    h=guidata(findall(0,'Tag','plotter2'));
    % Lets setup calibration values before plot
    h=calibration(h);
    g = h.g;
catch
    errordlg(sprintf('E-field plotting parameters are \ncoming from Plotter2. \n\nPlease run plotter2 first.\n'),'Plotter2 not found')
    return
end

% creating struct to store things
a=struct;


%% Creating the directory path for all files
% Opening sensor settings
settings=open('sensor_setting.mat');
% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

%Is manual Time shift for each sensor activated
if settings.man_tshiftOn==1
    tshift=settings.t_shift;
else
    tshift=zeros(1,60);
end

%Is offset voltages should be included?
if settings.vshiftOn==1
    vshift=settings.vshift;
else
    vshift=zeros(1,60);
end

%Is manual gain should be included?
if settings.gainOn==1
    gain=settings.gain;
else
    gain=zeros(1,60)+1;
end

% Is calibration should be included
if settings.plot_calibrated == 1
    factor = g.factor;
else
    factor = zeros(1,60)+1;
end

% Generating LP and CH Legend names
lp_legend=cell(1,60);
ch_legend=cell(1,60);

% Variable to store tshift info
tshift_str =[];

% Witch sensors has time shift turn on?
tmp = g.lpgraphs + g.chgraphs;

for i=1:20;
    lp_legend{i*3-2}=[settings.sen_IDs{i} ':lp1'];
    lp_legend{i*3-1}=[settings.sen_IDs{i} ':lp2'];
    lp_legend{(i*3)}=[settings.sen_IDs{i} ':lp3'];
    
    ch_legend{i*3-2}=[settings.sen_IDs{i} ':ch1'];
    ch_legend{i*3-1}=[settings.sen_IDs{i} ':ch2'];
    ch_legend{i*3}=[settings.sen_IDs{i} ':ch3'];
    
    if sum(tmp(i*3-2:i*3))> 0 && tshift(i*3)~= 0
        tmps = sprintf(' %s:%0.6fs', settings.sen_IDs{i},tshift(i*3));
        tshift_str = [ tshift_str tmps];
    end
end

% the variable for real legend
lg={};

%% Generating file names for ch (with header files) graphs and check those files are availble

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
        
        % Check whether the ch file and h file are exists
        a.absant_fn = {};
        if exist(filename,'file')==0
            a.absant_fn=[a.absant_fn filename];
        end
        
        if exist(hfilename,'file')==0
            a.absant_fn=[a.absant_fn hfilename];
        end
        
        % If file is not exist don't store the file name
        if exist(filename,'file')==0 || exist(hfilename,'file')==0
            a.ch_fn{i}='';
            a.h_fn{i}='';
        else
            a.ch_fn{i}=filename;
            a.h_fn{i}=hfilename;
        end
    else
        a.ch_fn{i}='';
        a.h_fn{i}='';
    end
end
wbh = waitbar(0,'test');
fg = figure;
set(fg,'visible','off')
hold all




%% Generating file name for LDAR data

if g.mm < 30
    ext=0;
else
    ext=30;
end

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);

if exist(ldar_fn,'file')~=0
    % Check whether the file is exists
    a.ldar_fn=ldar_fn;
else
    % If file is not exist don't store the file name
    a.ldar_fn='';
    % Absant file
    a.absant_fn=[a.absant_fn ldar_fn];
end

%% Time shift for LDAR and Linet
if settings.ldar_tshiftOn==1
    sn=settings.ldar_tshift_sn;
    x=settings.x;
    y=settings.y;
    z=settings.z;
    x0=x(sn)-settings.x0;
    y0=y(sn)-settings.y0;
    z0=z(sn)-settings.z0;
else
    x0=0;
    y0=0;
    z0=0;
end

mz = 3;



%% Loading and plotting ch data

t=[];   % Variable for time
y=[];   % Variable for voltage
empty_trigs=[]; % Variable for empty triggers

counter = 1;

% Array to compare
y_comp =[];
t_comp =[];

% Variables to store peaks
peak_t = [];
peak_v = [];
nOfPks = [];
sen_num = [];
delta_t = [];

for i=1:60
    
    try
        waitbar(0.3+i*0.005,wbh,'Loading Fast Antenna data','Name','Plotter Busy')
    catch
        return
    end
    
    t=[];   % Variable for time
    y=[];   % Variable for voltage
    % Plotting up to 60 fast antenna files
    if strcmp(a.ch_fn{i},'')==0
        % Load corresponding header file times
        fId = fopen(a.h_fn{i}, 'r');
        trigs  = fread(fId, inf, 'double') ;
        fclose( fId );
        % Introduce time shift before filtering time range
        trigs = trigs + tshift(i);
        % Finding the triggers between the given time range
        lol=length(trigs)- nnz(trigs>g.t1)+1;       % Index of the lower matrix element
        ul=nnz(trigs<g.t2);                      % Index of the upper matrix element
        % Let's find two more triggers from both ends
        if lol>1
            lol=lol-1;
        end
        
        if ul < length(trigs)
            ul=ul+1;
        end
        % Triggers between given time range and remove time shift
        trigs=trigs(lol:ul)-tshift(i);
        
        % Loading fast antenna data
        for j=1:length(trigs)
            sFa = epp_load_trigfile_time(a.ch_fn{i},trigs(j));
            %%%%%%%%%%%%%%%%% Reducing Sampling Rate %%%%%%%%%%%%%%%%%%%%%
            switch settings.chSumMethod
                case 1 % User needs downsampling
                    % find the current frequency
                    f=round(1/(sFa.t_i(2)-sFa.t_i(1)));
                    % Number of intervals
                    ni = round(f/settings.chfreq);
                    
                    if ni > 1
                        sFa.t_i = downsample(sFa.t_i,ni);
                        sFa.y_i = downsample(sFa.y_i,ni);
                    end
                    
                case 2 % User needs summing
                    % find the current frequency
                    f=round(1/(sFa.t_i(2)-sFa.t_i(1)));
                    % Number of intervals
                    ni = round(f/settings.chfreq);
                    
                    
                    if ni > 1
                        len = floor(length(sFa.t_i)/ni);
                        
                        temp1 = nan(ni,len);
                        temp2 = temp1;
                        
                        for indn = 1:ni
                            tempt= downsample(sFa.t_i,ni,indn-1);
                            tempy= downsample(sFa.y_i,ni,indn-1);
                            temp1(indn,:)=tempt(1:len);
                            temp2(indn,:)=tempy(1:len);
                            clear tempt tempy
    
                        end
                        
                        sFa.t_i=mean(temp1);
                        sFa.y_i=mean(temp2);
                    end
                    
                otherwise
                    % User needs nothing
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Keep a track of empty trigs
            if ~isempty(sFa.y_i)
                y=[y,NaN,sFa.y_i];
                t=[t,NaN,sFa.t_i];
            end
        end
        
        
            
        tshift(i);
        t=t+tshift(i);
        % Finding data in given time range
        lol=length(t)- nnz(t>g.t1)+1 ;      % Index of the lower matrix element
        ul=nnz(t<g.t2);                     % Index of the upper matrix element
        y=y(lol:ul);
        t=t(lol:ul);
        
        % Replace NaNs with zero so that filtfilt will work fine
        y(isnan(y))=0;
        
        if numel(t)>2
            
            % filter out low frequencies
            Fs = 1/(t(2)-t(1));
            [z,p] = butter(5,10000/(Fs/2),'high'); % Create a low-pass butterworth filter;
            % [z,p,k] = butter(n,Wn) designs an designs an order n lowpass digital Butterworth filter with normalized
            % cutoff frequency Wn. It returns the zeros and poles in length n column
            % vectors z and p, and the gain in the scalar k
            
            %y = filtfilt(z,p,y);    % filter the data.
            
                  
            if isempty(y_comp)
                y_comp = y;
                t_comp = t;
                dt = 0;

            else
                
                % finding convolution
                dt = t_comp(2) - t_comp(1);
                
                cy = xcorr(y,y_comp);
                
                [mm mmn ] = max(cy);
                
                lag = median([1:length(cy)]) - mmn;
                
                dt = lag*dt;
                
            end
            
            % Hilbert transform of voltages
%             y_hil = hilbert(y);
%             y_hil = sqrt(imag(y_hil).^2+real(y_hil).^2);
%             y_hil = y_hil / max(y_hil);
            
            % Normalize
            y = y/(max(y)-min(y));

            plot(t+dt,y)
            
            % Let's find peaks store those for future use
            [ pksLocs,  pks] = peakfinder(y,0.05,0.05,-1);
            
            sen_num = [sen_num i];
            nOfPks  = [nOfPks length(pksLocs)];
            peak_t  = [peak_t t(pksLocs) + dt NaN];
            peak_v  = [peak_v pks NaN];
            delta_t = [delta_t dt];
                        
            lg=[lg ch_legend{i}];
                       
        end
        
    end
end

box on
grid on
legend(lg)
delete(wbh)
%plot(peak_t,peak_v,'o')
set(fg,'visible','on')


%% find PBFA positions

% Let's find the longest record of peaks
[L maxInd] = max(nOfPks);
m = length(sen_num);

pk_t = NaN(L,m);
pk_v = pk_t;

start = 1;

% Make line array to squire matrix for  easy search
for j = 1:m
    pk_t(1:nOfPks(j),j) = peak_t(start:start+nOfPks(j)-1);
    pk_v(1:nOfPks(j),j) = peak_v(start:start+nOfPks(j)-1);
    start = start+nOfPks(j)+1;
end

% time window +-tw
tw = 2e-6;

arg.method = 2;
arg.outer = [];      % Outer sensors 
arg.t_out = [];

[AX,H1,H2]=plotyy(nan,nan,nan,nan);
hold(AX(2), 'on')



for i = 1:L
    
    temp = pk_t(i,maxInd);
    
    arg.inner = [];      % Inner sensors      
    arg.t_in  = [];
        
    for j=1:m
        [mini ind] = min(abs(pk_t(:,j)-temp));
        
        % is it within the time window
        if mini < tw
            arg.inner = [arg.inner ceil(sen_num(j)/3)];
            arg.t_in  = [arg.t_in  pk_t(ind,j) + delta_t(j)];
        end        
    end
   
    n = length(arg.t_in);
    
    if n > 5
        [x1,y1,z1,t1]=pbfa_finder(arg);
        fprintf('\n%.7f\t\t%0.1f\t\t%0.1f\t\t%0.1f\t\t%i',t1,x1,y1,z1,n)
        delt = sqrt((x1-x0)^2+(y1-y0)^2+(z1-z0)^2)/3e8;
        plot(AX(2),t1-delt,z1,'bo','MarkerFaceColor','b','MarkerSize',mz)
    end      
end

fprintf('\n')


[CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
    x0,y0,z0,0);
if ~isnan(DLS(:,10))
    plot (AX(2),DLS(:,10),DLS(:,8),'ko','MarkerFaceColor','k','MarkerSize',mz)
end

if ~isnan(CG(:,10))
    plot (AX(2),CG(:,10),CG(:,8),'go','MarkerFaceColor','g','MarkerSize',mz)
end

set(AX,'xlim',[g.t1 g.t2])


mmm=18000;
set(AX(2),'yLim',[0 mmm],'YTick',0:2000:mmm);
set(AX(2),'xtick',[])
set(AX, 'YColor', [0 0 0])


    


    
