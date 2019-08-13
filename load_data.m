function a = load_data(g,sen_set,a)
% This function will try to load all data that user selected.

%% Creating the directory path for all files
% Opening sensor settings

% creating struct to store things
if nargin < 2
    a=struct;
end

wbh = waitbar(0,'Loading data','name','Please wait');

settings=sen_set;

% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

% Location data folder for local access
loc_dir = settings.loc_dir;

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

a.factor = factor;
a.vshift = vshift;
a.gain = gain;


% Generating LP and CH Legend names
lp_legend=cell(1,60);
ch_legend=cell(1,60);

% Variable to store tshift info
tshift_str =[];

% Witch sensors has time shift turn on?
tmp = g.lpgraphs + g.chgraphs;


% Do we need to time shift graphs to points?
if ~settings.plot_tshiftOn
    settings.plot_tshiftType = 0;
end

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



%% Generating file names for lp graphs and check those files are availble

% checking witch graphs are on
lp_on=g.lpgraphs;

% Save absant file extentions
a.absant_fn={};

for i=1:60
    if lp_on(i)==1
        % finding the file extention number
        ext=mod(i,3);
        if ext==0
            ext=3;
        end
        
        % Finding the stattion ID
        sid=settings.sen_IDs{ceil(i/3)};
        
        
        % finding the file name
        filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.lp%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        % Check whether the file is exists
        if exist(filename,'file')~=0
            % Check whether the file is exists
            a.lp_fn{i}=filename;
        else
            % If file is not exist don't store the file name
            a.lp_fn{i}='';
            % Absant file
            a.absant_fn=[a.absant_fn filename];
        end
    else
        a.lp_fn{i}='';
    end
end


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

if g.mm < 30
    ext=0;
else
    ext=30;
end

%% Generating file name for Field mill data
for i=1:34
    if g.KSC_FM(i)
        
        bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});
        
        % finding the file name
        filename=sprintf('%s/FM/%s/%s/%s/KSC%2.2i%s%s%s_%2.2i%2.2i.txt', ...
            settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},i,g.YYYY{:},g.MM{:},g.DD{:},g.hh,ext);
        % Check whether the file is exists
        if exist(filename,'file')~=0
            % Check whether the file is exists
            a.fm_fn{i}=filename;
        else
            % If file is not exist don't store the file name
            a.fm_fn{i}='';
            % Absant file
            a.absant_fn=[a.absant_fn filename];
        end
    else
        a.fm_fn{i}='';
    end
end


%% Generating file name for LDAR data
if g.ldar==1 || g.cglss == 1 || (settings.plot_tshiftOn && settings.plot_tshiftType == 1)
    
    dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
        -datenum(str2double(g.YYYY),0,0));
    
    ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);
    
    if exist(ldar_fn,'file')~=0
        % Check whether the file is exists
        a.ldar_fn=ldar_fn;
    else
        % If file is not exist don't store the file name
        a.ldar_fn='';
        % Absant file
        if g.ldar==1 || g.cglss == 1
            a.absant_fn=[a.absant_fn ldar_fn];
        end
    end
else
    a.ldar_fn='';
end

%% Generating file name for LINET

if g.linet==1 || (settings.plot_tshiftOn && settings.plot_tshiftType == 3)
    
    linet_fn=sprintf('%s/LINET/%s/%s/%s/linet_%s%s%s_%2.2i%2.2i.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
    if exist(linet_fn,'file')~=0
        % Check whether the file is exists
        a.linet_fn=linet_fn;
    else
        % If file is not exist don't store the file name
        a.linet_fn='';
        % Absant file
        if g.linet
            a.absant_fn=[a.absant_fn linet_fn];
        end
    end
else
    a.linet_fn='';
end
 
%% Generating file name for NLDN 
% Temprary this would be PBFA2 (Auto PBFA)

if g.nldn==1 || (settings.plot_tshiftOn && settings.plot_tshiftType == 4)
    
    nldn_fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
    if exist(nldn_fn,'file')~=0
        % Check whether the file is exists
        a.nldn_fn=nldn_fn;
    else
        % If file is not exist don't store the file name
        a.nldn_fn='';
        % Absant file
        if g.nldn
            a.absant_fn=[a.absant_fn nldn_fn];
        end
    end
else
    a.nldn_fn='';
end


%% Generating file name for PBFA
% Did not completed yet!!!

if g.pbfa ==1 || (settings.plot_tshiftOn && settings.plot_tshiftType == 2)
   
    pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
    if exist(pbfa_fn,'file')~=0
        % Check whether the file is exists
        a.pbfa_fn=pbfa_fn;
    else
        % If file is not exist don't store the file name
        a.pbfa_fn='';
        % Absant file
        if g.pbfa 
            a.absant_fn=[a.absant_fn pbfa_fn];
        end
   end
else
    a.pbfa_fn='';
end

%% Generating file name for OLD PBFA

if g.pbfa_old ==1 || (settings.plot_tshiftOn && settings.plot_tshiftType == 5)
   
    pbfa_old_fn=sprintf('%s/PBFA_old/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
    if exist(pbfa_old_fn,'file')~=0
        % Check whether the file is exists
        a.pbfa_old_fn=pbfa_old_fn;
    else
        % If file is not exist don't store the file name
        a.pbfa_old_fn='';
        % Absant file
        if g.pbfa 
            a.absant_fn=[a.absant_fn pbfa_old_fn];
        end
   end
else
    a.pbfa_old_fn='';
end


%% Generating File name for dEdt
% This is just the trigger file name
if g.dedt==1
    
    dedt_fn=sprintf('%s/dedt/%s/%s/%s/%2.2i%2.2i/BCC_%s%s%s_%2.2i%2.2i.hh3',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.hh,ext,g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
     
   if exist(dedt_fn,'file')~=0
        % Check whether the file is exists
        a.dedt_fn=dedt_fn;
    else
        % If file is not exist don't store the file name
        a.dedt_fn='';
        % Absant file
        a.absant_fn=[a.absant_fn dedt_fn];
    end
else
    a.dedt_fn='';
end

%% Generating File name for dEdt_trig
% This is just the trigger file name
if g.dedt_trig==1
    
    dedt_trig_fn=sprintf('%s/dedt/%s/%s/%s/%2.2i%2.2i/BCC_%s%s%s_%2.2i%2.2i.hh2',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.hh,ext,g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
     
   
    if exist(dedt_trig_fn,'file')~=0
        % Check whether the file is exists
        a.dedt_trig_fn=dedt_trig_fn;
    else
        % If file is not exist don't store the file name
        a.dedt_trig_fn='';
        % Absant file
        a.absant_fn=[a.absant_fn dedt_trig_fn];
    end
else
    a.dedt_trig_fn='';
end


%% Generating file name for NLDN2 

if g.nldn2==1 || (settings.plot_tshiftOn && settings.plot_tshiftType == 6)
    
    nldn2_fn=sprintf('%s/NLDN2/%s/%s/NLDN2_%s%s%s.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.YYYY{:},...
        g.MM{:},g.DD{:});
    
    if exist(nldn2_fn,'file')~=0
        % Check whether the file is exists
        a.nldn2_fn=nldn2_fn;
    else
        % If file is not exist don't store the file name
        a.nldn2_fn='';
        % Absant file
        if g.nldn2
            a.absant_fn=[a.absant_fn nldn2_fn];
        end
    end
else
    a.nldn2_fn='';
end

%% Letting user know about missing file

a.absant_fn=sort(a.absant_fn);

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


%% Loading LDAR2 data
ld_lin='';

if strcmp(a.ldar_fn,'')==0
    ld_lin=[ld_lin '-LDAR-'];
    try
        waitbar(0.1,wbh,'Loading LDAR data','Name','Plotter Busy')
    catch
        return
    end
    % Load ldar data
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0,0);
   
    %For future use to extract ldar data in to the workspace
%     assignin('base', 'lt', x1)
%     assignin('base', 'lx', [DLS(:,6);CG(:,6)])
%     assignin('base', 'ly', [DLS(:,7);CG(:,7)])
%     assignin('base', 'lz', y1)
else
    CG=NaN(1,10);    
    DLS=NaN(1,10);       
end

a.CG = CG;
a.DLS = DLS;

if strcmp(a.linet_fn,'')==0
    ld_lin=[ld_lin '-LINET-'];
    try
        waitbar(0.15,wbh,'Loading LINET data','Name','Plotter Busy')
    catch
        return
    end
    LINET=linetExtract2(linet_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0);
    % LINET Column DISCRIPTION
    % 1 - time
    % 2 - Distance from the time correcting sensor
    % 3 - time shifted time
    % 6 - x distance
    % 7 - y distance
    % 8 - z distance
    %x1=[x1; LINET(:,3)];
    %y1=[y1; LINET(:,8)];
else
    LINET=NaN(1,8);
end

a.LINET = LINET;

%%%%%%%%This suppose to be nldn but PBFA2 will plot temporrary
if strcmp(a.nldn_fn,'')==0
    ld_lin=[ld_lin '-PBFA2-'];
    try
        waitbar(0.2,wbh,'Loading LINET2 data','Name','Plotter Busy')
    catch
        return
    end
    PBFA2=pbfaExtract(a.nldn_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0);
else
    PBFA2=NaN(1,6);    
end

a.PBFA2 = PBFA2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(a.pbfa_fn,'')==0
    ld_lin=[ld_lin '-PBFA-'];
    try
        waitbar(0.25,wbh,'Loading PBFA data','Name','Plotter Busy')
    catch
        return
    end
    PBFA=pbfaExtract(pbfa_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0);
    
    % PBFA column Discription
    %   1 - distance to the pulse from the given sensor
    %   2 - Occuring time
    %   3 - x
    %   4 - y
    %   5 - z
    %   6 - Detection time at the sensor
    
    
    %x1=[x1; PBFA(:,6)];
    %y1=[y1; PBFA(:,5)];
else
    PBFA=NaN(1,6);
end

a.PBFA = PBFA;

if strcmp(a.pbfa_old_fn,'')==0
    ld_lin=[ld_lin '-PBFAO-'];
    try
        waitbar(0.245,wbh,'Loading PBFA-O data','Name','Plotter Busy')
    catch
        return
    end
    PBFA_old=pbfaExtract(a.pbfa_old_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0);
else
    PBFA_old=NaN(1,6);    
end

a.PBFA_old = PBFA_old;

%%%% Loading NLDN2 data 
if strcmp(a.nldn2_fn,'')==0
    ld_lin=[ld_lin '-NLDN2-'];
    try
        waitbar(0.25,wbh,'Loading NLDN2 data','Name','Plotter Busy')
    catch
        return
    end
    [NLDN2c NLDN2g] = nldnExtract(a.nldn2_fn,0,0,g.t1,g.t2,...
        str2double(settings.ldar_r),x0,y0,z0,1);
else
    NLDN2c=NaN(1,6);
    NLDN2g= NLDN2c;
end

a.NLDN2c = NLDN2c;
a.NLDN2g = NLDN2g;

% Do we need to time shift graphs
if settings.plot_tshiftOn
    % Time shifts according to LDAR points
    switch settings.plot_tshiftType
        case 1 % Plots were time shift to LDAR2
            xts = nanmean([CG(:,6);DLS(:,6)]);
            yts = nanmean([CG(:,7);DLS(:,7)]);
            zts = nanmean([CG(:,8);DLS(:,8)]);
        case 2 % Plots were time shift to PBFA
            xts = nanmean(PBFA(:,3));
            yts = nanmean(PBFA(:,4));
            zts = nanmean(PBFA(:,5));
        case 3 % plots were time shifted to linet
            xts = nanmean(LINET(:,6));
            yts = nanmean(LINET(:,7));
            zts = nanmean(LINET(:,8));
        case 4 % plots were time shifted to
            xts = nanmean(PBFA2(:,3));
            yts = nanmean(PBFA2(:,4));
            zts = nanmean(PBFA2(:,5));
        case 5 % plots were time shifted to
            xts = nanmean(PBFA_old(:,3));
            yts = nanmean(PBFA_old(:,4));
            zts = nanmean(PBFA_old(:,5));
        case 6 % plots were time shifted to NLDN2
            xts = nanmean([NLDN2g(:,2);NLDN2c(:,2)]);
            yts = nanmean([NLDN2g(:,3);NLDN2c(:,3)]);
            zts = nanmean([NLDN2g(:,4);NLDN2c(:,4)]);
    end
    
    if isnan(zts)
        zts = 0;
    end
    
   
    if isnan(xts)
        tshift1 = zeros(1,20);
    else
        tshift1 = sqrt((settings.x - xts).^2 + ...
                       (settings.y - yts).^2+ ...
                       (settings.z - zts).^2)/299792458.0;
    end
else
    tshift1 = zeros(1,20);
end

tDecay = [1 100e-6 10];   % ch3, ch1, and ch2 time decays

%% Loading and plotting LP data

% Store one data set for plot yy
for i=1:60
    % for loop used becuase there may be 15 lp plots
    if strcmp(a.lp_fn{i},'')==0
        [tn,vn]=SA_Extract1(a.lp_fn{i},g.t1,g.t2,tshift(i)-tshift1(ceil(i/3)));
        %soundsc(vn,10000,8,[-1 1])
        %[tn,vn,filtered] = lp_high_pass(tn,vn,200);
        
        %vn = remove_60Hz(tn,vn);
        
        if strcmp(settings.remove_time_decay,'on') 
            vn = removeDecay(vn,tDecay(mod(i,3)+1),tn(2)-tn(1));
            lp_legend{i} = sprintf('%s C', lp_legend{i});
        end
        
        %plot(tn,vn*gain(i)*factor(i)+vshift(i))
        a.lp_data(i).t = tn;
        a.lp_data(i).y = vn*gain(i)*factor(i)+vshift(i);
        
        empty_plot=false;
        try
            waitbar(0.25+i*0.005,wbh,'Loading slow antenna data','Name','Plotter Busy')
        catch
            return
        end
    end
end




%% Loading and plotting ch data

a.trigrd = [];
lg = [];

for i=1:60
    
    try
        waitbar(0.55+i*0.005,wbh,'Loading Fast Antenna data','Name','Plotter Busy')
    catch
        return
    end
    
    
    if strcmp(a.ch_fn{i},'')==0
        
        [t, y, ch_freq_str] = FA_Extract1(a.ch_fn{i},a.h_fn{i},g.t1,g.t2,tshift(i)-tshift1(ceil(i/3)),settings,i);
        
        if isempty(t)==0
            
            if isempty(a.trigrd)
                a.trigrd = sprintf('%i',ceil(i/3));
            else
                a.trigrd = sprintf('%s,%i',a.trigrd,ceil(i/3));
            end
              
            
            ch_legend_str = [ch_legend{i} ch_freq_str];
            
            
            % Do we need to remove time decay
            if strcmp(settings.remove_time_decay,'on')
                
                for i3 = 1:10000;
                    dttt = t(i3+1)-t(i3);
                    
                    if dttt > 0
                        break;
                    end
                end

                y = removeDecay(y,tDecay(mod(i,3)+1),dttt);
                ch_legend_str = [ch_legend{i} ch_freq_str ' C'];
            end
            
            % Do we need to filter the data
            if g.temp.ch_hpf > 0
                
                [t,y,filtered] = ch_high_pass(t,y,g.temp.ch_hpf);
                
                switch filtered
                    case 0; frintf('%s data has not filtered\n', ch_legend_str)
                    case 2; fprintf('%s data has not filtered. Consider increasing freq.\n', ch_legend_str)
                    case 3; fprintf('%s data has not properly filtered. Consider increasing freq.\n', ch_legend_str)
                end
            end
            
            % Reduce bit depty to eight
            %I = round(255*(y+5)/10);
            %y = 10*I/255 - 5;
            
            
            %plot(t,y*gain(i)*factor(i)+vshift(i))
            
            a.ch_data(i).y = y*gain(i)*factor(i)+vshift(i);
            a.ch_data(i).t = t;
                
            %soundsc(y,15000,8,[-1 1])
            lg=[lg ch_legend_str];
            ch_legend{i} = ch_legend_str;
            empty_plot=false;
        end
    else
        a.ch_data(i).y = [];
        a.ch_data(i).t = [];
    end
end

%% Loading and plotting dE/dt data
if ~strcmp(a.dedt_fn,'')
    ch_number = 3;
    [dedt_t,dedt_y] = dedtExtract(dedt_fn,g.t1,g.t2,ch_number);
    
    if ~isnan(dedt_t(1))
        plot(dedt_t,dedt_y)
        lg=[lg 'dE/dt'];
        empty_plot = false;
    end
end

%% Loading and plotting dE/dt trig data
if ~strcmp(a.dedt_trig_fn,'')
    ch_number = 2;
    [dedt_t,dedt_y] = dedtExtract(a.dedt_trig_fn,g.t1,g.t2,ch_number);
    
    if ~isnan(dedt_t(1))
        plot(dedt_t,dedt_y)
        lg=[lg 'dE/dt Trig'];
        empty_plot = false;
    end
end


%% Loading and plotting feild mill data

% Read time shift data file
fm_tshifts = csvread('fm_time_shifts.csv');
ind = find(fm_tshifts(:,1) == str2double([g.YYYY{:} g.MM{:} g.DD{:}]));

if isempty(ind)
    fm_tshifts = zeros(1,34);
else
    fm_tshifts = fm_tshifts(ind,2:end);
end

for i=1:34

    if strcmp(a.fm_fn{i},'')==0
        [tn,vn]=fieldMillExtract3(a.fm_fn{i},g.t1,g.t2,fm_tshifts(i));
        
        if ~isempty(tn)
            lg=[lg sprintf('FM-%2.2i',i)];       
            %plot(tn,-vn)
            a.fm_data(i).t = tn;
            a.fm_data(i).y = -vn;
        
            empty_plot=false;
        end
        
        try
            waitbar(0.89,wbh,'Loading FM data')
        catch
            return
        end
    end
end

a.lp_legend = lp_legend;
a.ch_legend = ch_legend;

delete(wbh);



