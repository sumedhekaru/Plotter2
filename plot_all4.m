%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot all3(g)
%
%   This function writtent to specially work with plotter GUI
%   The necessary input argument g contains all the controll to the
%   plotter.
%
%   Modification History
%       2010-05-13 Created by Sumedhe Karunarathne
%       2011-07-04 Renamed as plotter 2 and add ability to plot 20 sensors
%       data
%       2011-07-04 Time shift info will dispaly if it is used. Voltage info
%       will not dilplayed
%       2011-08-31 dedt data and dedt_trig data will be plotted
%
%       % Changed between plot_all3 - can handle LINET2 from NLDN button
%                 plot_all4 - can handle PBFA2 (PBFA from positive peaks)
%                 from NLDN button
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = plot_all4(g,verbose)
%clc
if nargin < 2
    % turn on verbose default (meaning error massages will be turn on
    verbose = 1;
end
   

% creating struct to store things
a=struct;

global AX

%LDAR LINET PBFA marker size
mz = 4;


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


% Do we need to time shift graphs to points?
if ~settings.plot_tshiftOn
    settings.plot_tshiftType == 0;
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

% the variable for real legend
lg={};



% is plot empty
empty_plot=true;

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

%% Generating file name for LDAR data
if g.ldar==1 || g.cglss == 1 || (settings.plot_tshiftOn && settings.plot_tshiftType == 1)
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
        if g.ldar==1 || g.cglss == 1
            a.absant_fn=[a.absant_fn ldar_fn];
        end
    end
else
    a.ldar_fn='';
end

%% Generating file name for LINET

if g.linet==1 || (settings.plot_tshiftOn && settings.plot_tshiftType == 3)
    
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
    
    linet_fn=sprintf('%s/LINET/%s/%s/%s/linet_%s%s%s_%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
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
    
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
    
    nldn_fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
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
   
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
    
    pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
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
   
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
    
    pbfa_old_fn=sprintf('%s/PBFA_old/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
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
    
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
     
    dedt_fn=sprintf('%s/dedt/%s/%s/%s/%2.2i%2.2i/BCC_%s%s%s_%2.2i%2.2i.hh3',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.hh,ext,g.YYYY{:},...
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
    
    if g.mm < 30
        ext=0;
    else
        ext=30;
    end
     
    dedt_trig_fn=sprintf('%s/dedt/%s/%s/%s/%2.2i%2.2i/BCC_%s%s%s_%2.2i%2.2i.hh2',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.hh,ext,g.YYYY{:},...
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


%% Letting user know about missing file

a.absant_fn=sort(a.absant_fn);

if isempty (a.absant_fn)==0 && verbose
    answer=questdlg(a.absant_fn,'Files not found!','OK','Stop!','OK');
    if strcmp(answer,'Stop!')
        return
    end
end


%% If there is nothing to plot Exit
if sum(strcmp(a.lp_fn,''))==60 && ...
        sum(strcmp(a.ch_fn,''))== 60 && strcmp(a.ldar_fn,'')== 1 && ...
        strcmp(a.linet_fn,'')==1     && strcmp(a.pbfa_fn,'')== 1 && ...
        strcmp(a.dedt_fn,'') ==1     && strcmp(a.dedt_trig_fn,'')== 1
    if verbose
        errordlg('No data files were found. All Plot commads failed!', ...
            'Plotter Error','model')
    end
    return
end

wbh= waitbar(0,'Please wait...','Name','Plotter Busy');

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

if strcmp(a.linet_fn,'')==0
    ld_lin=[ld_lin '-LINET-'];
    try
        waitbar(0.15,wbh,'Loading LINET data','Name','Plotter Busy')
    catch
        return
    end
    LINET=linetExtract(linet_fn,g.t1,g.t2,str2double(settings.ldar_r),...
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


if strcmp(a.pbfa_old_fn,'')==0
    ld_lin=[ld_lin '-PBFAO-'];
    try
        waitbar(0.25,wbh,'Loading PBFA-O data','Name','Plotter Busy')
    catch
        return
    end
    PBFA_old=pbfaExtract(a.pbfa_old_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0);
else
    PBFA_old=NaN(1,6);    
end


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


%% Loading and plotting LP data

fg=figure;
set(fg,'visible','off')
hold all

% Store one data set for plot yy

for i=1:60
    % for loop used becuase there may be 15 lp plots
    if strcmp(a.lp_fn{i},'')==0
        [tn,vn]=SA_Extract1(a.lp_fn{i},g.t1,g.t2,tshift(i)-tshift1(ceil(i/3)));
        %soundsc(vn,10000,8,[-1 1])
        %[tn,vn,filtered] = lp_high_pass(tn,vn,200);
        
        plot(tn,vn*gain(i)*factor(i)+vshift(i))
        lg=[lg lp_legend{i}];
        empty_plot=false;
        try
            waitbar(0.25+i*0.005,wbh,'Loading slow antenna data','Name','Plotter Busy')
        catch
            return
        end
    end
end




%% Loading and plotting ch data

output.trigrd = [];

for i=1:60
    
    try
        waitbar(0.55+i*0.005,wbh,'Loading Fast Antenna data','Name','Plotter Busy')
    catch
        return
    end
    
    
    if strcmp(a.ch_fn{i},'')==0
        
        [t y ch_freq_str] = FA_Extract1(a.ch_fn{i},a.h_fn{i},g.t1,g.t2,tshift(i)-tshift1(ceil(i/3)),settings,i);
        
        if isempty(t)==0
            
            if isempty(output.trigrd)
                output.trigrd = sprintf('%i',ceil(i/3));
            else
                output.trigrd = sprintf('%s,%i',output.trigrd,ceil(i/3));
            end
              
            
            ch_legend_str = [ch_legend{i} ch_freq_str];
            lg=[lg ch_legend_str];
            
            % Do we need to filter the data
            if g.temp.ch_hpf > 0
                
                [t,y,filtered] = ch_high_pass(t,y,g.temp.ch_hpf);
                
                switch filtered
                    case 0; frintf('%s data has not filtered\n', ch_legend_str)
                    case 2; fprintf('%s data has not filtered. Consider increasing freq.\n', ch_legend_str)
                    case 3; fprintf('%s data has not properly filtered. Consider increasing freq.\n', ch_legend_str)
                end
            end
            
            plot(t,y*gain(i)*factor(i)+vshift(i))
            %soundsc(y,15000,8,[-1 1])

            empty_plot=false;
        end
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

%%
if settings.plot_calibrated == 1
    ylabel('E (V/m)')
else
    ylabel('Voltage (V)')
end


%% %%%%%%%%%%%%%%   PLOTTING LINET / LDAR / PBFA %%%%%%%%%%%%%%%%%%%%%%%%%


if ~isnan(DLS(1,10)) || ~isnan(CG(1,10)) || ~isnan(LINET(1,3)) ...
        || ~isnan(PBFA(1,6)) || ~isnan(PBFA2(1,6)) || ~isnan(PBFA_old(1,6))
    try
        waitbar(0.90,wbh,'Plotting graphs','Name','Plotter Busy')
    catch
        return
    end
    
    [AX,H1,H2]=plotyy(nan,nan,nan,nan);
    hold(AX(2), 'on')
    sec_y_on = false;
    
    if sum(~isnan(DLS(:,10))) > 0 && g.ldar
        plot (AX(2),DLS(:,10),DLS(:,8),'ko','MarkerFaceColor','k','MarkerSize',mz)
        lg = [lg 'LDAR2'  ];
        sec_y_on = true;
        empty_plot=false;
    end
    
    if sum(~isnan(CG(:,10)))>0 && g.cglss
        plot (AX(2),CG(:,10),CG(:,8),'gs','MarkerFaceColor','g','MarkerSize',mz+2)
        lg = [lg 'CGLSS'  ];
        empty_plot=false;
    end
    
    if sum(~isnan(LINET(:,3)))>0 && g.linet
        plot (AX(2),LINET(:,3),LINET(:,8),'co','MarkerFaceColor','c','MarkerSize',mz)
        lg = [lg 'LINET'  ];
        sec_y_on = true;
        empty_plot=false;
    end
    
    %     if sum(~isnan(PBFA2(:,6)))>0 && g.nldn
    %         plot (AX(2),PBFA2(:,6),PBFA2(:,5),'bo','MarkerFaceColor','b','MarkerSize',mz)
    %         lg = [lg 'PBFA-A'  ];
    %         sec_y_on = true;
    %         empty_plot=false;
    %     end
    if sum(~isnan(PBFA(:,6)))>0 && g.pbfa
        if settings.pbfa_ebar_on
            errorbar(AX(2),PBFA(:,6),PBFA(:,5),PBFA(:,12),'ro','MarkerFaceColor','r','MarkerSize',mz)
        else
            plot (AX(2),PBFA(:,6),PBFA(:,5),'ro','MarkerFaceColor','r','MarkerSize',mz)
        end
        
        lg = [lg 'PBFA'  ];
        sec_y_on = true;
        empty_plot=false;
    end


    if sum(~isnan(PBFA2(:,6)))>0 && g.nldn
        if settings.pbfa_ebar_on
            errorbar(AX(2),PBFA2(:,6),PBFA2(:,5),PBFA2(:,12),'bo','MarkerFaceColor','b','MarkerSize',mz)
        else
            plot (AX(2),PBFA2(:,6),PBFA2(:,5),'bo','MarkerFaceColor','b','MarkerSize',mz)
        end
        
        lg = [lg 'PBFA-A'  ];        
        sec_y_on = true;
        empty_plot=false;
    end
    
    
    if sum(~isnan(PBFA_old(:,6)))>0 && g.pbfa_old
        if settings.pbfa_ebar_on
            errorbar(AX(2),PBFA_old(:,6),PBFA_old(:,5),PBFA_old(:,12),'bs','MarkerFaceColor','r','MarkerSize',mz)
        else
            plot (AX(2),PBFA_old(:,6),PBFA_old(:,5),'bs','MarkerFaceColor','b','MarkerSize',mz)
        end
        
        lg = [lg 'PBFA-O'  ];
        sec_y_on = true;
        empty_plot=false;
    end
    
    
    if sec_y_on
        
        set(get(AX(2),'Ylabel'),'String','Altitude (m)')
        mmm  = settings.max_z;
        mmm2 = settings.min_z;
        yinc = (mmm-mmm2)/5; 
        set(AX(2),'yLim',[mmm2 mmm],'YTick',round(mmm2:yinc:mmm),...
            'YTickLabel',num2str(round(mmm2:yinc:mmm)'));          

    else
        set(AX(2),'YTick',[])        
    end
        
        
    set(AX,'xlim',[g.t1 g.t2])
    linkaxes(AX,'x')
    % set(AX(1),'ylimmode','auto')
    % mmm=ceil(max([y1;y2])/1000)*1000;
        

    set(AX(2),'xtick',[])
    set(AX, 'YColor', [0 0 0])
    
    
    %xlimits = get(ax1,'XLim');
    ylimits = get(AX(1),'YLim');
    %xinc = (xlimits(2)-xlimits(1))/5;
    yinc = (ylimits(2)-ylimits(1))/5;
    
    set(AX(1),'YTick',ylimits(1):yinc:ylimits(2))
    
   
    
    %lg=[lg 'LDAR'];
else
    AX=get(gcf,'CurrentAxes');
end

try
    waitbar(0.95,wbh,'All plotting finished!','Name','Plotter Busy')
catch
    return
end

% Create additional tool menu in the plot window

try
    waitbar(0.97,wbh,'Creating Extra Plotting Tools','Name','Plotter Busy')
catch
    return
end

% Add aditional toolls to figure
tools2fig;


try
    
    waitbar(0.99,wbh,'Final Preparation','Name','Plotter Busy')
catch
    return
end

if empty_plot==false
    box on
    grid on
%     title_str=sprintf('%s-%s-%s    %s   UT: %2.2i:%2.2i:%2.2i ', ...
%         g.YYYY{:},g.MM{:},g.DD{:},ld_lin,g.hh,g.mm,g.ss);

        title_str=sprintf('%s-%s-%s     UT: %2.2i:%2.2i:%2.2i ', ...
        g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss);
    
         
    if ~isempty(tshift_str)
        %title_str = [title_str  tshift_str];
        xlabel(['Time (s)   ' tshift_str])
    else
        xlabel('Time (s)')    
    end
    
    if g.temp.ch_hpf ~= 0
        title_str = [title_str  sprintf('   Filter = %i Hz',g.temp.ch_hpf)];
    end
    
    if settings.ldar_tshiftOn
        title_str=[title_str '  LDAR:' settings.sen_IDs{sn}];
    elseif settings.plot_tshiftOn
        ts_str = {'LDAR' 'PBFA' 'LINET' 'PBFA-A' 'PBFA-O'};
        title_str=[title_str '  Plots:' ts_str{settings.plot_tshiftType}];        
    end
    
    try
        waitbar(0.99,wbh,'Final Preparation','Name','Plotter Busy')
    catch
        return
    end
    
    
    if isempty(lg)==0
        legend(lg)
    end
    
    title(title_str)    
        
    
    try
        waitbar(1,wbh,'Final Preparation','Name','Plotter Busy')
    catch
        return
    end
    % set(AX,'PlotBoxAspectRatio',[1 1 1])
    set(AX,'xlim',[g.t1 g.t2])
    
    % Delete wait bar just before plot is visible
    try
        delete(wbh)
    catch
        % Do nothing
    end
    set(fg,'visible','on')
    plotter_tools(15)
    
else
    try
        delete(wbh)
    catch
        % Do nothing
    end
    if verbose
        errordlg('Plot was empty. May be no data in the given time range!', ...
            'Plotter Error','model')
    end
    
end



