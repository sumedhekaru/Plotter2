%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot all(g)
% 
%   This function writtent to specially work with plotter GUI
%   The necessary input argument g contains all the controll to the
%   plotter.
%
%   Modification History
%       2010-05-13 Created by Sumedhe Karunarathne
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_all_hilbert(g)
% clc

% creating struct to store things
a=struct;

global AX



%% Creating the directory path for all files
% Opening sensor settings
settings=open('sensor_setting.mat');
% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
                settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

%Is manual Time shift for each sensor activated
if settings.man_tshiftOn==1
    ts=settings.t_shift;
    tshift=[ts(1),ts(1),ts(1),ts(2),ts(2),ts(2),ts(3),ts(3),ts(3),...
        ts(4),ts(4),ts(4),ts(5),ts(5),ts(5)];
else
    tshift=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
end

%Is offset voltages should be included?
if settings.vshiftOn==1
    vshift=settings.vshift;
else
    vshift=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
end

%Is manual gain should be included?
if settings.gainOn==1
    gain=settings.gain;
else
    gain=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
end

% Generating LP Legend names
lp_legend={[settings.sen_IDs{1} ':lp1'],[settings.sen_IDs{1} ':lp2'],[settings.sen_IDs{1} ':lp3'],...
    [settings.sen_IDs{2} ':lp1'],[settings.sen_IDs{2} ':lp2'],[settings.sen_IDs{2} ':lp3'],...
    [settings.sen_IDs{3} ':lp1'],[settings.sen_IDs{3} ':lp2'],[settings.sen_IDs{3} ':lp3'],...
    [settings.sen_IDs{4} ':lp1'],[settings.sen_IDs{4} ':lp2'],[settings.sen_IDs{4} ':lp3'],...
    [settings.sen_IDs{5} ':lp1'],[settings.sen_IDs{5} ':lp2'],[settings.sen_IDs{5} ':lp3']};

% Generating Ch Legend names
ch_legend={[settings.sen_IDs{1} ':ch1'],[settings.sen_IDs{1} ':ch2'],[settings.sen_IDs{1} ':ch3'],...
    [settings.sen_IDs{2} ':ch1'],[settings.sen_IDs{2} ':ch2'],[settings.sen_IDs{2} ':ch3'],...
    [settings.sen_IDs{3} ':ch1'],[settings.sen_IDs{3} ':ch2'],[settings.sen_IDs{3} ':ch3'],...
    [settings.sen_IDs{4} ':ch1'],[settings.sen_IDs{4} ':ch2'],[settings.sen_IDs{4} ':ch3'],...
    [settings.sen_IDs{5} ':ch1'],[settings.sen_IDs{5} ':ch2'],[settings.sen_IDs{5} ':ch3']};

% the variable for real legend
lg={};

% is plot empty
empty_plot=true;

%% Generating file names for lp graphs and check those files are availble

% checking witch graphs are on
lp_on=g.lpgraphs;

% Save absant file extentions
a.absant_fn={};

for i=1:15
    if lp_on(i)==1
        % finding the file extention number
        ext=mod(i,3);
        if ext==0
            ext=3;
        end
        
        % Finding the stattion ID
        if i <=3
            sid=settings.sen_IDs{1};
        elseif i <= 6
            sid=settings.sen_IDs{2};
        elseif i <= 9
            sid=settings.sen_IDs{3};
        elseif i <= 12
            sid=settings.sen_IDs{4};
        else
            sid=settings.sen_IDs{5};
        end
    
        
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

for i=1:15
    if ch_on(i)==1
        % finding the file extention number
        ext=mod(i,3);
        if ext==0
            ext=3;
        end
        
        % Finding the stattion ID
        if i <=3
            sid=settings.sen_IDs{1};
        elseif i <= 6
            sid=settings.sen_IDs{2};
        elseif i <= 9
            sid=settings.sen_IDs{3};
        elseif i <= 12
            sid=settings.sen_IDs{4};
        else
            sid=settings.sen_IDs{5};
        end
        
        
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
if g.ldar==1
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
else
    a.ldar_fn='';
end

%% Generating file name for LINET
% Did not completed yet!!!

if g.linet==1
    
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
        a.absant_fn=[a.absant_fn linet_fn];
    end
else
    a.linet_fn='';
end


%% Generating file name for PBFA
% Did not completed yet!!!

if g.pbfa==1
    
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
        a.absant_fn=[a.absant_fn pbfa_fn];
    end
else
    a.pbfa_fn='';
end


%% Letting user know about missing file

a.absant_fn=sort(a.absant_fn);

if isempty (a.absant_fn)==0
    answer=questdlg(a.absant_fn,'Files not found!','OK','Stop!','OK'); 
    if strcmp(answer,'Stop!')
        return
    end
end


%% If there is nothing to plot Exit
if sum(strcmp(a.lp_fn,''))==15 && ...
    sum(strcmp(a.ch_fn,''))== 15 && ...
        strcmp(a.ldar_fn,'')== 1 && ...
            strcmp(a.linet_fn,'')==1 && ...
                strcmp(a.pbfa_fn,'')==1
  errordlg('No data files were found. All Plot commads failed!', ...
      'Plotter Error','model')
    return
end



%% Calculating time corrections according to sensor positions
% Loading position setting file
%  try
%      b=open('sensor_setting.mat');
%  catch
%      answer=questdlg('Sensor Position Settings were Not Found!','Error','OK','Stop!','OK');
%      if strcmp(answer,'Stop!')
%          return/home/daqop
%      else
%          tCorrection(1:4)=0;
%      end
%  end 
     

%% Loading and plotting LP data
% Line color
% gc=['r','r','r',...
%     'g','g','g',...
%     'b','b','b',...
%     'k','k','k',.../home/daqop
%     'c','c','c'];

wbh= waitbar(0,'Please wait...','Name','Plotter Busy');

fg=figure;
subplot(2,2,1);
set(fg,'visible','off')
hold all

% Store one data set for plot yy
ldar_t=[];
ldar_v=[];

for i=1:15
% for loop used becuase there may be 15 lp plots
    if strcmp(a.lp_fn{i},'')==0
        [tn,vn]=SA_Extract1(a.lp_fn{i},g.t1,g.t2,tshift(i));
        ldar_t=tn;
        ldar_v=vn;
        subplot(2,2,1)
        plot(tn,vn*gain(i)+vshift(i))
        
%         % Hilbert transform of voltages
%         vn_hil = hilbert(vn*gain(i)+vshift(i));
%         
%         % Plotting Real Part
%         subplot(2,2,2)
%         plot(tn,real(vn_hil))
%         title('Real')
%         
%         
%         % Plotting Imaginary Part
%         subplot(2,2,3)
%         plot(tn,imag(vn_hil))
%         title('Imaginary')
%         
%          % Plotting Imaginary Part
%         subplot(2,2,4)
%         plot(tn,imag(vn_hil).^2+real(vn_hil).^2)
%         title('Magnitude')        
        
        lg=[lg lp_legend{i}];
        empty_plot=false;
        
        waitbar(i*0.02,wbh,'Loading slow antenna data','Name','Plotter Busy')
    end    
end




%% Loading and plotting ch data

t=[];   % Variable for time
y=[];   % Variable for voltage
empty_trigs=[]; % Variable for empty triggers

for i=1:15
    
    waitbar(0.3+i*0.02,wbh,'Loading Fast Antenna data','Name','Plotter Busy')
    
    t=[];   % Variable for time
    y=[];   % Variable for voltage
        
    % Plotting up to 15 fast antenna files   
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
            % Keep a track of empty trigs
            if isempty(sFa.y_i)
                empty_trigs=[empty_trigs,trigs(j)]
            else
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
        if isempty(t)==0
            ldar_t=t;
            ldar_v=y;
            
            subplot(2,2,1)
            hold all
            plot(t,y*gain(i)+vshift(i))
            title('Real Data')
            grid on
            box on
            
            
            % filter out low frequencies
            Fs = 1/(t(2)-t(1));
            [z,p] = butter(5,1000/(Fs/2),'high'); % Create a low-pass butterworth filter;
            % [z,p,k] = butter(n,Wn) designs an designs an order n lowpass digital Butterworth filter with normalized
            % cutoff frequency Wn. It returns the zeros and poles in length n column
            % vectors z and p, and the gain in the scalar k
            
            smoothy = filtfilt(z,p,y*gain(i));    % filter the data.
                        
            % Hilbert transform of voltages
            y_hil = hilbert(smoothy);
                        
            % Plotting Real Part
            subplot(2,2,2)
            hold all
            plot(t,real(y_hil))
            box on
            grid on
            title('Real Part of Hilbert Transformation')
          
            % Plotting Imaginary Part
            subplot(2,2,3)
            hold all
            plot(t,imag(y_hil))
            title('Imaginary part of Hilbert Transformation')
            box on
            grid on
                        
            % Plotting Magnitude
            subplot(2,2,4)
            hold all
            plot(t,sqrt(imag(y_hil).^2+real(y_hil).^2)) 
            box on
            grid on
            
            lg=[lg ch_legend{i}];
            empty_plot=false;
        end
    end
end

subplot(2,2,1)

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

%% %%%%%%%%%%%%%%   PLOTTING LINET / LDAR / PBFA %%%%%%%%%%%%%%%%%%%%%%%%%
ld_lin='';
x1=[];
y1=[];

if strcmp(a.ldar_fn,'')==0
    ld_lin=[ld_lin '-LDAR-'];
    waitbar(0.7,wbh,'Loading LDAR data','Name','Plotter Busy')
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0,0);
    
    x1=[DLS(:,10);CG(:,10)];
    y1=[DLS(:,8);CG(:,8)];

end

if strcmp(a.linet_fn,'')==0
    ld_lin=[ld_lin '-LINET-'];
    waitbar(0.73,wbh,'Loading LINET data','Name','Plotter Busy')
    LINET=linetExtract(linet_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0);
    % LINET Column DISCRIPTION
    % 1 - time
    % 2 - Distance from the time correcting sensor
    % 3 - time shifted time
    % 6 - x distance
    % 7 - y distance
    % 8 - z distance
    x1=[x1; LINET(:,3)];
    y1=[y1; LINET(:,8)];

end

if strcmp(a.pbfa_fn,'')==0
    ld_lin=[ld_lin '-PBFA-'];
    waitbar(0.76,wbh,'Loading PBFA data','Name','Plotter Busy')
    PBFA=pbfaExtract(pbfa_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0);
    
    % PBFA column Discription
    %   1 - distance to the pulse from the given sensor     
    %   2 - Occuring time
    %   3 - x
    %   4 - y
    %   5 - z
    %   6 - Detection time at the sensor
    
    
    x1=[x1; PBFA(:,6)];
    y1=[y1; PBFA(:,5)];
 
end





if isempty(x1)==0
    waitbar(0.78,wbh,'Plotting graphs','Name','Plotter Busy')
    if isempty (ldar_t)==0
        [AX,H1,H2]=plotyy(ldar_t,ldar_v,x1,y1);
        empty_plot=false;
    else
        ldar_t=x1;
        ldar_v=x1-x1;
        [AX,H1,H2]=plotyy(ldar_t,ldar_v,x1,y1);
        empty_plot=false;
    end
    set(H1,'LineStyle','none')
    set(H2,'LineStyle','none')
    set(H2,'Marker','p','MarkerEdgeColor','g','MarkerFaceColor','g')
    set(get(AX(1),'Ylabel'),'String','Voltage (V)')
    set(get(AX(2),'Ylabel'),'String','Altitude (m)')
    set(AX,'xlim',[g.t1 g.t2])
    linkaxes(AX,'x')
    set(AX(1),'ylimmode','auto')
    %mmm=ceil(max([y1;y2])/1000)*1000;
    mmm=12500;
    set(AX(2),'yLim',[0 mmm],'YTick',0:2000:mmm);
    set(AX(2),'xtick',[])
    set(AX, 'YColor', [0 0 0])
    
    
    %xlimits = get(ax1,'XLim');
    ylimits = get(AX(1),'YLim');
    %xinc = (xlimits(2)-xlimits(1))/5;
    yinc = (ylimits(2)-ylimits(1))/5;
    
    set(AX(1),'YTick',ylimits(1):yinc:ylimits(2))
    
    ylimits = get(AX(2),'YLim');
    %xinc = (xlimits(2)-xlimits(1))/5;
    yinc = (ylimits(2)-ylimits(1))/5;
    
    set(AX(2),'YTick',ylimits(1):yinc:ylimits(2))
    
    %#############################################
    subplot(2,2,4)
    if isempty (ldar_t)==0
        [AX,H1,H2]=plotyy(ldar_t,ldar_v,x1,y1);
        empty_plot=false;
    else
        ldar_t=x1;
        ldar_v=x1-x1;
        [AX,H1,H2]=plotyy(ldar_t,ldar_v,x1,y1);
        empty_plot=false;
    end
    set(H1,'LineStyle','none')
    set(H2,'LineStyle','none')
    set(H2,'Marker','p','MarkerEdgeColor','g','MarkerFaceColor','g')
    set(get(AX(1),'Ylabel'),'String','Voltage (V)')
    set(get(AX(2),'Ylabel'),'String','Altitude (m)')
    set(AX,'xlim',[g.t1 g.t2])
    linkaxes(AX,'x')
    set(AX(1),'ylimmode','auto')
    %mmm=ceil(max([y1;y2])/1000)*1000;
    mmm=12500;
    set(AX(2),'yLim',[0 mmm],'YTick',0:2000:mmm);
    set(AX(2),'xtick',[])
    set(AX, 'YColor', [0 0 0])
    
    
    %xlimits = get(ax1,'XLim');
    ylimits = get(AX(1),'YLim');
    %xinc = (xlimits(2)-xlimits(1))/5;
    yinc = (ylimits(2)-ylimits(1))/5;
    
    set(AX(1),'YTick',ylimits(1):yinc:ylimits(2))
    
    ylimits = get(AX(2),'YLim');
    %xinc = (xlimits(2)-xlimits(1))/5;
    yinc = (ylimits(2)-ylimits(1))/5;
    
    set(AX(2),'YTick',ylimits(1):yinc:ylimits(2))
    
    %#############################################
    
    %lg=[lg 'LDAR'];
else
    AX=get(gcf,'CurrentAxes');
end

waitbar(0.8,wbh,'All plotting finished!','Name','Plotter Busy')


% Create additional tool menu in the plot window

waitbar(0.85,wbh,'Creating Extra Plotting Tools','Name','Plotter Busy')

f = uimenu('Label','Plotter');
    uimenu(f,'Label','Update Y Grids','Callback','plotter_tools(1)','Accelerator','Y');
    uimenu(f,'Label','Update X Grids','Callback','plotter_tools(2)','Accelerator','X');    
    uimenu(f,'Label','Update Time Range','Callback','plotter_tools(3)','Accelerator','T');
    uimenu(f,'Label','Find Delta t','Callback','plotter_tools(4)','Accelerator','D');
    uimenu(f,'Label','Find Delta y','Callback','plotter_tools(10)','Accelerator','F');
    uimenu(f,'Label','hh:mm:ss Data Curser','Callback','plotter_tools(5)','Accelerator','K');
    uimenu(f,'Label','ss.ssssss Data Curser','Callback','plotter_tools(6)','Accelerator','M');
    uimenu(f,'Label','X Zoom!','Callback','plotter_tools(7)','Accelerator','G');
    uimenu(f,'Label','Y Zoom!','Callback','plotter_tools(8)','Accelerator','H');
    uimenu(f,'Label','XY Zoom!','Callback','plotter_tools(9)','Accelerator','J');
    uimenu(f,'Label','Find Position V5','Callback','location_v5');
    uimenu(f,'Label','Find Position V7','Callback','location_v7','Accelerator','L');
    uimenu(f,'Label','Pulse Modeling','Callback','pulse4','Accelerator','B');
    uimenu(f,'Label','Vedio Framing','callback','plotter_tools(11)');
    uimenu(f,'Label','Fix for Publishing','callback','fix_fig');





waitbar(0.9,wbh,'Final Preparation','Name','Plotter Busy')

if empty_plot==false
    box on
    grid on
    title_str=sprintf('%s-%s-%s    %s   UT: %2.2i:%2.2i:%2.2i  \nTime Shifts: %0.1f -- %0.1f -- %0.1f -- %0.1f -- %0.1f uS', ...
        g.YYYY{:},g.MM{:},g.DD{:},ld_lin,g.hh,g.mm,g.ss,tshift(1)*1e6,...
        tshift(4)*1e6,tshift(7)*1e6,tshift(10)*1e6,tshift(13)*1e6);
   
    if settings.ldar_tshiftOn==1
        title_str=[title_str '  ' settings.sen_IDs{sn}];
    end
    waitbar(0.95,wbh,'Final Preparation','Name','Plotter Busy')
    
%     title_str=sprintf('%s \n V shifts : %.2f:%.2f:%.2f -:- %.2f:%.2f:%.2f -:- %.2f:%.2f:%.2f -:- %.2f:%.2f:%.2f -:- %.2f:%.2f:%.2f',...
%         title_str,vshift(1),vshift(2),vshift(3),vshift(4),vshift(5),vshift(6),vshift(7),vshift(8),vshift(9),vshift(10),vshift(11),vshift(12),vshift(13),vshift(14),vshift(15)); 
%     
%     title_str=sprintf('%s \n Gains : %.1f:%.1f:%.1f -:- %.1f:%.1f:%.1f -:- %.1f:%.1f:%.1f -:- %.1f:%.1f:%.1f -:- %.1f:%.1f:%.1f',...
%         title_str,gain(1),gain(2),gain(3),gain(4),gain(5),gain(6),gain(7),gain(8),gain(9),gain(10),gain(11),gain(12),gain(13),gain(14),gain(15)); 
%     
    if isempty(lg)==0
        legend(lg)       
    end
    title(title_str)
    
    xlab=sprintf('Time (s) \n V shifts : %.2f:%.2f:%.2f -:- %.2f:%.2f:%.2f -:- %.2f:%.2f:%.2f -:- %.2f:%.2f:%.2f -:- %.2f:%.2f:%.2f',...
        vshift(1),vshift(2),vshift(3),vshift(4),vshift(5),vshift(6),vshift(7),vshift(8),vshift(9),vshift(10),vshift(11),vshift(12),vshift(13),vshift(14),vshift(15)); 
    xlab=sprintf('%s \n Gains : %.1f:%.1f:%.1f -:- %.1f:%.1f:%.1f -:- %.1f:%.1f:%.1f -:- %.1f:%.1f:%.1f -:- %.1f:%.1f:%.1f',...
        xlab,gain(1),gain(2),gain(3),gain(4),gain(5),gain(6),gain(7),gain(8),gain(9),gain(10),gain(11),gain(12),gain(13),gain(14),gain(15)); 
    
    xlabel(xlab)
    
    waitbar(1,wbh,'Final Preparation','Name','Plotter Busy')
    % set(AX,'PlotBoxAspectRatio',[1 1 1])   
    set(AX,'xlim',[g.t1 g.t2])
    
    % Delete wait bar just before plot is visible
    delete(wbh)
    
    set(fg,'visible','on')
    
else
    delete(wbh)
    errordlg('Plot was empty. May be no data in the given time range!', ...
      'Plotter Error','model')
  
end




function fix_yticks()
   
    ylimits = get(AX1,'YLim');    
    yinc = (ylimits(2)-ylimits(1))/5;    
    set(AX1,'YTick',ylimits(1):yinc:ylimits(2))
    
    ylimits = get(AX2,'YLim');    
    yinc = (ylimits(2)-ylimits(1))/5;
    
    set(AX2,'YTick',ylimits(1):yinc:ylimits(2))
   