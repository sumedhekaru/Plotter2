function NBP_finder
% This function is inteneded to find NBP using fast antenna data.

%% User inputs
t1 = 0;
t2 = 86400;
sn = 9;
ch = 1;
pre = 0.2;
post = 0.2;
% Base folder to save stuffs
bf = 'C:\Users\Sumedhe\Desktop\NBP\';

%% Setup
count = 0;
% get Plotter2 data
try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

g=h.g;
settings = h.sen_set;
date = sprintf('%s%s%s',g.YYYY{:},g.MM{:},g.DD{:});


% save data in txt file
sv_fn = sprintf('%s/%s/FLT/%s-NBP_info.txt',bf,date,date)

% Create file if not exist and add the header
% Write the header if not exit
if exist(sv_fn,'file')
    fID = fopen(sv_fn,'a+');
else
    fID = fopen(sv_fn,'a+');
    fprintf(fID,'count\tt\tt\tLDAR-x\tLDAR-y\tLDAR-z\tLDAR-r\tTriggered\n');
end

% less convert t1 in to 5mins file
t1_5 = floor(t1/300)*300;
t2_5 = floor(t2/300)*300;

% turn off all the plots
g.lpgraphs = zeros(1,60);
g.chgraphs = zeros(1,60);
ind0 = (sn-1)*3+ch;
g.chgraphs(ind0) = 1;
g.ldar = 1;
g.cglss = 1;
g.pbfa = 1;
g.nldn = 1;


tshift = settings.t_shift(ind0);


while t1_5 <= t2_5 - 300;
    
    g.t1 = t1_5;
    g.t2 = t1_5 + 300;
    [hhmmss hh mm] = sec2hhmmss(g.t1);
    g.hh = hh;
    g.mm= floor(mm/5)*5;
    g.ss = 0;
    
    % generate file names
    [ch_fn h_fn] = generate_ch_fn(settings,g);
    
    % Get triggers
    if ~isempty(h_fn{ind0})
        trigs = get_triggers(h_fn{ind0},tshift,g.t1,g.t2);
        L = length(trigs);
        
        % check each trigger
        for i = 9:L
            t1 = trigs(i);
            t2 = trigs(i) + pre + post;
            
            [t,y]=FA_Extract1(ch_fn{ind0}, h_fn{ind0} ,t1,...
                t2,tshift,settings,sn);
            
            if ~isempty(t)
                
                % Apply a filter to get rid of DC
                [t,yf] = ch_high_pass(t,y,3000);
                
                m1 = max(yf);
                m2 = abs(min(yf));
                
                if m1 > abs(m2)
                    disp('Positive')
                    peakLoc = peakfinder(yf,m1/4,m1/4,1);
                else
                    disp('Negative')
                    peakLoc = peakfinder(yf,m2/4,m2/4,-1);
                end
                
                
                L1 = length(peakLoc);
                disp(L1)
                
                
                if L1 == 1
                    count = count + 1;
                    t0 = t(peakLoc);
                    [hhmmss hh,mm,ss,ss_dec] = sec2hhmmss(t0);
                    
                    g.t1 = t1;
                    g.t2 = t2;
                    
                    % Closest ldar2 point
                    ldar_fn = generate_ldar_fn(settings,g);
                    [CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),0,0,0,0);
                    
                    ts=DLS(:,10);         % time
                    xmin = nan;
                    ymin = nan;
                    zmin = nan;
                    rmin = nan;
                    
                    if ~isempty(ts)
                        % find closest ldar point to time "t"
                        [minVal,minInd] = min(abs(ts - t0));
                        xmin = DLS(minInd,6)/1000;
                        ymin = DLS(minInd,7)/1000;
                        zmin = DLS(minInd,8)/1000;
                        rmin = sqrt(xmin^2+ymin^2);
                        
                    end
                    
                    triggered = find_triggered(settings,g);
                    
                    fprintf(fID,'%4.4i\t%s\t%6.7f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%s\n',...
                        count,hhmmss,t0,xmin,ymin,zmin,rmin,triggered);
                    
                    
                    
                    
                    % LDAR plot
                    ldarColorTime1(ldar_fn,'','','','','',g.t1,g.t2,...
                        str2double(settings.ldar_r),0,0,0,1,[0,0,1,0,0],settings,g);
                    xy_fig_name = sprintf('%s%s/%4.4i-1-%s--%2.2i%2.2i%2.2i-%3.3i-XY.png',...
                        bf,date,count,date,hh,mm,ss,floor(ss_dec*1000));
                    saveas(gcf,xy_fig_name)
                    close(gcf)
                    
                    % Overall Plot
                    plot_all4(g,0)
                    ovrl_fig_name = sprintf('%s%s/%4.4i-2-%s--%2.2i%2.2i%2.2i-%3.3i-Overall.png',...
                        bf,date,count,date,hh,mm,ss,floor(ss_dec*1000));
                    saveas(gcf,ovrl_fig_name)
                    close(gcf)
                    
                    g.t1 = t(peakLoc)-0.0001;
                    g.t2 = t(peakLoc)+0.0001;
                    plot_all4(g,0)
                    zmd_fig_name = sprintf('%s%s/%4.4i-3-%s--%2.2i%2.2i%2.2i-%3.3i-Zoomed.png',...
                        bf,date,count,date,hh,mm,ss,floor(ss_dec*1000));
                    saveas(gcf,zmd_fig_name)
                    close(gcf)
                    
                end
            end
        end
    end
    
    
    t1_5 = t1_5 + 300;
    
end

fclose fID;

function trigs = get_triggers(hfn,tshift,t1,t2)

% Load corresponding header file times
fId = fopen(hfn, 'r');
trigs  = fread(fId, inf, 'double') ;
fclose( fId );

% Introduce time shift before filtering time range
trigs = trigs + tshift;
% Finding the triggers between the given time range
lol=length(trigs)- nnz(trigs>t1)+1;       % Index of the lower matrix element
ul=nnz(trigs<t2);                      % Index of the upper matrix element
% Let's find two more triggers from both ends

if lol>1 ;    lol=lol-1; end

if ul < length(trigs) ;   ul=ul+1; end

% Triggers between given time range and remove time shift
trigs=trigs(lol:ul)-tshift;


function [ch_fn h_fn] = generate_ch_fn(settings,g)
% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

% checking witch graphs are on
ch_on=g.chgraphs;
ch_fn = cell(1,60);
h_fn = cell(1,10);

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
        
        
        % If file is not exist don't store the file name
        if exist(filename,'file')==0 || exist(hfilename,'file')==0
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

function ldar_fn = generate_ldar_fn(sen_set,g)

       
        if g.mm < 30
            ext=0;
        else
            ext=30;
        end
        
        
        dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
            -datenum(str2double(g.YYYY),0,0));
        
        ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
            sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);
        
    function triggered = find_triggered(settings,g)        
        
        % turn on all ch1 data
        g.chgraphs = zeros(1,60);
        g.chgraphs = [1 0 0 1 0 0 1 0 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 ...
                      1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
        [ch_fn h_fn] = generate_ch_fn(settings,g);
        
        tshift=settings.t_shift;
        
        triggered = 'none';
        wbh2 = waitbar(0,'Checking data... ');
        % Check whether sensors have triggered data
        for j=1:60;
            waitbar(j/60,wbh2,sprintf('Checking data... (%0.1f%%)',j/6*10))
            if ~isempty(ch_fn{j})
                tch = FA_Extract1(ch_fn{j},h_fn{j},g.t1,g.t2,...
                    tshift(j),settings,j);
                if ~isempty(tch)
                    if strcmp(triggered,'none')
                        triggered = sprintf('%i',ceil(j/3));
                    else
                        triggered = sprintf('%s,%i',triggered,ceil(j/3));
                    end
                end
                
            end
        end
        delete(wbh2)
