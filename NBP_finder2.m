function NBP_finder2
% This function is inteneded to find NBP using fast antenna data.
% History 
%       2014-04-16 -- Modified to find negative NBPs from NBP_finder1
%

%% User inputs
t1 = 0;
t2 = 86400;
sn = 9;
ch = 1;
pre = 0.35;
post = 0.15;
v_thresh = -0.1; % Threshold volatage for peak findings.
% Base folder to save stuffs
%bf = 'C:\Users\Sumedhe\Desktop\NBP-test\';
bf = 'C:\Users\Sumedhe\Desktop\NBP\20110801\RawData\';

%% Setup
count = 1;
% get Plotter2 data
try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

g=h.g;
settings = h.sen_set;
date = sprintf('%s%s%s',g.YYYY{:},g.MM{:},g.DD{:});


% save data in txt file
sv_fn = sprintf('%s/%s/%s/%s-NBP_info.txt',bf,date,settings.sen_IDs{sn},date)



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
figure(100)
hold all
tools2fig
box on

while t1_5 <= t2_5 - 300;

    g.t1 = t1_5;
    g.t2 = t1_5 + 300;
    [hhmmss, hh, mm] = sec2hhmmss(g.t1);
    g.hh = hh;
    g.mm= floor(mm/5)*5;
    g.ss = 0;
    
    % generate file names
    [ch_fn, h_fn] = generate_ch_fn(settings,g);
    
    % Get CGLSS data for this 5 mins period
    [cgt, cgx, cgy] = get_cglss(settings,g,sn); 
    
    % Get triggers
    if ~isempty(h_fn{ind0})
        trigs = get_triggers(h_fn{ind0},tshift,g.t1,g.t2);
        L = length(trigs);
        
        % check each trigger
        for i = 9:L
            try
            t1 = trigs(i);
            t2 = trigs(i) + pre + post;
            fprintf('%13.6fs\n',t1)
            
            [t,y]=FA_Extract1(ch_fn{ind0}, h_fn{ind0} ,t1,...
                t2,tshift,settings,sn);
            
            if ~isempty(t)
                
                % Apply a filter to get rid of DC
                [t,yf] = ch_high_pass(t,y,1000);
                %yf = y;
                
                % Find negative peaks
                peakLoc = peakfinder(yf,-v_thresh,v_thresh,-1);
                
                L1 = length(peakLoc);
                
                % starting point is a peak beacuse of the filter 
                % let's ignore it
                
                if L1 > 0 && peakLoc(1) < 20
                    startInd = 2;
                else
                    startInd = 1;
                end
                
                for j = startInd:L1
                                       
                    
                    % NBP test 1 (zero cross time)                    
                    if peakLoc(j) > 20
                        ntest = 20;
                    else
                        ntest = peakLoc(j) - 1;
                    end
                    
                    % get the begining time
                    for k = 1:ntest
                        if yf(peakLoc(j)-k) > 0
                            % found it
                            break
                        end
                    end
                    
                    beg_t = t(peakLoc(j)-k);
                    beg_v = yf(peakLoc(j)-k);
                    
                    if peakLoc(j) > length(t)-20
                        ntest = 20;
                    else
                        ntest = length(t) - peakLoc(j)-1;
                    end
                    
                    % get the xero cross time
                    for k = 1:ntest
                        if yf(peakLoc(j)+k) > 0
                            % found it
                            break
                        end
                    end
                    
                    zc_t = t(peakLoc(j)+k);
                    zc_v = yf(peakLoc(j)+k);
                    
                    % zero cross time
                    zeroCrTime = zc_t - beg_t;
                    
                    if zeroCrTime < 15e-6
                        
                        % Possible -NBP but let's do the second test
                        fprintf('Possible -NBP found.')
                        
                        % Test 2: Sorrounding peaks
                        % If this is a NBP, there shouldn't be any other
                        % sorrounding peaks
                        
                        lol = sum(t < (t(peakLoc(j)) - 200e-6))+1;
                        ul  = sum(t < (t(peakLoc(j)) + 200e-6));
                        
                        if length(peakfinder(yf(lol:ul),-v_thresh,v_thresh,-1)) == 1
                            % Ok this is the minimum, means this could be a
                            % NBP
                            fprintf(' 2nd test passed.')
                            
                            
                            % Test 3: Is it a RS?
                            if min(abs(cgt - t(peakLoc(j)))) > 25e-6
                                % this is potentially not a CGLSS
                                fprintf(' Confirmed!\n')                                
                                
                                cla
                                plot(t(lol:ul), yf(lol:ul))
                                temp = peakfinder(yf(lol:ul),-v_thresh,v_thresh,-1);
                                plot(t(temp+lol-1),yf(temp+lol-1),'ko')
                                plot(t(peakLoc(j)),yf(peakLoc(j)),'ro')
                                plot(beg_t,beg_v,'go')
                                plot(zc_t,zc_v,'go')
                                xlim([t(lol) t(ul)])
                                
                                count = count + 1;
                                sv_fn = sprintf('%s/%s/%s/%4.4i.png',...
                                    bf,date,settings.sen_IDs{sn},count);
                                saveas(gcf,sv_fn)
                                
                            else
                                cla
                                plot(t(lol:ul), yf(lol:ul))
                                plot(t(peakLoc(j)),yf(peakLoc(j)),'ro')
                                plot(beg_t,beg_v,'go')
                                plot(zc_t,zc_v,'go')
                                xlim([t(lol) t(ul)])
                                fprintf(' Sorry. It was a RS.\n')
                            end
                        else
                                                        
                            cla
                            plot(t(lol:ul), yf(lol:ul))
                            temp = peakfinder(yf(lol:ul),-v_thresh,v_thresh,-1);
                            plot(t(temp+lol-1),yf(temp+lol-1),'ko')                          
                            plot(t(peakLoc(j)),yf(peakLoc(j)),'ro')
                            plot(beg_t,beg_v,'go')
                            plot(zc_t,zc_v,'go')
                            xlim([t(lol) t(ul)])
                            %title(length(temp))
                            
                            
                                
                            fprintf(' Sorry. 2nd test failed\n')
                        end
                        
                    end
                    
                end
                
            end
            catch
                % this trigger has problem just goes to next
                
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
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:})

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
        
function [cgt, cgx, cgy] = get_cglss(settings,g,sn)
    
cg_fn = sprintf('%s/cglss/%s/%s/KSCCGLSS%s%s%s.dat',...
            settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},g.MM{:},g.DD{:});
x0 = settings.x(sn);
y0 = settings.y(sn);
z0 = settings.z(sn);


data=CGLSS_extract(cg_fn,g.t1,g.t2,200000,x0,y0,z0);

cgt = data(:,13);
cgx = data(:,2);
cgy = data(:,3);
    
