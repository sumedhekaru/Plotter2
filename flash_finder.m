function flash_finder
% This function will use LDAR2 data to find IC and CG flashes and try to
% catogorize (will use Lu et al. 2012 for CG catogorization.


%% User Inputs
ts.t1 = 74920.329;          % start time to load data
ts.t2 = 86400;          % End time to load data
ts.dT = 0.1;            % Maximum time between flashes
ts.dR = 30000;          % Maximum spacial seperation between flashes
ts.fCGH = 5000;         % Average LDAR2 max height for mailed CGs
ts.index = 8319;        % start from zero
ts.bf = 'C:\Users\sumedhe\Desktop\Flash_categorizer\20110814\';
ts.fn = 'flash_finder_test.txt'; % File name to save data




%% Start working

handles=guidata(findall(0,'Tag','plotter2'));

if exist([ts.bf ts.fn],'file')
    ts.fID = fopen([ts.bf ts.fn],'a+');
else
    ts.fID = fopen([ts.bf ts.fn],'a+');
    % write the header
    fprintf(ts.fID,'Index\tFileName\tt1(hhmmss)\tt1(s)\tt2(s)\tN-DLS\tN-CGLSS\tAvg-X\tAvg-Y\tSize\tTriggred\tFlashType\n');
end
   


% Load LDAR2 points
wbh = waitbar(0.1,'Loading LDAR2 data...','name','Long Range LDAR2 plot');



g = handles.g;
settings = handles.sen_set;

t = floor(ts.t1/1800)*1800;

% Time shift for LDAR and Linet
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
DLSt = [];  DLSx = [];  DLSy = [];  DLSz = [];
CGt  = [];  CGx  = [];  CGy  = [];  CGz   = [];


while t < ts.t2
   
    % Generating file name for LDAR data
    [hms hh mm] = sec2hhmmss(t);
    if mm < 30
        ext=0;
    else
        ext=30;
    end
    
    dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
        -datenum(str2double(g.YYYY),0,0));
    
    ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,hh,ext);
    
    % Exit if the file does not exsist
    if ~exist(ldar_fn,'file')
        errordlg(ldar_fn,'File not found')
    end
    
    %fprintf('%s\n',ldar_fn)
    
    % Load LDAR2 data
    
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,ts.t1,ts.t2,str2double(settings.ldar_r),...
        x0,y0,z0,0);
   
   DLSt = [DLSt; DLS(:,10)];
   DLSx = [DLSx; DLS(:,6)];
   DLSy = [DLSy; DLS(:,7)];
   DLSz = [DLSz; DLS(:,8)];
   
   CGt  = [CGt; CG(:,10)];
   CGx  = [CGx ; CG(:,6)];
   CGy  = [CGy ; CG(:,7)];
   CGz  = [CGz ; CG(:,8)];
       
    t = t+1800;
end


% nothing to plot? exit
if sum(isnan(DLSt)) == length(DLSt) && ...
        sum(isnan(CGt)) == length(CGt)
    disp('No data in given time range')
    delete(wbh)
    return
end

data(:,1) = [DLSt; CGt];
data(:,2) = [DLSx; CGx];
data(:,3) = [DLSy; CGy];
data(:,4) = [DLSz; CGz];
data(:,5) = [zeros(size(DLSt)); zeros(size(CGt))+1];

data = sortrows(data,1);


L = length(data(:,1));


% Let's count the first point
if data(1,5) == 1
    type = 1;
    nCGLSS = 1;
    nDLSS = 0;
else
    type = 0;
    nDLSS = 1;
    nCGLSS = 0;
end

up = 0;
down = 0;

ts.nCG = 0;
ts.nIC = 0;



for i = 2:L;
    try
        waitbar(0.1+0.9*i/L,wbh,'Analyzing flashes')
    catch
        disp('User stopped the programm')
        fclose(ts.fID);
        return
    end
          
    
    dt = data(i,1) - data(i-1,1);
    
    if dt > ts.dT;
        % what is the previos flash type?
        
        
        max_dR = maxDistance(data(i-nDLSS-nCGLSS:i-1,:));
        tt1 = floor(data(i-nDLSS-nCGLSS,1)*1000)/1000;
        tt2 = ceil(data(i-1,1)*1000)/1000;
        tt1hms = sec2hhmmss(tt1);
        
        ts.tt1 = tt1;
        ts.tt2 = tt2;
        
        x = nanmean(data(i-nDLSS-nCGLSS:i-1,2))/1000;
        y = nanmean(data(i-nDLSS-nCGLSS:i-1,3))/1000;        
        
        if max_dR > ts.dR
            comment = 'Multi Flash';
            
            ts.index = ts.index + 1;
            
            %file name to save
            ms = round((ts.tt1 - floor(ts.tt1))*1000);            
            ts.fig_fn = sprintf('%5.5i-%s%s%s_%5.0f_%3.3i-%s',...
                ts.index,g.YYYY{:},g.MM{:},g.DD{:},ts.tt1,ms,comment(1:2));
            
           
            ts = check_triggers(handles,ts);           
            
                       
            fprintf(ts.fID,'%5.5i\t%s\t%s\t%.3f\t%0.3f\t%i\t%i\t%0.1f\t%0.1f\t%0.1f\t%s\t%s\n',...
                ts.index,ts.fig_fn,tt1hms,tt1,tt2,nDLSS,nCGLSS,x,y,max_dR/1000,ts.trigrd,comment);
            
          
            ts = resolve_multi_flashes(data(i-nDLSS-nCGLSS:i-1,:),ts,handles);
            
        elseif type
               if up > down
                   comment = 'CG : Hybrid';             
               else
                   comment = 'CG : -';
               end
               ts.nCG = ts.nCG + 1;
        else
            if up == down
                comment = 'IC';
                ts.nIC = ts.nIC + 1;
            elseif up > down
                  comment = 'IC : +';   
                  ts.nIC = ts.nIC + 1;
            else
                if mean(data(i-nDLSS-nCGLSS:i-1,4)) < ts.fCGH
                    comment = 'CG : failed';
                    ts.nCG = ts.nCG + 1;
                else
                    comment = 'IC : -';
                    ts.nIC = ts.nIC + 1;
                end
   
            end
            
        end
      
                
        if ~strcmp(comment,'Multi Flash')
            
            ts.index = ts.index + 1;
            
            %file name to save
            ms = round((ts.tt1 - floor(ts.tt1))*1000);            
            ts.fig_fn = sprintf('%5.5i-%s%s%s_%5.0f_%3.3i-%s',...
                ts.index,g.YYYY{:},g.MM{:},g.DD{:},ts.tt1,ms,comment(1:2));
            ts = check_triggers(handles,ts);
                       
            fprintf(ts.fID,'%5.5i\t%s\t%s\t%.3f\t%0.3f\t%i\t%i\t%0.1f\t%0.1f\t%0.1f\t%s\t%s\n',...
                ts.index,ts.fig_fn,tt1hms,tt1,tt2,nDLSS,nCGLSS,x,y,max_dR/1000,ts.trigrd,comment);
           
        end

       % Let's count the first point of the next flash
        if data(i,5) == 1
           type = 1;
           nCGLSS = 1;
           nDLSS = 0;
       else
           type = 0;
           nDLSS = 1;
           nCGLSS = 0;
        end
       
        up = 0;
        down = 0;
      
       
    else

        % Now check whether it is CG by checking CGLSS points
        if data(i,5) == 1
            type = 1;
            nCGLSS = nCGLSS + 1;
        else
            nDLSS = nDLSS + 1;
        end
        
        % Going up or down?
        if data(i,4) - data(i-1,4) > 0
            up = up + 1;
        else
            down = down + 1;
        end
        
    end
    
        
end

delete(wbh)
fclose(ts.fID);
ts

function maxR = maxDistance(p)

[N ~] = size(p);

% generate N*(N-1)/2 x 2 array of pair indices
[x y]=meshgrid(1:N);
i=find(triu(ones(N)-eye(N)));
i=[x(i) y(i)];

if ~isempty(i)
    
    % compute all pairwise distances
    d2=(p(i(:,1),2)-p(i(:,2),2)).^2 + ...
        (p(i(:,1),3)-p(i(:,2),3)).^2;
    
    % find minimum distance and corresponding points
    [d2min ~] = max(d2);
    maxR=sqrt(d2min);
else
    maxR = 0;
end

function ts = resolve_multi_flashes(d,ts,handles)

[L ~] = size(d);

A(1,:) = d(1,:);
B = [];
C = [];
D = [];
E = [];
F = [];



for i = 2:L
    dA = sqrt((d(i,2)-A(end,2))^2+(d(i,3)-A(end,3))^2);
    
    if dA < ts.dR
        A(end + 1,:) = d(i,:);
    else
        try dB = sqrt((d(i,2)-B(end,2))^2+(d(i,3)-B(end,3))^2);
        catch; dB = 0; end
        
        if dB < ts.dR
            B(end + 1,:) = d(i,:);
        else
            try dC = sqrt((d(i,2)-C(end,2))^2+(d(i,3)-C(end,3))^2);
            catch; dC = 0; end
            if dC < ts.dR
                C(end + 1,:) = d(i,:);
            else
                try dD = sqrt((d(i,2)-D(end,2))^2+(d(i,3)-D(end,3))^2);
                catch; dD = 0; end
                
                if dD < ts.dR
                    D(end + 1,:) = d(i,:);
                else
                    try dE = sqrt((d(i,2)-E(end,2))^2+(d(i,3)-E(end,3))^2);
                    catch; dE = 0; end
                    
                    if dE < ts.dR
                        E(end + 1,:) = d(i,:);
                    else                  
                        try dF = sqrt((d(i,2)-F(end,2))^2+(d(i,3)-F(end,3))^2);
                        catch; dF = 0; end
                        
                        if dF < ts.dR
                            F(end + 1,:) = d(i,:);
                        end
                    end
                end
            end
        end
    end
end
    
ts = find_type(A,ts,handles);
if ~isempty(B); ts = find_type(B,ts,handles); end
if ~isempty(C); ts = find_type(C,ts,handles); end
if ~isempty(D); ts = find_type(D,ts,handles); end
if ~isempty(E); ts = find_type(E,ts,handles); end
if ~isempty(F); ts = find_type(F,ts,handles); end


function ts = find_type(d,ts,handles)

    g = handles.g;
        
    [L ~] = size(d);

    if d(1,5) == 1
        type = 1;
        nCGLSS = 1;
        nDLSS = 0;
    else
        type = 0;
        nDLSS = 1;
        nCGLSS = 0;
    end

    up = 0;
    down = 0;
    
    x = nanmean(d(:,2));
    y = nanmean(d(:,3));
    
    for i = 2:L
        
        % Now check whether it is CG by checking CGLSS points
        if d(i,5) == 1
            type = 1;
            nCGLSS = nCGLSS + 1;
        else
            nDLSS = nDLSS + 1;
        end
        
        % Going up or down?
        if d(i,4) - d(i-1,4) > 0
            up = up + 1;
        else
            down = down + 1;
        end
        
    end
    
    if type
        if up > down
            comment = 'CG(M) : Hybrid';
        else
            comment = 'CG(M) : -';
        end
        ts.nCG = ts.nCG + 1;
    else
        if up == down
            comment = 'IC(M)';
            ts.nIC = ts.nIC + 1;
        elseif up > down
            comment = 'IC(M) : +';
            ts.nIC = ts.nIC + 1;
        else
            if mean(d(:,4)) < ts.fCGH
                comment = 'CG(M) : failed';
                ts.nCG = ts.nCG + 1;
            else
                comment = 'IC(M) : -';
                ts.nIC = ts.nIC + 1;
            end
            
        end
    end
    tt1 = floor(d(1,1)*1000)/1000;
    tt2 = ceil(d(end,1)*1000)/1000;
    tt1hms = sec2hhmmss(tt1);
    
    max_dR = maxDistance(d);
    
    ts.index = ts.index + 1;
    
    %file name to save
    ms = round((ts.tt1 - floor(ts.tt1))*1000);
    ts.fig_fn = sprintf('%5.5i-%s%s%s_%5.0f_%3.3i-%s',...
        ts.index,g.YYYY{:},g.MM{:},g.DD{:},ts.tt1,ms,comment(1:5));
    
    ts.tt1 = tt1;
    ts.tt2 = tt2;
    ts = check_triggers(handles,ts);   
    
    
    fprintf(ts.fID,'%5.5i\t%s\t%s\t%.3f\t%0.3f\t%i\t%i\t%0.1f\t%0.1f\t%0.1f\t%s\t%s\n',...
        ts.index,ts.fig_fn,tt1hms,tt1,tt2,nDLSS,nCGLSS,x,y,max_dR/1000,ts.trigrd,comment);
    

    
function ts = check_triggers(handles,ts)
sen_set = handles.sen_set;
g = handles.g;

hms = sec2hhmmss(ts.tt1);
g.hh = str2double(hms(1:2));
g.mm = floor(str2double(hms(4:5))/5)*5;
g.ss = 0;

g.chgraphs = [0 0 1 0 0 1 0 0 1 0 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 ...
              1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

          
g.t1 = ts.tt1;
g.t2 = ts.tt2;        

output = plot_all4(g,0);
saveas(gcf,[ts.bf ts.fig_fn '.png'])
close(gcf)


ts.trigrd = output.trigrd;

% Plot LDAR2
if g.mm < 30
    ext=0;
else
    ext=30;
end


dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
    sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);

ldarColorTime1(ldar_fn,'','','',g.t1,g.t2,...
    str2double(sen_set.ldar_r),0,0,0,1,[0,0,1,0,0],sen_set);
saveas(gcf,[ts.bf ts.fig_fn '-XY' '.png'])
close(gcf)





    

