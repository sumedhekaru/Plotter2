function flash_finder_easy
% This function will use LDAR2 data to find IC and CG flashes and try to
% catogorize (will use Lu et al. 2012 for CG catogorization.


%% User Inputs
ts.t1 = 0;      % start time to load data
ts.t2 = 86400;          % End time to load data
ts.dT = 0.1;            % Maximum time between flashes
ts.dR = 15000;          % Maximum spacial seperation between flashes
ts.fCGH = 5000;         % Average LDAR2 max height for mailed CGs
ts.index = 0;        % start from zero
ts.bf = 'C:\Users\sumedhe\Desktop\';
ts.fn = 'flash_catogories_test-20140814-75km.txt'; % File name to save data

ts.LDAR2_R = 75000; % LDAR2 radius

tic

%% Start working

handles=guidata(findall(0,'Tag','plotter2'));

if exist([ts.bf ts.fn],'file')
    ts.fID = fopen([ts.bf ts.fn],'a+');
else
    ts.fID = fopen([ts.bf ts.fn],'a+');
    % write the header
    fprintf(ts.fID,'Index\tt1\tt2\txAvg\tyAvg\tzAvg\tnCG\tnIC\tType\n\n');
end
   


% Load LDAR2 points
wbh = waitbar(0.1,'Loading LDAR2 data...','name','Long Range LDAR2 plot');



g = handles.g;
settings = handles.sen_set;

t = floor(ts.t1/1800)*1800;


DLSt = [];  DLSx = [];  DLSy = [];  DLSz = [];
CGt  = [];  CGx  = [];  CGy  = [];  CGz   = [];

nCG = 0;
nIC = 0;
nTot = 0;

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
    
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,ts.t1,ts.t2,ts.LDAR2_R,...
        0,0,0,0,0);
   
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
data(:,6) = sqrt(data(:,2).^2 + data(:,3).^2)/1000;

data = sortrows(data,1);


L = length(data(:,1));

% check the first point
if data(1,5) == 1
    isIC = 0;
    isCG = 1;
    nCG = 1;
else
    isIC = 1;
    isCG = 0;
    nIC = 1;
end

tot_flashes = 1;

dt = data(2:end,1)-data(1:end-1,1);
dt = [0; dt];

flashStart = 1;
flashEnd = 0;

for i = 2:L
    try
        msg = sprintf('nCG = %6.6i    nIC = %6.6i   nTot = %6.6i',nCG,nIC,nTot);
        waitbar(i/L,wbh,msg)
    catch
        disp('User stopped processing.')
        return
    end
    
    if dt(i) > ts.dT
        
        flashEnd = i;
       
             
        % is this a multi flash
        if  range(data(flashStart:flashEnd,2)) > ts.dR || ...
            range(data(flashStart:flashEnd,3)) > ts.dR
                                  
            [flashNo] = resolve_multiflash(data,flashStart,flashEnd,ts);
            
            maxF = max(flashNo);
            
            for k = 1:maxF
                indx = find(flashNo == k);
                
                indx = indx + flashStart - 1;
                nTot = nTot + 1;
                
                type = catogorize_flash(nTot,data,indx,ts);
                
                if strcmp(type,'CG'); nCG = nCG + 1;
                else nIC = nIC + 1;
                end
            end                
            
        else
            indx = flashStart:flashEnd;
            nTot = nTot + 1;
            type = catogorize_flash(nTot,data,indx,ts);
            
            if strcmp(type,'CG'); nCG = nCG + 1;
            else nIC = nIC + 1;
            end            
        end
                        
        flashStart = i+1; 
        
    end
end

nCG
nIC
nTot

fclose(ts.fID);

delete(wbh)
toc


function type = catogorize_flash(counter, data,indx,ts)

% Assume the flash is IC
type = 'IC';

% Is it CG?
nCGLSS = sum(data(indx,5));

if nCGLSS > 1
    type = 'CG';
end

nLDAR2 = length(indx) - nCGLSS;

X = mean(data(indx,2))/1000;
Y = mean(data(indx,3))/1000;
D = sqrt(X^2+Y^2);

fprintf('%6.0i%13.6f\t%13.6f\t%5.1f\t%5.1f\t%4.1f\t%5.1f\t%3.2i\t%4.2i\t%s\n',...
    counter,data(indx(1),1),data(indx(end),1),...
    X,Y,mean(data(indx,4))/1000,D, ...
    nCGLSS,nLDAR2,type)

fprintf(ts.fID,'%6.0i%13.6f\t%13.6f\t%5.1f\t%5.1f\t%4.1f\t%5.1f\t%3.2i\t%4.2i\t%s\n',...
    counter,data(indx(1),1),data(indx(end),1),...
    X,Y,mean(data(indx,4))/1000,D, ...
    nCGLSS,nLDAR2,type);


function flashNo = resolve_multiflash(data,flashStart,flashEnd,ts)

f = flashStart;
kindsX = data(f,2);
kindsY = data(f,3);

L = flashEnd - flashStart + 1;

flashNo = zeros(1,L);


for i = flashStart:flashEnd
    
    found = 0;
    
    for j = 1:length(kindsX)        
        if sqrt((kindsX(j) - data(i,2))^2 +(kindsY(j) - data(i,3))^2)...
                <= ts.dR
    
           flashNo(i-flashStart+1) = j;
           found = 1;
           break
        end
    end
    
    if ~found
        kindsX = [kindsX data(i,2)];
        kindsY = [kindsY data(i,3)]; 
        flashNo(i-flashStart+1) = j+1;
    end
end






       
       
       
    






