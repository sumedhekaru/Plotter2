function rhessi_plot
cccc
sen_set = open('sensor_setting.mat');

% Open RHESSI satelite location data from Thomas
rh_fn = 'C:\Users\Sumedhe\Desktop\Rhessi_analysis\RHESSI_near_KSC_2011.txt';

fID = fopen(rh_fn);
data = textscan(fID,'%s %s %f %f %f');
fclose(fID);

%fg = figure;
%set(fg,'position',[200   100   1200   1000])
%hold all
%xlabel('East (km)')
%ylabel('West (km)')

% Manual Starting point
start_date = '20110716'; % yyymmmdd format -- keep empty to start from the begining
nTot = 4266; % Default 0


dates = data{1};
times = data{2};


if ~isempty(start_date)
    ind0 = sum(datenum(dates,'dd-mmm-yyyy') < datenum(start_date,'yyyymmdd'))+1;
    
    dates = dates(ind0:end);
    times = times(ind0:end);
  
else
    ind0 = 1;   
end



for i = 1:5
    d_temp{i} = data{i}(ind0:end,:);
end


data = d_temp;
data{1} = data{3};
data{2} = data{3};
data = cell2mat(data);

% First data point
t = hhmmss2sec(times{1},0)+1.8e-3; % 1.8 ms time correction for Rhessi data
[x,y] = latlon2xy(data(1,4),360-data(1,3));
st = 1;
L = length(times);

sat_d.fID = fopen('C:\Users\Sumedhe\Desktop\Rhessi_analysis\2011_E-change_data_for_rhessi2.txt','a+');

for i = 2:L
    
    tc = hhmmss2sec(times{i},0)+1.8e-3;
    [xc,yc] = latlon2xy(data(i,4),360-data(i,3));
       
       
    if abs(tc - t(i-1)) > 40
        en = i-1;
        hold all
        plot(x(st:en)/1000,y(st:en)/1000,'r-');
        plot(0,0,'ro')
        for j = st:en; text(x(j)/1000,y(j)/1000,sprintf('   %i s',t(j))); end
        
             
        sat_d.t = t(st:en);
        sat_d.x = x(st:en)/1000;
        sat_d.y = y(st:en)/1000;
        sat_d.lat = data(st:en,4)';
        sat_d.lon = 360 - data(st:en,3)';
        
        [nLDAR nTot] = plot_ldar2(nTot+1,dates{st},t(st:en),x(st:en)/1000,y(st:en)/1000,sen_set,sat_d);
        
        %fn = sprintf('%s %6.6i--%6.6is nLDAR = %5.5i',datestr(dates{st},'yyyymmdd'),t(st),t(en),nLDAR);
        %title(sprintf('%s     %i -- %i s   nLDAR = %i',datestr(dates{st},'yyyy-mm-dd'),t(st),t(en),nLDAR))
        
        %florida_map
        
        %daspect([1 1 1]);
        %box on
        
        %saveas(gcf,['C:\Users\Sumedhe\Desktop\Rhessi_analysis\' fn],'fig')
        %saveas(gcf,['C:\Users\Sumedhe\Desktop\Rhessi_analysis\' fn],'jpg')
        

                   
        st = i;
       
    end
    
    t = [t tc];
    x = [x xc];
    y = [y yc];
   
    clf
   
end

fclose(sat_d.fID);

% data = [ ...
%     86060       276.40125       23.899845 
%     86080       277.50611       24.532662 
%     86100       278.62280       25.155983
%     86120       279.75086       25.769074
%     86140       280.89224       26.372405
%     86160       282.04623       26.965116
%     86180       283.21231       27.546476
%     86200       284.39242       28.116889
%     86220       285.58576       28.675510];
% 
% 
% [x,y] = latlon2xy(data(:,3),360-data(:,2));
% 
% x = x/1000;
% y = y/1000;
% r = sqrt(x.^2 + y.^2);
% 
% 
% 
% figure
% hold all
% box on
% daspect([1 1 1])
% plot(x,y,'r-o')
% scatter(x,y,50,data(:,1),'filled')
% plot(0,0,'ro')
% text(x,y,{'  86060 s' '  86080 s' '  86100 s' '  86120 s' '  86140 s' '  86160 s' '  86180 s' '  86200 s' '  86220 s'})
% st = 1;
% en = length(x);
% dates = {'14-Aug-2011'};
% t = data(:,1)';
% x = x';
% y = y';
% 
% sat_d.t = t;
% sat_d.x = x;
% sat_d.y = y;
% sat_d.lat = data(:,3);
% sat_d.lon = 360 - data(:,2);
% 
% nLDAR = plot_ldar2(dates{st},t(st:en),x(st:en)/1000,y(st:en)/1000,sen_set,sat_d);
% 
% florida_map
% xlabel('East (km)')
% ylabel('West (km)')


function [nLDAR nTot] = plot_ldar2(nTot,day,t,x,y,sen_set,sat_d)


day = datestr(day,'ddmmyyyy');
dnum = datenum(day,'ddmmyyyy')- datenum(day(5:8),'yyyy')+1;
hms = sec2hhmmss(t(1));
mm = floor(str2double(hms(4:5))/30)*30;


ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%s%2.2i.txt',...
        sen_set.base_dir,day(5:8),day(3:4),day(1:2),day(5:8),dnum,hms(1:2),mm);
    
try    
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,t(1),t(end),1000000,...
            0,0,0,0,0);
catch
    DLS = nan(1,8);
    CG  = nan(1,8);
end

nTot = IC_analyize(nTot,CG,DLS,sen_set,day,sat_d);

%scatter([DLS(:,6)/1000 ; x'] ,[DLS(:,7)/1000; y'],25,[DLS(:,1); t'],'filled')

nLDAR = sum(~isnan(DLS(:,6)));


function nTot = IC_analyize(nTot,CG,DLS,sen_set,day,sat_d)


% User inputs (Flash thresholds)
ts.dT = 0.1;
ts.dR = 15000;


nCG = 0;
nIC = 0;


data(:,1) = [DLS(:,1); CG(:,1)];
data(:,2) = [DLS(:,6); CG(:,6)];
data(:,3) = [DLS(:,7); CG(:,7)];
data(:,4) = [DLS(:,8); CG(:,8)];
data(:,5) = [zeros(size(DLS(:,1))); zeros(size(CG(:,1)))+1];
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
                
                type = catogorize_flash(nTot,data,indx,sen_set,day,sat_d);
                
                if strcmp(type,'CG'); nCG = nCG + 1;
                else nIC = nIC + 1;
                end
            end                
            
        else
            indx = flashStart:flashEnd;
            nTot = nTot + 1;
            type = catogorize_flash(nTot,data,indx,sen_set,day,sat_d);
            
            if strcmp(type,'CG'); nCG = nCG + 1;
            else nIC = nIC + 1;
            end            
        end
                        
        flashStart = i+1; 
        
    end
end

%nCG
%nIC
%nTot



function type = catogorize_flash(counter, data,indx,sen_set,day,sat_d)

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

if nLDAR2 > 4 && strcmp(type,'IC')
    
    g.chgraphs = [0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 ...
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    
    hms = sec2hhmmss(data(indx(1),1));
    
    [ch_fn, h_fn] = generate_ch_fn(sen_set,g,day,hms);
    
    % Get apropriate time shifts
    tshiftd = csvread('time_shifts_2011.csv',1,0);
    tshifts = get_tshifts(day,hms,tshiftd);
    
    %figure    
    %xlim([data(indx(1),1),data(indx(end),1)])
    %hold all
    
    % find out wihch sensor has triggered data
    trigd = [];
    
    % Find out closest sensor number, amplitude of the biggest pulse,
    % closest LDAR2 point
    cl_dTxyzr = [nan nan nan nan inf];
    cl_sn = 20;
    cl_pkV = nan;
    cl_pkVn = nan;
    cl_pkT = nan;
    pl_oc_t = nan;
    
    
    trigIp = [];
    Vps    = [];
    
        
    
    for i=1:30
        
        % Sensor number
        sn = ceil(i/3);
        
        if ~isempty(ch_fn{i})

             [tt,yy]=FA_Extract1(ch_fn{i}, h_fn{i},data(indx(1),1),...
                 data(indx(end),1),tshifts(i),sen_set,i);
             
           
             if length(tt) ~= sum(isnan(tt))
                 
                 
                 % Get the correct calibration factors
                 factor = calibration_local(day);
               
                 yy = yy*factor(i);
                 
                 %[yy yH] = hill_tra(tt,yy,2000);
                 [tt,yy] = ch_high_pass(tt,yy,2000);
                                  

                 
                 % Find out the biggest positive (IC) pulse
                 [pk_v pkInd] = max(yy);              
                 pk_t = tt(pkInd);
                 
                 % Keep sensor numbers of triggered data
                 if isempty(trigd)
                     trigd = sprintf('%i',sn);
                 else
                     trigd = sprintf('%s,%i',trigd,sn);
                 end 
                 
                 % For peak current calculations
                 trigIp = [trigIp sn];
                 Vps    = [Vps pk_v/factor(i)];
                 
                 x0 = sen_set.x(sn);
                 y0 = sen_set.y(sn);
                 z0 = sen_set.z(sn);
                 
                 % find closest LDAR2 point in time for the peak
                 t_ar = data(:,1) + sqrt((data(:,2)-x0).^2 + ...
                                         (data(:,3)-y0).^2 + ...
                                         (data(:,4)-z0).^2)/3e8; % Arrival time
                 dts = t_ar - pk_t;
                 
                 % Closest LDAR2 location
                 [cl_dt cl_ind] = min(abs(dts));
                 
                 cl_dt = dts(cl_ind);
                 
                 r = sqrt((data(cl_ind,2)-x0)^2+(data(cl_ind,3)-y0)^2 + (data(cl_ind,4)-z0)^2);
                 
                 % Find out closest sensor number, amplitude of the biggest pulse,
                 % closest LDAR2 point
                 if r < cl_dTxyzr(5)
                         cl_sn = sn;
                         cl_pkV = pk_v;
                         cl_pkVn = pk_v*r/100000;
                         cl_pkT = pk_t;
                         cl_dTxyzr = [cl_dt data(cl_ind,2:4) r];
                                             
                         
                         % Pulse occured time (assuming it occured at LDAR2 location)
                         pl_oc_t = pk_t - sqrt((data(cl_ind,2)-x0).^2 + ...
                                               (data(cl_ind,3)-y0).^2 + ...
                                               (data(cl_ind,4)-z0).^2)/3e8;
                         
%                          cla
%                          box on
%                          plot(tt,yy)
%                          plot(pk_t,pk_v,'ro','markerfacecolor','r')
%                          if ~isnan(pk_t)
%                             xlim([pk_t - 100e-6, pk_t + 100e-6])
%                          end
%                          legend({sen_set.sen_IDs{sn},sprintf('r = %0.1f km',cl_dTxyzr(5)/1000)})
                 end                 
             end
        end
    end

   
     % Find out aproximate sat location when receiving the sig.
     
     % location of the Satelite when the pulse happening
     dts = sat_d.t -pl_oc_t;
     dts = dts + (dts > 0)*50000;     
     [mmm ind] = min(abs(dts));
     
         
     xs = sat_d.x(ind):(sat_d.x(ind+1)-sat_d.x(ind))/1e6:sat_d.x(ind+1);
     ys = sat_d.y(ind):(sat_d.y(ind+1)-sat_d.y(ind))/1e6:sat_d.y(ind+1);
     zs = zeros(size(xs))+593;
     ts = sat_d.t(ind):(sat_d.t(ind+1)-sat_d.t(ind))/1e6:sat_d.t(ind+1);
     lats = sat_d.lat(ind):(sat_d.lat(ind+1)-sat_d.lat(ind))/1e6:sat_d.lat(ind+1);
     lons = sat_d.lon(ind):(sat_d.lon(ind+1)-sat_d.lon(ind))/1e6:sat_d.lon(ind+1);
     
     
     % Pulse signal arrival time at each sat location
     ts2 = pl_oc_t + sqrt((cl_dTxyzr(2)-xs*1000).^2 + ...
                          (cl_dTxyzr(3)-ys*1000).^2 + ...
                          (cl_dTxyzr(4)-zs*1000).^2)/3e8;
                      
     [mmm ind] = min(abs(ts - ts2));
     
     % Location and time of the satelite when it receive the signal
     sat_x = xs(ind);
     sat_y = ys(ind);
     %sat_z = zs(ind);
     sat_t = ts(ind);
     sat_lat = lats(ind);
     sat_lon = lons(ind);
     
     % Horizontal distance from pulse to the sat
     sat_D = sqrt((cl_dTxyzr(2)/1000-sat_x).^2 + (cl_dTxyzr(3)/1000-sat_y).^2);
     
     % Filter out data when we don't trigger anything              
     if isnan(cl_pkT)
        sat_x = NaN;
        sat_y = NaN;
        sat_t = NaN;
        sat_lat = NaN;
        sat_lon = NaN;
        sat_D = NaN;
     end
               
    
%     fprintf('%6.0i%13.6f\t%13.6f\t%5.1f\t%5.1f\t%4.1f\t%5.1f\t%3.2i\t%4.2i\t%s\n',...
%     counter,data(indx(1),1),data(indx(end),1),...
%     X,Y,mean(data(indx,4))/1000,D, ...
%     nCGLSS,nLDAR2,type)

        % index date t1 t2 trigd 
        % closestSn maxPulseTime maxPulseAmplitude maxPulseNormAmpltude
        % dr X Y Z R of closest LDAR2 point
        % t X Y D of the satelite when signal arrive to the satelite
        % LAT LOn of the satelite when signal arrive to the sat
        
        Ip = calculate_Ip(Vps,cl_dTxyzr(2),cl_dTxyzr(3),cl_dTxyzr(4),trigIp,sen_set,1);
    
        
        
        fprintf('%6.0i\t%s\t%13.6f\t%13.6f\t%16s\t%s\t%4i\t%13.7f\t%5.2f\t%5.2f\t%7.1f\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%13.7f\t%6.1f\t%6.1f\t%6.1f\t%9.5f\t%9.5f\t%9.1f\n',...
            counter,datestr(datenum(day,'ddmmyyyy'),'yyyymmdd'),data(indx(1),1),data(indx(end),1),...
            trigd,sen_set.sen_IDs{cl_sn},nLDAR2,cl_pkT,cl_pkV,cl_pkVn,...
            cl_dTxyzr(1)*1e6, cl_dTxyzr(2:5)/1000,...
            sat_t,sat_x,sat_y,sat_D,...
            sat_lat,sat_lon,Ip)
        
        
        fprintf(sat_d.fID,'%6.0i\t%s\t%13.6f\t%13.6f\t%16s\t%s\t%4i\t%13.7f\t%5.2f\t%5.2f\t%7.1f\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%13.7f\t%6.1f\t%6.1f\t%6.1f\t%9.5f\t%9.5f\t%9.1f\n',...
            counter,datestr(datenum(day,'ddmmyyyy'),'yyyymmdd'),data(indx(1),1),data(indx(end),1),...
            trigd,sen_set.sen_IDs{cl_sn},nLDAR2,cl_pkT,cl_pkV,cl_pkVn,...
            cl_dTxyzr(1)*1e6, cl_dTxyzr(2:5)/1000,...
            sat_t,sat_x,sat_y,sat_D,...
            sat_lat,sat_lon,Ip);
end

% Load Ch1 data


% fprintf(ts.fID,'%6.0i%13.6f\t%13.6f\t%5.1f\t%5.1f\t%4.1f\t%5.1f\t%3.2i\t%4.2i\t%s\n',...
%     counter,data(indx(1),1),data(indx(end),1),...
%     X,Y,mean(data(indx,4))/1000,D, ...
%     nCGLSS,nLDAR2,type);


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

function [ch_fn, h_fn] = generate_ch_fn(settings,g,day,hms)
% Base folder
ddd = datenum(day,'ddmmyyyy');

if ddd <= 734369 % 2010-08-19
    bf = 'C:\data\2010-07-01 -- 2010-08-19\data';
elseif ddd <= 734700 % 2011-07-16
    bf = 'C:\data\2011-07-07 -- 2011-07-16\data';
elseif ddd <= 734719
    bf = 'C:\data\2011-07-17 -- 2011-08-04\data';
else
    bf = 'C:\data\2011-08-05 -- 2011-08-16\data';
end
    
% ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%s%2.2i.txt',...
%         sen_set.base_dir,day(5:8),day(3:4),day(1:2),day(5:8),dnum,hms(1:2),mm);

mm = floor(str2double(hms(4:5))/5)*5;

bfolder=sprintf('%s/%s/%s/%s/', ...
    bf,day(5:8),day(3:4),day(1:2));

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
        filename=sprintf('%s%s_%s%s%s_%s%2.2i%2.2i.ch%1.1i', ...
            bfolder,sid,day(5:8),day(3:4),day(1:2),hms(1:2),mm,0,ext);
        
        hfilename=sprintf('%s%s_%s%s%s_%s%2.2i%2.2i.h%1.1i', ...
            bfolder,sid,day(5:8),day(3:4),day(1:2),hms(1:2),mm,0,ext);
        
        
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


function factor = calibration_local(day)

data = csvread('calibration_data.csv',2,0);

year = day(5:8);
month = day(3:4);
date   = day(1:2);

%date value
dvalue=str2double(sprintf('%s%s%s',year,month,date));

% Lets calculate date range average anyway because daily average may not
% availbale everyday.
%lower limit

if strcmp(year,'2010')
    a.cal_start = 20100701;
    a.cal_end   = 20100830;
elseif strcmp(year,'2011')
    a.cal_start = 20110701;
    a.cal_end   = 20110830;
else
    % do nothing
end


lol =  sum(data(:,1) < a.cal_start)+1;
ul  =  sum(data(:,1) <= a.cal_end);

if lol < ul
    range_avg = nanmean(data(lol:ul,2:61));
else
    range_avg = data(lol:ul,2:61);
end

% Use daiy average value
factor=NaN(1,60);


% finding calibrations
index = find(data(:,1)==dvalue);

% found the day?
if isempty(index) == 0
    factor=data(index,2:61);
end


% Now lets replace the values of NaN with Time Range Averge
for i=1:60
    if isnan(factor(i))
        %if TRA not available lets use manual value
        if isnan(range_avg(i))
            factor(i)=data(1,i+1);
        else
            factor(i)=range_avg(i);
        end
        
    end
end

% ch1 and ch3 numbers are just gains. We need real calibration numbers.
factor(1:3:60) = factor(2:3:60)./factor(1:3:60);
factor(3:3:60) = factor(2:3:60)./factor(3:3:60);



function tshifts = get_tshifts(day,hms,thsiftd)


hh = str2double(hms(1:2));

if hh > 12
    % Current date
    date = str2double(sprintf('%s%s%s',day(5:8),day(3:4),day(1:2)));
          
else
    % Previous date
    date = sprintf('%s%s%s',day(5:8),day(3:4),day(1:2));
    daten = datenum(date,'yyyymmdd');
    date = str2double(datestr(daten-1,'yyyymmdd'));
end


index = find(thsiftd(:,1)== date);

if isempty(index)
    tshifts = zeros(1,60);
else
    % replace NaNs with zeros
    thsiftd(isnan(thsiftd)) = 0;
    
    thsiftd(index,1+1)
    
    for i=1:20
        
        tshifts(i*3-1) = thsiftd(index,i+1);
        tshifts(i*3-2) = thsiftd(index,i+1);
        tshifts(i*3)   = thsiftd(index,i+1);
    end
end
