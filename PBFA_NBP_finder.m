function PBFA_NBP_finder
% % This function is writtent to find the return stroke positions using 1s
% % fast antenna data.
% 

try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

fn = 'C:\Users\sumedhe\Desktop\PBFA_NBP_TOAM5.txt';
NBP_fn = 'C:\data\Share\Nadeeka\NBP\20110814\20110814-NBP_info.txt';

% Write the header if not exit
if exist(fn,'file')
    fID = fopen(fn,'a+');
else
    fID = fopen(fn,'a+');
    fprintf(fID,'count\ttpbfa\txpbfa\typbfa\tzpbfa\tNpbfa\tIpbfa\tsen_used\tki-sqrd\ttcg\txcg\tycg\tIp\tNcg\n');
end

fprintf(fID,'\n');

fID2 = fopen([fn(1:end -4) '_peaks.txt'],'a+');

settings = h.sen_set;
g = h.g;

x0 = 0;
y0 = 0;
z0 = 0;

% Load NBP data
fID3 = fopen(NBP_fn,'r');
NBPdata = textscan(fID3,'%f %s %f %f %f %f %f %s','HeaderLines',1);
fclose(fID3);


tshift=settings.t_shift;

factor = g.factor;

% ch_legend=cell(1,60);
% 
% for i=1:20;   
%     ch_legend{i*3-2}=[settings.sen_IDs{i} ':ch1'];
%     ch_legend{i*3-1}=[settings.sen_IDs{i} ':ch2'];
%     ch_legend{i*3}=[settings.sen_IDs{i} ':ch3'];    
% end

% Total number of CGLSS sources
tot_CG = length(NBPdata{1});


try
    wbh = waitbar(0,'Loading Fast Antenna data','Name','Collecting Peak Currents');
catch
    return
end

count = 0;

cglss = [NBPdata{3},NBPdata{4}, NBPdata{5}, NBPdata{6}];


for k = 1:tot_CG
  
    tcg = cglss(k,1);
    xcg = cglss(k,2);
    ycg = cglss(k,3);
    zcg = cglss(k,4);
    %Ip = cglss(k,11);
    %Ncg = cglss(k,17); % number of sensors used by CGLSS
  
    
    
    
    % Generate ch file names according to the current CGLSS poing
    hhmmss = sec2hhmmss(tcg);
    g.hh = str2double(hhmmss(1:2));
    g.mm = floor(str2double(hhmmss(4:5))/5)*5;
    
    % Turn on all ch2 and ch3 graphs (avoid WSB and BCC)
    g.chgraphs = [0 1 1 0 1 1 0 1 1 0 0 0 0 0 0 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 ...
                  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    
    [ch_fn h_fn] = generate_ch_fn(settings,g);
    
    % time shifts for each sensor
    t_shifts =sqrt((xcg - settings.x).^2 + (ycg - settings.y).^2 + (zcg - settings.z).^2)/3e8;
    
       
    % Peak times and sensor numbers to find PBFA
  
    ts = NaN(1,10);
    sns = ts;
    sns2 = ts;
    Eps = ts;
    indx1=1;
    indx2=1;
    
    
   figure
    hold all
    
    for i=3:3:60 % only get ch3 graphs
         
                       
        try
            val = (k-1)/tot_CG + 1/tot_CG*i/60;
            waitbar(val,wbh,sprintf('Loading Fast Antenna data. (%.2f%%)',val*100))
        catch
            return
        end
        
        
        if strcmp(ch_fn{i},'')==0
            
            sn = ceil(i/3);
            
            %g.t1 = tcg + t_shifts(sn) - 1e-3;
            %g.t2 = tcg + t_shifts(sn) + 1e-3;
            %tshift(i) = tshift(i) - t_shifts(sn);
            
            g.t1 = tcg - 10e-3;
            g.t2 = tcg + 10e-3;
            
            
            [t y ch_freq_str] = FA_Extract1(ch_fn{i},h_fn{i},g.t1,g.t2,tshift(i)- t_shifts(sn),settings,i);
            %plot(t,y)
            
            % distance            
            r = sqrt((settings.x(sn)-xcg)^2 + (settings.y(sn)-ycg)^2+ settings.z(sn)^2)/1000;
            
            if ~isempty(t)
                
                % Replace NaNs with zero so that filtfilt will work fine
                yn = y;
                yn(isnan(yn))=0;
                
                % filter out low frequencies
                Fs = 1/(t(2)-t(1));
                
                if isnan(Fs)
                    Fs = 1/(t(5) - t(4));
                end
                
                
                [z,p] = butter(5,5000/(Fs/2),'high'); % Create a High-pass butterworth filter;
                
                yn = filtfilt(z,p,yn);    % filter the data.
                
                % find maximum y time
                [miny tind] = max(yn);

                % if is a saturated sensor then try to load ch2 data
                if tind == 1
                    
                    [t y ch_freq_str] = FA_Extract1(ch_fn{i-1},h_fn{i-1},g.t1,g.t2,tshift(i)-t_shifts(sn),settings,i);
                    
                    if ~isempty(t)
                        
                        % Replace NaNs with zero so that filtfilt will work fine
                        yn = y;
                        yn(isnan(yn))=0;
                        
                        % filter out low frequencies
                        Fs = 1/(t(2)-t(1));
                        
                        if isnan(Fs)
                            Fs = 1/(t(5) - t(4));
                        end
                        
                        
                        [z,p] = butter(5,5000/(Fs/2),'high'); % Create a High-pass butterworth filter;
                        
                        yn = filtfilt(z,p,yn);    % filter the data.
                        
                        
                        plot(t,yn)
                        
                        % find maximum  time
                        [miny tind] = max(yn);
                    end
                    
                    if tind~=1
                        %ts = [ts t(tind)+ t_shifts(sn)];
                        ts(indx1) = t(tind)+ t_shifts(sn);
                        %dts = [dts dt];
                        %sns = [sns sn];
                        sns(indx1) = sn;
                        indx1 = indx1 + 1;
                        %plot(t,yn)
                        %plot(t(tind),yn(tind),'ro','markerfacecolor','r')
                    end
                    
                    
                else
                    
                    plot(t,yn)
                    %plot(t(tind),yn(tind),'ro','markerfacecolor','r')
                    
                    %ts = [ts t(tind)+t_shifts(sn)];
                    ts(indx1) = t(tind)+ t_shifts(sn);
                    %dts = [dts dt];
                    %sns = [sns sn];
                    sns(indx1) = sn;
                    indx1 = indx1 + 1;
                    %plot(t(tind),miny,'ro')
                                        
                    % Down sample to 0.5MHz to get peak e-field
                    
                    % Number of intervals
                    ni = round(Fs/0.5e6);
                    
                    try
                        ysum = y(tind-100:tind+100);
                    catch
                        ysum = y;
                    end
                    
                    
                    if ni > 1
                       
                        len = floor(length(ysum)/ni);
                        
                        temp1 = nan(ni,len);
                        temp2 = temp1;
                        
                        for indn = 1:ni
                            tempy= downsample(ysum,ni,indn-1);
                            temp2(indn,:)=tempy(1:len);
                            clear tempt tempy
                            
                        end
                        
                        temp2=mean(temp2);
                        
                        
                    end
                    
                    % Derivative of voltage (or the difference of the voltage)
                    dy = temp2(2:end)-temp2(1:end-1);
                    
                    %plot(t,y*factor(i))
                    %plot(t(1:end-1)+dt/2,dy*factor(i))
                    
                    maxdE = max(dy);
                    
                    if abs(maxdE*factor(i)) > 3 && r > 30
                        %Eps = [Eps maxdE];
                        %sns2 = [sns2 sn];
                        Eps(indx2) = maxdE;
                        sns2(indx2) = sn;
                        indx2 = indx2 + 1;
                    end
                end                 
            end              
        end
    end
    
    % Lets remove NaNs       
    if indx1 < 11
        sns(indx1:10) =[];
        ts(indx1:10) =[];
    end
    
    if indx2 < 11
        Eps(indx2:10) = [];
        sns2(indx2:10) = [];
    end
    
    arg.inner = sns;
    arg.outer = [];
    arg.t_in = ts;
    arg.t_out  = [];
    arg.method = 5;
    
    % Total numgber of sensors used
    Npbfa = length(sns)
    
    % request more than 4 sensors to calculate pbfa points
    if Npbfa > 4
        arg
        try
        [xpbfa ypbfa zpbfa tpbfa] = pbfa_finder(arg)
        end
        %fprintf('%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.6f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
        %        count,tpbfa,xpbfa,ypbfa,zpbfa,0,Npbfa,tcg,xcg,ycg,Ip,Ncg)
        
%         if ~isreal(zpbfa)
%             zpbfa = 0;
%         end
%         
%         ki_sqrd = cal_ki_sqrd(tpbfa,xpbfa,ypbfa,zpbfa,arg,settings)
      
%        
%         if ki_sqrd < 5 && zpbfa < 3000
%             
%             Ipbfa = findPeakCurr(settings,xpbfa,ypbfa,zpbfa,Eps,sns2);
%             sen_used = sprintf('%i,',sns);
%             count = count + 1;
%             fprintf('%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.1f\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
%                 count,tpbfa,xpbfa,ypbfa,zpbfa,Npbfa+arg.method*10,Ipbfa,sen_used(1:end-1),ki_sqrd,tcg,xcg,ycg,Ip,Ncg)
%             
%             fprintf(fID,'%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.1f\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
%                 count,tpbfa,xpbfa,ypbfa,zpbfa,Npbfa+arg.method*10,Ipbfa,sen_used(1:end-1),ki_sqrd,tcg,xcg,ycg,Ip,Ncg);
%             
%             store_peaks_info(count,[ts', sns'],fID2)
%         
%         elseif Npbfa > 5 && (ki_sqrd > 5 || zpbfa > 3000);
%             removes = remove_check(arg,settings);
%             
%             removes  = [removes.s_ki, removes.s_zs, sns'];
%             
%             removes( removes(:,1) > 5,:) = [];
%             removes( removes(:,2) > 3000,:) = [];
%             
%             if ~isempty(removes)
% 
%                 [mm indx] = min(removes(:,1));                
%                 indx = find(sns == removes(indx,3));
%                 sns(indx) = [];
%                 ts(indx) = [];
%       
%                 arg.t_in = ts;
%                 arg.inner = sns;
%                 
%                 [xpbfa ypbfa zpbfa tpbfa] = pbfa_finder(arg);
%                 ki_sqrd = cal_ki_sqrd(tpbfa,xpbfa,ypbfa,zpbfa,arg,settings);
%                 Ipbfa = findPeakCurr(settings,xpbfa,ypbfa,zpbfa,Eps,sns2);
%                 sen_used = sprintf('%i,',sns);
%                 count = count + 1;
%                 Npbfa = length(sns);
%                 
%                 fprintf('%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.1f\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
%                     count,tpbfa,xpbfa,ypbfa,zpbfa,Npbfa+arg.method*10,Ipbfa,sen_used(1:end-1),ki_sqrd,tcg,xcg,ycg,Ip,Ncg)
%                 
%                 fprintf(fID,'%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.1f\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
%                     count,tpbfa,xpbfa,ypbfa,zpbfa,Npbfa+arg.method*10,Ipbfa,sen_used(1:end-1),ki_sqrd,tcg,xcg,ycg,Ip,Ncg);
%             
%                 store_peaks_info(count,[ts', sns'],fID2)
%                 
%             end
%             
%            
%         end
            
        
%         count = count + 1;
%         
%         fprintf('%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.1f\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
%             count,tpbfa,xpbfa,ypbfa,zpbfa,Npbfa+arg.method*10,Ipbfa,sen_used(1:end-1),ki_sqrd,tcg,xcg,ycg,Ip,Ncg)
        
%         % remove possible outliars
%         r = sqrt(xcg^2 + ycg^2);
%         if ((r < 75000 && (~isreal(zpbfa) || zpbfa < 5000)) || ...
%                          (~isreal(zpbfa) || zpbfa < 8000))  ...
%              && isrange(xpbfa,xcg) && isrange(ypbfa,ycg)
%             
%             count = count + 1;
%             
%             % Get sensors used
%             sen_used = sprintf('%i',sns(1));
%             
%             for i=2:length(sns)
%                 sen_used = sprintf('%s,%i',sen_used,sns(i));
%             end
%             
%                                    
%             Ipbfa = findPeakCurr(settings,xpbfa,ypbfa,zpbfa,Eps,sns2);          
%                    
%             
%             fprintf('%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%i\t%s\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
%                 count,tpbfa,xpbfa,ypbfa,zpbfa,Ipbfa,Npbfa,sen_used,tcg,xcg,ycg,Ip,Ncg)
%             
%             fprintf(fID,'%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%i\t%s\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
%                 count,tpbfa,xpbfa,ypbfa,zpbfa,Ipbfa,Npbfa,sen_used,tcg,xcg,ycg,Ip,Ncg);
% 
% 
%         end
    end
end

delete(wbh);
fclose(fID);
fclose(fID2);
%toc


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


function yes = isrange(x1,x2)
% this function will test x1 is within some percentage of x2
% percentage 

if x1 > 100000
    percent = 20;
elseif x1 > 40000
    percent = 30;
elseif x1 > 20000
    percent = 70;
else
    percent = 100;
end
    

delta = abs(x2)*percent/100;

R1 = x2 - delta;
R2 = x2 + delta;

if x1 > R1 && x1 < R2
    yes = true;
else
    yes = false;
end


function Ip = findPeakCurr(sen_set,x,y,z,maxdV,sns)

%slope
m = [ 4.8930 3.9805 4.1938  NaN     2.1098 ...
      4.2255 NaN    2.0474  1.6661  2.0610  NaN];
%intercept
c = [0.0382 0.0638  0.2847  NaN     1.2001 ...
     0.5778 NaN     1.8081  1.0528  1.1174 NaN];

 
r = sqrt((sen_set.x - x).^2 + ...
         (sen_set.y - y).^2 + ...
         (sen_set.z - z).^2)/1000;

Ip = nanmean( m(sns).* maxdV .* (r(sns).^1.13) + c(sns));


function ki_sqrd = cal_ki_sqrd(t,x,y,z,arg,sen_set)
        
        % Ki -sqrd
        tt1 = arg.t_in';
        
        deg_free = length(arg.t_in);
        
        tt2 = t+sqrt((x - sen_set.x(arg.inner)).^2 + ...
                     (y - sen_set.y(arg.inner)).^2 + ...
                     (z - sen_set.z(arg.inner)).^2)'/3e8;
                          
        ki_sqrd = 1/deg_free*sum(((tt1 - tt2)/0.2e-6).^2);

function store_peaks_info(index,peaks,fID2)
        
    peaks = sortrows(peaks,2);
    
    fprintf(fID2,'\n%i',index);    

    j = 1;
    
    for i = 1:20
        try
            if i == peaks(j,2)
                fprintf(fID2,'\t%0.8f',peaks(j,1));
                j = j + 1;
            else
                fprintf(fID2,'\tNaN');
            end
        catch
            fprintf(fID2,'\tNaN');
        end
    end

function removes = remove_check(arg,sen_set);
    
    L = length(arg.t_in);
    
    if L > 5
        removes.s_ki = NaN(L,1);
        removes.s_zs = NaN(L,1);
        
        arg.method = 5;
        ts = arg.t_in;
        sns = arg.inner;
        
        for i=1:L
            
            arg.t_in = ts;
            arg.inner = sns;
            
            arg.t_in(i) = [];
            arg.inner(i) = [];
            
            [x,y,z,t]=pbfa_finder(arg);
            removes.s_ki(i,1) = cal_ki_sqrd(t,x,y,z,arg,sen_set);
            removes.s_zs(i,1) = z;
        end
    end
    

