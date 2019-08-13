function PBFA_CG_finder2
% % This function is writtent to find the return stroke positions using 1s
% % fast antenna data.
% 

try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

bf = 'C:\Users\sumedhe\Desktop\PBFA_CG_auto_comp\';
fn = [bf 'PBFA_CG_data_TOAM6-new.txt'];

if ~exist(bf,'dir');    mkdir(bf); end

% Write the header if not exit
if exist(fn,'file')
    fID = fopen(fn,'a+');
else
    fID = fopen(fn,'a+');
    fprintf(fID,'count\ttpbfa\txpbfa\typbfa\tzpbfa\tNpbfa\tIpbfa\tsen_used\tki-sqrd\ttcg\txcg\tycg\tIp\tNcg\n');
end

fprintf(fID,'\n');

fID2 = fopen([fn(1:end -4) '_peaks.txt'],'a+');
fID3 = fopen([fn(1:end -4) '_real_peaks.txt'],'a+');
fID4 = fopen([fn(1:end -4) '_hilb_peaks.txt'],'a+');
fID5 = fopen([fn(1:end -4) '_dEdt_peaks.txt'],'a+');

settings = h.sen_set;
g = h.g;

cglss_fn = sprintf('%s/cglss/%s/%s/KSCCGLSS%s%s%s.dat',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},g.MM{:},g.DD{:});
x0 = 0;
y0 = 0;
z0 = 0;

% Load CGLSS data
cglss = CGLSS_extract(cglss_fn,43200,86400,2000000,x0,y0,z0);

tshift=settings.t_shift;

factor = g.factor;

% Total number of CGLSS sources
tot_CG = length(cglss(:,1));

wbh = waitbar(0,'Loading Fast Antenna data','Name','Collecting Peak Currents');

count = 1;

% Turn on all ch2 and ch3 graphs (avoid WSB and BCC)
g.chgraphs = [0 1 1 0 1 1 0 1 1 0 0 0 0 0 0 0 1 1 0 1 1 0 1 1 0 1 1 0 1 1 ...
                  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

for k = count:tot_CG
  
    tcg = cglss(k,1);
    xcg = cglss(k,2);
    ycg = cglss(k,3);
    Ip = cglss(k,11);
    Ncg = cglss(k,17); % number of sensors used by CGLSS
  
    % Generate ch file names according to the current CGLSS poing
    hhmmss = sec2hhmmss(tcg);
    g.hh = str2double(hhmmss(1:2));
    g.mm = floor(str2double(hhmmss(4:5))/5)*5;   
        
    [ch_fn, h_fn] = generate_ch_fn(settings,g);
    
    % time shifts for each sensor
    t_shifts =sqrt((xcg - settings.x).^2 + (ycg - settings.y).^2)/3e8;
           
    % Peak times and sensor numbers to find PBFA  
    ts = NaN(1,10);
    sns = ts;
    sns2 = ts;
    sns3 = ts;
    Eps = ts;
    Eps_real = ts;
    Eps_hilb = ts;
    Eps_dEdt = ts;
    indx1=1;
    indx2=1;
    indx3 = 1;
    
    g.t1 = tcg - 0.1e-3;
    g.t2 = tcg + 0.1e-3;
   
    for i=3:3:60 % only get ch3 graphs        
        
        ii = i;
        
        try
            val = (k-1)/tot_CG + 1/tot_CG*i/60;
            waitbar(val,wbh,sprintf('Loading Fast Antenna data. (%.2f%%)',val*100))
        catch
            fclose(fID); fclose(fID2); fclose(fID3); fclose(fID4); fclose(fID5);
            return
        end
        
        sn = ceil(i/3); 
        
        if strcmp(ch_fn{i},'')==0
             
            % Let's load ch3 data
            saturated = 0;
            [t, y] = FA_Extract1(ch_fn{i},h_fn{i},g.t1,g.t2,tshift(i)- t_shifts(sn),settings,i);
            %plot(t,y)
            
            % If ch3 is saturated, let's load ch2 data 
            if range(y) < 0.001
                [t, y] = FA_Extract1(ch_fn{i-1},h_fn{i-1},g.t1,g.t2,tshift(i)-t_shifts(sn),settings,i);
                ii = i -1;
                saturated = 1;
            end
            
            % distance to sensor           
            r = sqrt((settings.x(sn)-xcg)^2 + (settings.y(sn)-ycg)^2+ settings.z(sn)^2)/1000;
             
            % if there is data that are not saturated, let's continue
            if ~isempty(t) && range(y) > 0.001
                
                % filter out low frequencies and get hilbert transformation
                [yn, yH] = hill_tra(t,y,2000);
                yH = abs(yH);
                
                % find minimum y time
                [~, tind] = min(yn);              

                ts(indx1) = t(tind)+ t_shifts(sn); 
                Eps_real(indx1) = range(yn(1:tind))*factor(ii);
                
                if ~saturated
                    Eps_hilb(indx3) = max(yH);
                    sns3(indx3) = sn;
                    indx3 = indx3 + 1;
                end
                
                sns(indx1) = sn;
                
                % Down sample to 0.5MHz to get peak e-field
                % frequency
                Fs = 1/(t(2)-t(1));
                
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
                                
                maxdE = min(dy);
                Eps_dEdt(indx1) = maxdE*factor(ii);
                
                if abs(maxdE*factor(ii)) > 3 && r > 30
                    Eps(indx2) = maxdE;
                    sns2(indx2) = sn;
                    indx2 = indx2 + 1;
                end
                
                indx1 = indx1 + 1;
                
            end
        end
    end
    
    % Lets remove NaNs
    if indx1 < 11
        sns(indx1:10) =[];
        ts(indx1:10) =[];
        Eps_real(indx1:10) =[];        
        Eps_dEdt(indx1:10) =[];
    end
    Eps_hilb(indx3:10) =[];
    sns3(indx3:10) = [];
    
    if indx2 < 11
        Eps(indx2:10) = [];
        sns2(indx2:10) = [];
    end
    
    arg.inner = sns;
    sns_org = sns'; % To save peak info
    arg.outer = [];
    arg.t_in = ts;
    arg.t_out  = [];
    arg.method = 6;
    
    % Total numgber of sensors used
    Npbfa = length(sns);
    
    % request more than 4 sensors to calculate pbfa points
    if Npbfa > 4
        [xpbfa, ypbfa, zpbfa, tpbfa] = pbfa_finder(arg);
        %fprintf('%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.6f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
        %        count,tpbfa,xpbfa,ypbfa,zpbfa,0,Npbfa,tcg,xcg,ycg,Ip,Ncg)
        
        if ~isreal(zpbfa)
            zpbfa = 0;
        end
        
        ki_sqrd = cal_ki_sqrd(tpbfa,xpbfa,ypbfa,zpbfa,arg,settings);
        
        
        if ki_sqrd < 3 
            
            Ipbfa = findPeakCurr(settings,xpbfa,ypbfa,zpbfa,Eps_hilb,sns3);
            sen_used = sprintf('%i,',sns);
            count = count + 1;
            fprintf('%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.1f\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
                count,tpbfa,xpbfa,ypbfa,zpbfa,Npbfa+arg.method*10,Ipbfa,sen_used(1:end-1),ki_sqrd,tcg,xcg,ycg,Ip,Ncg)
            
            fprintf(fID,'%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.1f\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
                count,tpbfa,xpbfa,ypbfa,zpbfa,Npbfa+arg.method*10,Ipbfa,sen_used(1:end-1),ki_sqrd,tcg,xcg,ycg,Ip,Ncg);
            
            store_peaks_info(count,[ts', sns'],fID2)
            store_peaks_info(count,[Eps_real', sns_org],fID3)
            store_peaks_info(count,[Eps_hilb', sns3'],fID4)
            store_peaks_info(count,[Eps_dEdt', sns_org],fID5)
            
        elseif Npbfa > 5 && ki_sqrd > 3;
            removes = remove_check(arg,settings);
            
            removes  = [removes.s_ki, removes.s_zs, sns'];
            
            removes( removes(:,1) > 5,:) = [];
            removes( removes(:,2) > 3000,:) = [];
            
            if ~isempty(removes)
                
                [~, indx] = min(removes(:,1));
                indx = find(sns == removes(indx,3));
                sns(indx) = [];
                ts(indx) = [];
                
                arg.t_in = ts;
                arg.inner = sns;
                
                [xpbfa, ypbfa, zpbfa, tpbfa] = pbfa_finder(arg);
                ki_sqrd = cal_ki_sqrd(tpbfa,xpbfa,ypbfa,zpbfa,arg,settings);
                Ipbfa = findPeakCurr(settings,xpbfa,ypbfa,zpbfa,Eps_hilb,sns3);
                sen_used = sprintf('%i,',sns);
                count = count + 1;
                Npbfa = length(sns);
                
                fprintf('%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.1f\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
                    count,tpbfa,xpbfa,ypbfa,zpbfa,Npbfa+arg.method*10,Ipbfa,sen_used(1:end-1),ki_sqrd,tcg,xcg,ycg,Ip,Ncg)
                
                fprintf(fID,'%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.1f\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\n',...
                    count,tpbfa,xpbfa,ypbfa,zpbfa,Npbfa+arg.method*10,Ipbfa,sen_used(1:end-1),ki_sqrd,tcg,xcg,ycg,Ip,Ncg);
            
                store_peaks_info(count,[ts', sns'],fID2)
                store_peaks_info(count,[Eps_real', sns_org],fID3)
                store_peaks_info(count,[Eps_hilb', sns3'],fID4)
                store_peaks_info(count,[Eps_dEdt', sns_org],fID5)
            
            end
        else
            fprintf('Coldn''t located the PBFA\n')            
        end
    end
end

delete(wbh);
fclose(fID);
fclose(fID2);
fclose(fID3);
fclose(fID4);
fclose(fID5);
%toc


function [ch_fn, h_fn] = generate_ch_fn(settings,g)
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


function Ip = findPeakCurr(sen_set,x,y,z,maxdV,sns)

% %slope
% m = [ 4.8930 3.9805 4.1938  NaN     2.1098 ...
%       4.2255 NaN    2.0474  1.6661  2.0610  NaN];
% %intercept
% c = [0.0382 0.0638  0.2847  NaN     1.2001 ...
%      0.5778 NaN     1.8081  1.0528  1.1174 NaN];
% 
%  
% r = sqrt((sen_set.x - x).^2 + ...
%          (sen_set.y - y).^2 + ...
%          (sen_set.z - z).^2)/1000;
% 
% Ip = nanmean( m(sns).* maxdV .* (r(sns).^1.13) + c(sns));

maxdV1 = nan(1,11);
maxdV1(sns) = maxdV;
m =  [472.4 387.1 418.7 NaN NaN 411.9 676.4  204.2 165.6 209.0 NaN];
r = sqrt((sen_set.x(1:11) - x).^2 + ...
         (sen_set.y(1:11) - y).^2 + ...
         (sen_set.z(1:11) - z).^2)/1000;

Ip = m.*maxdV1.*((r/100).^1.13);
indx = find(r > 30);
%indx = intersect(indx,sns)
Ip = -nanmean(Ip(indx));

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

function removes = remove_check(arg,sen_set)
    
    L = length(arg.t_in);
    
    if L > 5
        removes.s_ki = NaN(L,1);
        removes.s_zs = NaN(L,1);
        
        %arg.method = 5;
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
    
    
   

