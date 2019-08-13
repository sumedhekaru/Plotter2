function find_max_dE_vs_Ip
% this function will find maxdE vs Ip for all 1s fast-slow antenna data

output_fn = 'C:\Users\sumedhe\Desktop\Ip_Ep_data_20110814_1800-2400_HT.txt';
fID = fopen(output_fn, 'a+');
fprintf(fID,'\n');

% Import plotter data
try; h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

settings = h.sen_set;
g = h.g;

    cglss_fn = sprintf('%s/cglss/%s/%s/KSCCGLSS%s%s%s.dat',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},g.MM{:},g.DD{:});
x0 = 0;
y0 = 0;
z0 = 0;

% Load CGLSS data
cglss = CGLSS_extract(cglss_fn,61200,86400,str2double(settings.ldar_r),x0,y0,z0);

tshift=settings.t_shift;

factor = g.factor;

ch_legend=cell(1,60);
for i=1:20;   
    ch_legend{i*3-2}=[settings.sen_IDs{i} ':ch1'];
    ch_legend{i*3-1}=[settings.sen_IDs{i} ':ch2'];
    ch_legend{i*3}=[settings.sen_IDs{i} ':ch3'];    
end




% Total number of CGLSS sources
tot_CG = length(cglss(:,1));


try
    wbh = waitbar(0,'Loading Fast Antenna data','Name','Collecting Peak Currents');
catch
    return
end

Ips = NaN(1,tot_CG*10);
rs  = Ips;
Eps = Ips;
count = 0;

for k = 1:tot_CG
   

    
    tcg = cglss(k,1);
    xcg = cglss(k,2);
    ycg = cglss(k,3);
    Ip = cglss(k,11);
 
    % Generate ch file names according to the current CGLSS poing
    hhmmss = sec2hhmmss(tcg);
    g.hh = str2double(hhmmss(1:2));
    g.mm = floor(str2double(hhmmss(4:5))/5)*5;
    
    [ch_fn h_fn] = generate_ch_fn(settings,g);
    
    
    % Load fast antenna data according to the current CGLSS point
    
    g.t1 = tcg - 0.5e-3;
    g.t2 = tcg + 0.5e-3;
    
    %     figure
    %     hold all
    
    
    
    for i=1:60
        
        try
            val = (k-1)/tot_CG + 1/tot_CG*i/60;
            waitbar(val,wbh,sprintf('Loading Fast Antenna data. (%.2f%%)',val*100))
        catch
            return
        end
        
        
        if strcmp(ch_fn{i},'')==0
            
            [t y ch_freq_str] = FA_Extract1(ch_fn{i},h_fn{i},g.t1,g.t2,tshift(i),settings,i);
            
            % distance
            sn = ceil(i/3);
            r = sqrt((settings.x(sn)-xcg)^2 + (settings.y(sn)-ycg)^2+ settings.z(sn)^2)/1000;
            
            if ~isempty(t)
                % Derivative of voltage (or the difference of the voltage)
                dy = y(2:end)-y(1:end-1);
                               
                %plot(t,y*factor(i))
                %plot(t(1:end-1)+dt/2,dy*factor(i))
                                
                maxdE=min(dy)*factor(i);                
                maxdE_name= [ch_legend{i} ch_freq_str];
                
                if abs(maxdE) > 1 && r > 5 && r < 100
                    
                    count = count + 1;
                    
%                     fprintf(fID,'%i\t%0.6f\t%s\t%0.2f\t%0.1f\t%0.1f\n',...
%                         count,tcg,maxdE_name, maxdE,r,Ip);
                    fprintf('%i\t%0.6f\t%s\t%0.2f\t%0.1f\t%0.1f\n',...
                        count,tcg,maxdE_name, maxdE,r,Ip)
                    
                    rs(count) = r;
                    Eps(count) = maxdE;
                    Ips(count) = Ip;
                    
                end         
                
            end
        end
    end
    
end

fclose(fID);
delete(wbh)
N = [1 1.3];
peakCurrBestFit(Ips(1:count),Eps(1:count),rs(1:count),N)


function [ch_fn h_fn] = generate_ch_fn(settings,g)
% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

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

