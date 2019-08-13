function IBP_to_RS_ratio
clc
% This function is intended to find IBP to RS ratio of flashes. This
% function will use CG data as reference to identify the first CG source.

% Load plotter2 data
try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end


% save files in this folder
bf = 'C:\Users\sumedhe\Desktop\IBP_to_RS_ratio_auto2\';

% starting index
indx = 0;  % if program stopped unexpectly, just put the last index you see in abouve bf folder

settings = h.sen_set;
g = h.g;
tshift=settings.t_shift;
factor = g.factor;

cglss_fn = sprintf('%s/cglss/%s/%s/KSCCGLSS%s%s%s.dat',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},g.MM{:},g.DD{:});
x0 = 0;
y0 = 0;
z0 = 0;

% Load CGLSS data
cglss = CGLSS_extract(cglss_fn,43200,86400,2000000,x0,y0,z0);
L = length(cglss);

% find first RS of flashes
dt = cglss(2:L,1)-cglss(1:L-1,1);
peakLoc = [];
dt_th = 1; % RS interval threshold
for i=1:length(dt)
   if dt(i) > dt_th
       peakLoc = [peakLoc i+1];
   end
end



L = length(peakLoc);

% hold all
% plot(cglss(:,1),cglss(:,1)-cglss(:,1)-650,'ro','markerfacecolor','r')
% plot(cglss(inds,1),cglss(inds,1)-cglss(inds,1)-675,'ro','markerfacecolor','r')


% Turn on all ch3 graphs (avoid WSB and BCC)
% g.chgraphs = [0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 ...
%      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

g.chgraphs = [1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0  ...
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

%g.chgraphs = zeros(1,60);
%g.chgraphs(4) = 1; % K14 ch1

wbh = waitbar(0,'Loading Fast Antenna data','Name','IBP to RS ratio finder');

% setup figure
% delete the previosu figure
try; delete(145231); end
figure(145231)
set(145231,'Position',get(0,'ScreenSize'))
tools2fig

for k = 2:L
    
    tcg = cglss(peakLoc(k)-1,1);
    xcg = cglss(peakLoc(k)-1,2);
    ycg = cglss(peakLoc(k)-1,3);
    
    % Generate ch file names according to the current CGLSS poing
    hhmmss = sec2hhmmss(tcg);
    g.hh = str2double(hhmmss(1:2));
    g.mm = floor(str2double(hhmmss(4:5))/5)*5;
    
    [ch_fn, h_fn] = generate_ch_fn(settings,g);
    
    t_shifts =sqrt((xcg - settings.x).^2 + (ycg - settings.y).^2)/3e8;
    
    g.t1 = tcg - 50e-3;
    g.t2 = tcg + 0.1e-3;
    
    for i=1:1:60 % only get ch3 graphs
        
        try
            val = (k-1)/L + 1/L*i/60;
            waitbar(val,wbh,sprintf('Loading Fast Antenna data. (%.2f%%)',val*100))
        catch
            return
        end
        
        sn = ceil(i/3);
        
        if strcmp(ch_fn{i},'')==0
            
            % Let's load ch3 data
            [t, y] = FA_Extract1(ch_fn{i},h_fn{i},g.t1,g.t2,tshift(i)- t_shifts(sn),settings,i);
            
            % calibration
            y = y*factor(i);
            
            % Sample length
            sL = length(t);
            
            if sL > 3000
                
                indx = indx + 1/60;
                
                % offset value
                ofst = mean(y(1:sL-2000));
                noice = std(y(1:sL-2000));
                
                % Obtaining RS peak
                tRS = t(sL - 150:sL);
                yRS = y(sL - 150:sL);
                
                [RS_minY , RS_min_ind] = min(yRS);
                RS_mint = tRS(RS_min_ind);
                if RS_min_ind > 10; RS_minYR =  range(yRS(RS_min_ind-10:RS_min_ind));
                else                RS_minYR =  range(yRS(1:RS_min_ind)); end;
                
                
                % Obtaining IBP peak
                tIBP = t(1:sL-2000);
                yIBP = y(1:sL-2000);
                
                [IBP_minY , IBP_min_ind] = min(yIBP);
                IBP_mint = tIBP(IBP_min_ind);
                if IBP_min_ind > 20;  IBP_minYR =  range(yIBP(IBP_min_ind-20:IBP_min_ind));
                else                  IBP_minYR =  range(yIBP(1:IBP_min_ind));   end
                
                % Obtaining IBP peak to peak
                if IBP_min_ind < 50;      ind1 = 1;
                else ind1 = IBP_min_ind - 50;  end;
                
                LIBP = length(tIBP);
                if IBP_min_ind + 50 > LIBP;      ind2 = LIBP;
                else ind2 = IBP_min_ind + 50;  end;
                
                
                if IBP_min_ind + 20 > LIBP;      ind3 = LIBP;
                else ind3 = IBP_min_ind + 20;  end;
                
                tP_IBP = tIBP(IBP_min_ind:ind3);
                yP_IBP = yIBP(IBP_min_ind:ind3);
                [IBP_maxY , IBP_max_ind] = max(yP_IBP);
                IBP_maxt = tP_IBP(IBP_max_ind);
                
                
                figure(145231)
                % main plot
                subplot(2,3,1:3); cla;
                plot(t,y,'b')
                legend(settings.sen_IDs{sn})
                hold all
                plot(tcg,RS_minY,'ro')
                plot(IBP_mint,IBP_minY,'ro')
                plot([t(1),t(end)],[ofst,ofst],'g')
                plot([t(1),t(end)],[ofst-noice,ofst-noice],'r')
                plot([t(1),t(end)],[ofst+noice,ofst+noice],'r')
                
                
                % let write everything in to the title
                str = sprintf('Indx = %05d-%02d   RS = %07.3f V/m   IBP-np = %07.3f V/m   IBP-pp = %07.3f V/m   Noise = %07.3f   dt = %06.3f ms   t = %12.6f s   x = %7.1f km   y = %7.1f km', ...
                    ceil(indx), sn, RS_minYR, IBP_minYR, range([IBP_minY,IBP_maxY]), noice, (RS_mint - IBP_mint)*1000,tcg, xcg/1000,ycg/1000);
                disp(str)
                title(str)
                
                % IBP plot
                subplot(2,3,5); cla; hold all; box on;
                plot(tIBP(ind1:ind2),yIBP(ind1:ind2),'b')
                plot(IBP_mint,IBP_minY,'ro')
                plot([IBP_mint, IBP_mint],[IBP_minY,IBP_minY + IBP_minYR],'r')
                plot(IBP_maxt,IBP_maxY,'ro')
                plot([IBP_maxt, IBP_maxt],[IBP_minY,IBP_maxY],'r')

                % Rs plot
                subplot(2,3,6); cla; hold all; box on;
                plot(tRS,yRS,'b')
                plot(RS_mint,RS_minY,'ro')
                plot([RS_mint, RS_mint],[RS_minY,RS_minY + RS_minYR],'r')
                
                % Rs location
                subplot(2,3,4); cla;
                daspect([1 1 1]); hold all; box on;
                plot(xcg/1000,ycg/1000,'gs','markerfacecolor','g')
                ns = find(settings.x ~= 0);
                plot(settings.x(ns)/1000,settings.y(ns)/1000,'kp','markerfacecolor','r')
                text(settings.x(ns)/1000,settings.y(ns)/1000,settings.sen_IDs(ns))
                
                fignm = sprintf('%s%05d-%02d',bf,ceil(indx),sn);
                saveas(145231,fignm,'png')
                saveas(145231,fignm,'fig')               
                
            end         
            
        end
    end
    
    indx = ceil(indx);
        
end

delete(wbh)


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

function out = get_peak_info(t,y);
% Obtatining RS peak value

out.RS_peak = min(y(end - 50000:end));

figure(2000)
cla
plot(y(end-2000:end))



