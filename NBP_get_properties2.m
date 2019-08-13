function NBP_get_properties2
% After loading press
%   q - to quit
%   c - to collect data (collect rise time, fall time, etc)
%   n - to go to next plot (ignore the current plot)
%   t - to try next available sensor for this pulse
%   z - correct zero crossing time
%   e - correct end time
%   r - correct rise peak
%   d - Done! record data in the file
%   f - correct fall peak

% Modification History
%   2015-03-23 Copied from NBP_get_properties and changed to collect total
%   NBP duration from different sensors.

% Open the NBP file
%data = xlsread('C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814-NBP_info2.xlsx');
data = xlsread('C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\20110814-NBP_info3.xlsx');

% Starting index acoording to column 1
index = 1;

% base foloder for saving plots
bf = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\TotalTime';

% data saving file name
dsfn = [bf '20110814-time_info.txt'];


%% Start the program
% Prepare writing file
if ~exist(dsfn,'file')
    fID = fopen(dsfn,'a+');
    fprintf(fID,'%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n\n',...
        'Index','K02','K14','K24','BCC','K17','EDW','STC','FLT','OVD','FFI');
else
    fID = fopen(dsfn,'a+');
end


% Load data
% get plotter2 data
h=guidata(findall(0,'Tag','plotter2'));
sen_set = h.sen_set;
g = h.g;

% Turn off all plots
g.lpgraphs = zeros(1,60);
g.chgraphs = zeros(1,60);

cc = '';

while ~strcmp(cc,'q')
    
    % to ingnore header we need to add 2 lines to raws
    t = data(index+2,10);
    INDEX = data(index+2,1)
    
    if ~isnan(t)
        
        plot_ind = 1;
        % Plot data
        [sID, ch_data, a] = plot_data(g,t,sen_set,plot_ind,data(index+2,11:13));
        fg = gcf;
        
        
        while ~ismember(cc,{'q','c','n','t','z','e','r','d','f'})
            figure(fg)
            cc = get(fg,'CurrentCharacter');
            pause(0.1)
            
            % Try next sensor
            if strcmp(cc,'t')
                plot_ind = plot_ind + 1;
                delete(fg)
                [sID, ch_data,a] = plot_data(g,t,sen_set,plot_ind,data(index+2,11:13));
                fg = gcf;
                cc = '';
                
                % Correct zero cross time
            elseif strcmp(cc,'z')
                [zct,zcv] = get_point(ch_data.to, ch_data.yo);
                delete(a.h(3));
                a.h(3) = plot(zct,zcv,'go','markerfacecolor','g');
                a.zCrossTime = (zct - a.st_t)*1e6;
                a = printText(a);
                set(fg,'CurrentCharacter','a')
                cc = '';
                
                % Correct end time
            elseif strcmp(cc,'e')
                [zct,zcv] = get_point(ch_data.t, ch_data.y);
                delete(a.h(5));
                a.h(5) = plot(zct,zcv,'go','markerfacecolor','g');
                a.totTime = (zct - a.st_t)*1e6;
                a = printText(a);
                set(fg,'CurrentCharacter','a')
                cc = '';
                
                % Correct positive peak time
            elseif strcmp(cc,'r')
                [zct,zcv,ind] = get_point(ch_data.to, ch_data.yo);
                delete(a.h(2));
                a.h(2) = plot(zct,zcv,'go','markerfacecolor','g');
                a.riseTime = (zct - a.st_t)*1e6;
                % get new 10-90 rise time
                a.mx_t = zct;
                a.mx_y = zcv;
                a.mx_ind = ind;
                
                a = get_FWHM(ch_data.to, ch_data.yo,a);
                a = z90risetime(ch_data.to, ch_data.yo,a);
                
                a = printText(a);
                set(fg,'CurrentCharacter','a')
                cc = '';
                
                % Correct negative peak
            elseif strcmp(cc,'f')
                [zct,zcv] = get_point(ch_data.t, ch_data.y);
                delete(a.h(4));
                a.h(4) = plot(zct,zcv,'go','markerfacecolor','g');
                a.negPeakT = (zct - a.st_t)*1e6;
                a = printText(a);
                set(fg,'CurrentCharacter','a')
                cc = '';
                
                
                % Collect data
            elseif strcmp(cc,'c')
                collect_data(index,sID,ch_data);
                
                
                % Go to next plot means we do not record the data for this plot
            elseif strcmp(cc,'n')
                fprintf(fID,'%10.3i\t%10.1f\t%10.1f\t%10.1f\t%10.1f\t%10.1f\t%10.1f\t\n',...
                    index,NaN,NaN,NaN,NaN,NaN,NaN);
                
                
                % Done (means record the current data)
            elseif strcmp(cc,'d')
                % Print info to the file
                fprintf(fID,'%10.3i\t%10.1f\t%10.1f\t%10.1f\t%10.1f\t%10.1f\t%10.1f\t\n',...
                    index,a.riseTime,a.z90_riseTime,a.zCrossTime,a.negPeakT,a.totTime,a.FWHM);
                
                % save the figure
                saveas(fg,sprintf('%s%4.4i-timing.fig',bf,index))
                saveas(fg,sprintf('%s%4.4i-timing.png',bf,index))
                
                
            end
            
        end
        
        % reset the character
        if ~strcmp(cc,'q')
            cc = '';
        end
        
        % Delete the current figure
        delete(fg)
        
    end

    index = index + 1;
    
end

fclose(fID)
disp('done')




function [sID, ch_data, a] = plot_data(g,t,sen_set,sn_ind,xyz)

[hms h m] = sec2hhmmss(t);

g.hh = h;
g.mm = floor(m/5)*5;
g.ss = 0;

% Let's plot closest sensor available after 10 km
%sns = [1 2 3 6 8 9 10 11];

% Horizontal distance
r = sqrt((sen_set.x-xyz(1)).^2+(sen_set.y-xyz(2)).^2);

% Do not use anything close 15km
%r(r < 25000) = inf;

% select the closest {sn_ind}th sensor after 15 km
%[~,I] = sort(r);

try
    g.chgraphs(sn_ind*3) = 1;
catch
    disp('No more sensors')
    figure
    return
end

% find begining and end time
g.t1 = t + r(sn_ind)/3e8 - 0.00005;
g.t2 = g.t1 + 0.005;

sID = sn_ind;
sen_set.ldar_tshift_sn = sID;
save('sensor_setting.mat','-Struct','sen_set')

data = plot_all5(g,0);

ch_data.yo = [];
a.txt_h = [];

% collect some info
if isfield(data,'ch_data') && range(data.ch_data(sn_ind*3).y) > 0.1
    
    % Original data
    yo = data.ch_data(sn_ind*3).y;
    to = data.ch_data(sn_ind*3).t;
    
    ch_data.yo = yo;
    ch_data.to = to;
    
    % Let's plot 1 MHz data
    L = floor(length(yo)/5)*5;
    y = zeros(1,L/5);
    t = y;
    for i = 1:5
        y = y + yo(i:5:L+i-1);
        t = t + to(i:5:L+i-1);
    end
    y = y/5;
    t = t/5;
    
    ch_data.y = y;
    ch_data.t = t;
    plot(t,y)
    
    
    % Let's find out 0V horizontal line
    [mm ind] = max(y);
    zeroV1 = mean(y(1:ind-30));
    plot(t,t-t+zeroV1)
    
    zeroV2 = mean(y(ind+100:end));
    stdV2 = std(y(ind+100:end));
    
    plot(t,t-t+zeroV2)
    
    % + peak time
    [mx_y, ind] = max(yo);
    mx_t = to(ind);
    mx_ind = ind;
    
    a.h(2) = plot(mx_t,mx_y,'ko','markerfacecolor','k');
    
    %Get start time
    for i = 1:100
        if yo(ind-i) < zeroV1
            % found it
            break
        end
    end
    [st_t st_y] = crosspoint(to(ind-i),yo(ind-i),to(ind-i+1),yo(ind-i+1),...
        to(1),zeroV1,to(end),zeroV1);
    
    a.h(1) = plot(st_t,st_y,'ko','markerfacecolor','k');
    
    % Get zero crossing time
    for i = 1:100
        if yo(ind+i) < zeroV1
            % found it
            break
        end
    end
    
    [zc_t,zc_y] = crosspoint(to(ind+i),yo(ind+i),to(ind+i-1),yo(ind+i-1),...
        to(1),zeroV1,to(end),zeroV1);
    a.h(3) = plot(zc_t,zc_y,'ko','markerfacecolor','k');
    
    % - Peak time (use 1 us data)
    [mm ind] = min(y);
    mn_t = t(ind);
    mn_y = mm;
    a.h(4) = plot(mn_t,mn_y,'ko','markerfacecolor','k');
    
    
    % end time (use 1us data) (find point within 1STD)
    diffs = abs(y(ind:end)-zeroV2);
    
    for i = 1:length(diffs)
        if diffs(i) < stdV2/5
            % found it
            break
        end
    end
    
    en_t = t(ind+i-1);
    en_y = y(ind+i-1);
    a.h(5) = plot(en_t,en_y,'ko','markerfacecolor','k');
    
    
    
    
    
    % Basic info
    a.st_t = st_t;
    a.mn_t = mn_t;
    a.zeroV1 = zeroV1;
    a.mx_t = mx_t;
    a.mx_y = mx_y;
    a.mx_ind = mx_ind;
    
    % FWHM (full width at half maximum)
    a = get_FWHM(to,yo,a);
    
    % Obtain the 10-90% rise time
    a = z90risetime(to,yo,a);
    
    % Calculate quantities
    a.riseTime = (mx_t - st_t)*1e6;
    a.zCrossTime = (zc_t - st_t)*1e6;
    a.totTime = (en_t - st_t)*1e6;
    a.negPeakT = (mn_t - st_t)*1e6;
    
    a.txt_h = [];
    a = printText(a);
    
end

set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% Let's print the distance
fprintf('Horizontal distance = %0.2f km\n',r(sn_ind)/1000)


function a = printText(a)
delete(a.txt_h);
%str = sprintf('FWHM = %0.1f us\n10-90 Rise Time = %0.1f us\nRise Time = %0.1f us\n0-Cross Time = %0.1f us\nTotal Time = %0.1f us\n',...
%           a.FWHM,a.z90_riseTime,a.riseTime,a.zCrossTime,a.totTime);

str = sprintf('%s = %4.1f us\n%s = %4.1f us\n%s = %4.1f us\n%s = %4.1f us\n%s = %4.1f us\n%s = %4.1f us\n',...
    'RiseTime  ',a.riseTime,'10-90RiseT',a.z90_riseTime,'ZeroCrossT',a.zCrossTime,...
    'NegPeakT  ',a.negPeakT,'TotalTime ',a.totTime,'FWHM      ',a.FWHM);

a.txt_h = text(a.mn_t,a.fwhm_y1,str,'fontsize',20);


function a = get_FWHM(to,yo,a)

% delete previous rise time info
try
    delete(a.h(8:10))
end

% + peak time
mm = a.mx_y;
ind = a.mx_ind;
%[mm ind] = max(yo);

%Get first half time
for i = 1:100
    if yo(ind-i) < (a.zeroV1+mm)/2
        % found it
        break
    end
end

[fwhm_t1, fwhm_y1] = crosspoint(to(ind-i),yo(ind-i),to(ind-i+1),yo(ind-i+1),...
    to(1),(a.zeroV1+mm)/2,to(end),(a.zeroV1+mm)/2);

a.h(8) = plot(fwhm_t1,fwhm_y1,'ro');

% Get the second half time
for i = 1:100
    if yo(ind+i) < (a.zeroV1+mm)/2
        % found it
        break
    end
end

[fwhm_t2, fwhm_y2] = crosspoint(to(ind+i),yo(ind+i),to(ind+i-1),yo(ind+i-1),...
    to(1),(a.zeroV1+mm)/2,to(end),(a.zeroV1+mm)/2);

a.h(9) = plot(fwhm_t2,fwhm_y2,'ro');

a.h(10) = plot([fwhm_t1 fwhm_t2],[fwhm_y1 fwhm_y2],'linewidth',2);

a.fwhm_y1 = fwhm_y1;
a.FWHM = (fwhm_t2-fwhm_t1)*1e6;



function a = z90risetime(to,yo,a)

% delete previous rise time info
try
    delete(a.h(6:7))
end


ind = a.mx_ind;

% get the 10 % time
ten_tre = a.zeroV1 + (a.mx_y-a.zeroV1)*0.1;

for i = 1:100
    if yo(ind-i) < ten_tre
        % found it
        break
    end
end

[ten_t, ten_y] = crosspoint(to(ind-i),yo(ind-i),to(ind-i+1),yo(ind-i+1),...
    to(1),ten_tre,to(end),ten_tre);

a.h(6) = plot(ten_t,ten_y,'ro');

nin_tre = a.zeroV1 + (a.mx_y-a.zeroV1)*0.9;

% get the 90% time
for i = 1:100
    if yo(ind-i) < nin_tre
        % found it
        break
    end
end

[nin_t, nin_y] = crosspoint(to(ind-i),yo(ind-i),to(ind-i+1),yo(ind-i+1),...
    to(1),nin_tre,to(end),nin_tre);

a.h(7) = plot(nin_t,nin_y,'ro');

a.z90_riseTime = (nin_t-ten_t)*1e6;


function h = collect_data(index,sID,ch_data)


t = ch_data.to;
E = ch_data.yo;

% handles to points
h = nan(1,4);
ts = h;
Es = h;

for i = 1:4
    
    [t0 E0] = get_point(t,E);
    
    h(i) = plot(t0,E0,'ko');
    ts(i) = t0;
    Es(i) = E0;
end

% % Save the data in a file
% fID = fopen('C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\NBP_physical_properties.txt','a+');
% fprintf(fID,'%4.4i\t%3.1i\t%13.6f\t%13.6f\t%13.6f\t%13.6f\t%5.1f\t%5.1f\t%5.1f\t%5.1f\n',index,sID,ts,Es);
% fclose(fID);

% Let's print some data as some feedback
fprintf('\nRise time: %5.1f us\n',(ts(2)-ts(1))*1e6)
fprintf('Fall time: %5.1f us\n',(ts(4)-ts(2))*1e6)
fprintf('Totl time: %5.1f us\n\n',(ts(4)-ts(1))*1e6)


function [x y ind] = get_point(t,E)

[t0, E0] = ginput(1);

% Find closest point and snap it
D = sqrt(((t-t0)/range(xlim)).^2+((E-E0)/range(ylim)).^2);

[mm ind] = min(D);
x = t(ind);
y = E(ind);





function [x,y] = crosspoint(x1,y1,x2,y2,x3,y3,x4,y4)

x = ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
y = ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));




