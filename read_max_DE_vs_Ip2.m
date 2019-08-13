function read_max_DE_vs_Ip2

%% USER INPUTS
sn = 1;                 % Sensor number
r1 = 30000;             % lower r limit
r2 = 100000;            % heigher r limit
rd = 1;                % plot real data?
ht = 1;                 % Plot Hibert transform data?
bd = 1;                 % Plot Brushed data

% When bd = 1, this will just get the data from the current figure and
% replot. This is good for brushing and removing outliers. WHen this switch
% on, the program will not go beyond that.

%%

if bd 
    h = findobj(gca,'Type','line');
    x=get(h,'Xdata');
    y=get(h,'Ydata');
    tit = get(get(gca,'title'),'string');
    
    
    x = x{4};
    y = y{4};
    

    ind = isnan(y);
    x(ind) = [];
    y(ind) = [];
    
    hh = best_fit_line(x,y);
    
    tit = sprintf('%s%0.4f',tit(1:end-6),hh.corrCoeff);
    title(tit,'FontWeight','bold')
    xlabel('V_pR^{1.13}')
    ylabel('CGLSS I_p (kA)') 

    box on
    grid off

    return
end


%%
sen_set = open('sensor_setting.mat');

fn = 'C:\Users\sumedhe\Desktop\Ip_Ep_data_20110814_1800-2400_HT.txt';

fID = fopen(fn,'r');
data = textscan(fID,'%f %f %s %f %f %f %f','Headerlines',0);
fclose(fID);

% 1 - index
% 2 - time
% 3 - sensor name
% 4 - Peak Voltage
% 5 - peak HT
% 6 - distance in km
% 7 - peak Current

% Copy the sensor ID column
senID = data{3};
data{3} = data{1};
data = cell2mat(data);

% Get the data for the interesting sensor
switch sn
    case 1; index = strcmp(senID,'K02:ch3:S-1.00MHz');
    case 2; index = strcmp(senID,'K14:ch3:S-1.00MHz');
    case 3; index = strcmp(senID,'K24:ch3:S-1.00MHz');
    case 5; index = strcmp(senID,'BCC:ch3:O-1.00MHz');
    case 6; index = strcmp(senID,'K17:ch3:S-1.00MHz');
    case 7; index = strcmp(senID,'EDW:ch3:S-1.00MHz');
    case 8; index = strcmp(senID,'STC:ch3:S-1.00MHz');
    case 9; index = strcmp(senID,'FLT:ch3:S-1.00MHz');
    case 10; index = strcmp(senID,'OVD:ch3:S-1.00MHz');
    otherwise; disp('Wrong sensor number'); return
end

data = data(index,:);

% choose data between r1 and r2
data = sortrows(data,6);
lol = sum(data(:,6) < r1/1000) + 1;
ul  = sum(data(:,6) < r2/1000);
data = data(lol:ul,:);

if rd
    
    x = abs(data(:,4)).*data(:,6).^1.13;
    y = abs(data(:,7));
    
    %figure
    %plot(x,y,'ro','markerfacecolor','r','markersize',2)
    
    out = best_fit_line(x,y);
    title(sprintf('%s Real Data    %0.0f - %0.0fkm     CorCoeff = %0.4f',...
        sen_set.sen_IDs{sn},r1/1000,r2/1000,out.corrCoeff),'FontWeight','bold')
    box on
    grid off
end


if ht
    
    x = abs(data(:,5)).*data(:,6).^1.13;
    y = abs(data(:,7));
    
    %figure
    %plot(x,y,'ro','markerfacecolor','r','markersize',2)
    
    out = best_fit_line(x,y);
    title(sprintf('%s HT Data    %0.0f - %0.0fkm     CorCoeff = %0.4f',...
        sen_set.sen_IDs{sn},r1/1000,r2/1000,out.corrCoeff),'FontWeight','bold')
    box on
    grid off
end



