function battery_monitor_post

fn = 'C:\Users\sumedhe\Desktop\BatteryVoltageData05.txt';

% battery numbers
bat = [6 5 8 7];

fID = fopen(fn);
data = textscan(fID,'%f-%f-%f %f:%f:%f %f %f %f %f');
fclose(fID);

[fn fn] = fileparts(fn);

% Voltage
% Seconds
t = data{3}*86400+data{4}*3600+data{5}*60+data{6};
V1 = data{7};
V2 = data{8};
V3 = data{9};
V4 = data{10};

% Remove 0s
V1(V1 == 0) = NaN;
V2(V2 == 0) = NaN;
V3(V3 == 0) = NaN;
V4(V4 == 0) = NaN;

% Connected Resistor
R = 2.5;

% Published battery capacity in amphere hour
Q = 75;

% Strart time
ind1 = find(V1>0); tSt1 = t(ind1(1));
ind2 = find(V2>0); tSt2 = t(ind2(1));
ind3 = find(V3>0); tSt3 = t(ind3(1));
ind4 = find(V4>0); tSt4 = t(ind4(1));


%Remove < 10.5V
ind1 = find(V1(ind1(1):end) < 10.5);
ind2 = find(V2(ind2(1):end) < 10.5);
ind3 = find(V3(ind3(1):end) < 10.5);
ind4 = find(V4(ind4(1):end) < 10.5);

if isempty(ind1); ind1 = length(V1); end
if isempty(ind2); ind2 = length(V2); end
if isempty(ind3); ind3 = length(V3); end
if isempty(ind4); ind4 = length(V4); end

V1(ind1(1):end) = NaN;
V2(ind2(1):end) = NaN;
V3(ind3(1):end) = NaN;
V4(ind4(1):end) = NaN;

% Battery run time
tE1 = t(ind1(1))-tSt1;
tE2 = t(ind2(1))-tSt2;
tE3 = t(ind3(1))-tSt3;
tE4 = t(ind4(1))-tSt4;


% Currnt battery Capacity (Ah)
dts = t(2:end)-t(1:end-1);
Q1 = nansum( V1(2:end)/R.*dts/3600);
Q2 = nansum(V2(2:end)/R.*dts/3600);
Q3 = nansum(V3(2:end)/R.*dts/3600);
Q4 = nansum(V4(2:end)/R.*dts/3600);


% Plot data
figure
box on
hold all
plot((t-tSt1)/3600,V1)
plot((t-tSt2)/3600,V2)
plot((t-tSt3)/3600,V3)
plot((t-tSt4)/3600,V4)
ylim([0,13])
plot(xlim,[10.5 10.5],'k','linewidth',2)
title(sprintf('%s     %4.4i-%2.2i-%2.2i',fn,data{1}(1),data{2}(1),data{3}(1)))
xlabel('Time (hours)')
ylabel('Voltage (V)')
legend(sprintf('BAT%3.3i %0.1f Ah (%0.1f%%) %s h',bat(1),Q1,Q1*100/Q,sec2hhmm(tE1)),...
    sprintf('BAT%3.3i %0.1f Ah (%0.1f%%) %s h',bat(2),Q2,Q2*100/Q,sec2hhmm(tE2)),...
    sprintf('BAT%3.3i %0.1f Ah (%0.1f%%) %s h',bat(3),Q3,Q3*100/Q,sec2hhmm(tE3)),...
    sprintf('BAT%3.3i %0.1f Ah (%0.1f%%) %s h',bat(4),Q4,Q4*100/Q,sec2hhmm(tE4)),...
    'GoPower 175W Lowest input voltage',...
    'Location','SouthWest')

tools2fig

function hhmm = sec2hhmm(t)

hh = floor(t/3600);
mm = round((t-hh*3600)/60);
hhmm = sprintf('%2.2i:%2.2i',hh,mm);

