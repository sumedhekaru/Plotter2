function load_data2plotter(indx)
% This function will read time info from excel file and then load in to
% plotter2 inputs. Note that you should enter date and plottings manually.
% This will only enter hh, mm, ss and t1 and t2

% Open the NBP file
[data,TXT,RAW] = xlsread('C:\Users\sumedhe\Desktop\NBP-timing-JGR-2015\20110814-NBP_info2.xlsx');
t = data(indx+2,10);

% [data,TXT,RAW] = xlsread('C:\Users\sumedhe\Downloads\20110722-NBP_info.xlsx');
% indx = find(data(:,1)==indx);
% t = data(indx,4);

if isnan(t)
    disp('Data not loaded. No time info found')
    return
end

fprintf('%0.7f\n',t)
xlimits = [t-2e-3,t+2e-3];

% get plotter2 data
h=guidata(findall(0,'Tag','plotter2'));

%handle to t1 and t2 in plotter
ht1=findall(0,'Tag','t1');
ht2=findall(0,'Tag','t2');

h.g.t1=xlimits(1);
h.g.t2=xlimits(2);

t1=sprintf('%.6f',h.g.t1);
t2=sprintf('%.6f',h.g.t2);

set(ht1,'String',t1)
set(ht2,'String',t2)

% Updating hh, mm
hh = floor(h.g.t1/3600);
mm = floor((h.g.t1 - hh*3600)/300)*5;

hh = sprintf('%02i',hh);
l =get(h.hour,'String');
value = find(strcmp(hh,l));

h.g.hhn = value;
h.g.hh  = str2double(l(value));
set(h.hour,'Value',value)


mm = sprintf('%02i',mm);
l =get(h.minute,'String');
value = find(strcmp(mm,l));

h.g.mmn = value(1);
h.g.mm  = str2double(l(value(1)));
set(h.minute,'Value',value(1));

% Set time maximum and minimum values
t1minimum=h.g.hh*3600+60*h.g.mm;
t2maximum=t1minimum+300;

h.g.t1min=t1minimum;
h.g.t2max=t2maximum;

t1min_s=sprintf('t Min = %is',t1minimum);
t2max_s=sprintf('t Max = %is',t2maximum);
set(h.t1min,'String',t1min_s)
set(h.t2max,'String',t2max_s)

guidata(findall(0,'Tag','plotter2'), h)