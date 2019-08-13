function g2plotter(g)


% get plotter2 data
h=guidata(findall(0,'Tag','plotter2'));

%handle to t1 and t2 in plotter
ht1=findall(0,'Tag','t1');
ht2=findall(0,'Tag','t2');

h.g.t1=g.t1;
h.g.t2=g.t2;

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