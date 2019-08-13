function sum_plot_ff_data
% clc

fn = 'd:/data/2011/08/13/2300/BCC_20110813_230619.fa2';
d = load(fn,'-mat');

% Begin at 
t1 = 83175.74231142264;
% end tim
t2 = 83175.74234800489;
%np = d.act2.ActualLength;


timenow=d.timenow;
data=d.data;
%size(data)
time = timenow(4)*3600+timenow(5)*60+timenow(6)-(d.PPSt0-d.t0)*1e-6+d.act2.ActualStart/d.freq1:...
       1/d.freq1: ...
       timenow(4)*3600+timenow(5)*60+timenow(6)-(d.PPSt0-d.t0)*1e-6+(d.act2.ActualStart+d.act2.ActualLength-1)/d.freq1;
   
if t1<time(1) || t1 > time(end)
    t1 = time(1);
end

if t2<time(1) || t2 > time(end)
    t2 = time(end);
end

t1
t2

lol=sum(time<=t1);
ul=sum(time<=t2);



figure
plot(time(lol:ul),data(lol:ul))
xlabel('time(s)')
ylabel('Voltage(V)')
title('dE/dt data')
%plot(data)