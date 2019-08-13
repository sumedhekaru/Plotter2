function [t,y] = dedtExtract(fn,t1,t2,ch_number)
% clc
% Triggers
fn = sprintf('%shh%i',fn(1:end-3),ch_number);

trigs = load(fn,'-mat');

% Find triggers between t1 and t2
lol=sum(trigs.trigs<t1)+1 ;    % Index of the lower matrix element
ul =sum(trigs.trigs<t2);

trigs = trigs.trigs(lol:ul);

leng = length(trigs);

if leng > 1;
    
    msg = sprintf('Trigger at %f s\n',trigs);
    errordlg(msg,'Large Amount of data');

    uiwait    
    t=NaN;
    y=NaN;
    return
    
elseif ul == 0
    t=NaN;
    y=NaN;
    return
    
else
    if lol > 1
        lol= lol -1;
    end
    
    if ul < leng
        ul = ul +1;
    end
end


str = sprintf('%s*.fa%i', fn(1:end-21),ch_number);
list = ls(str);


t=[];
y=[];
% Find file names for data
% Directory name
for i=lol:ul

    str = sprintf('%s%s', fn(1:end-21),list(i,:));
    d = load(str,'-mat');
    
    timenow=d.timenow;
    data=d.data;
    %size(data)
    time = timenow(4)*3600+timenow(5)*60+timenow(6)-(d.PPSt0-d.t0)*1e-6+d.act2.ActualStart/d.freq1:...
        1/d.freq1: ...
        timenow(4)*3600+timenow(5)*60+timenow(6)-(d.PPSt0-d.t0)*1e-6+(d.act2.ActualStart+d.act2.ActualLength-1)/d.freq1;
    
    y=[y data NaN];
    t=[t time NaN];
    
end
    


lol=sum(time<t1)+1;
ul=sum(time<t2);

if lol < ul
    t=time(lol:ul);
    y=data(lol:ul);
else
    % Means do data in the time range
    y = NaN;
    t = NaN;
end

% figure
% plot(time(lol:ul),data(lol:ul))
% xlabel('time(s)')
% ylabel('Voltage(V)')
% title('dE/dt data')
% %plot(data)