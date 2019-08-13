function LDAR_Hist1(fn,linet_fn,t1,t2,rc,x0,y0,z0,sen_set)

%fn='F:/data/ldar2/2010/07/11/ldar2_20101920000.txt';            % File name for idar   
wbh=waitbar(0.00,'Starting..','Name','Plotter Busy');


%t1=0;
%t2=456890066666;
%rc=50000;

if sen_set.ldarOn==1 && sen_set.linetOn==1
    waitbar(0.10,wbh,'Loading LDAR data','Name','Plotter Busy');
    [CG,CAL,DLS]=ldarExtract2(fn,t1,t2,rc,x0,y0,z0,0);
    
    waitbar(0.2,wbh,'Loading LINET data','Name','Plotter Busy');
    
       
    LINET=linetExtract(linet_fn,t1,t2,rc,x0,y0,z0);
    waitbar(0.3,wbh,'Finish loading LINET data','Name','Plotter Busy');
    
    %t=[CG(:,10);DLS(:,10);LINET(:,3)];         % time
    %x=[CG(:,6);DLS(:,6);LINET(:,6)]/1000;      % East
    %y=[CG(:,7);DLS(:,7);LINET(:,7)]/1000;      % North
    data=[CG(:,8);DLS(:,8);LINET(:,8)]/1000;      % Altitude
    %r=[CG(:,11);DLS(:,11);LINET(:,2)]/1000;    % Distance from a sensor
    
    % Plot title
    tit='LDAR & LINET';
    
    waitbar(0.35,wbh,'Finish loading LINET data','Name','Plotter Busy');
    
elseif sen_set.ldarOn==1
    waitbar(0.1,wbh,'Loading LDAR data','Name','Plotter Busy');
    
    [CG,CAL,DLS]=ldarExtract2(fn,t1,t2,rc,x0,y0,z0,0);
    
    waitbar(0.2,wbh,'Finish loading LDAR data','Name','Plotter Busy');
    
    %t=[CG(:,10);DLS(:,10)];         % time
    %x=[CG(:,6);DLS(:,6)]/1000;      % East
    %y=[CG(:,7);DLS(:,7)]/1000;      % North
    data=[CG(:,8);DLS(:,8)]/1000;      % Altitude
    %r=[CG(:,11);DLS(:,11)]/1000;    % Distance from a sensor
    
    %LINET=[NaN NaN NaN NaN NaN NaN NaN NaN];
     % Plot title
    tit='LDAR';
    
    waitbar(0.3,wbh,'Finish loading LDAR data','Name','Plotter Busy');
    
elseif sen_set.linetOn==1
    waitbar(0.1,wbh,'Loading LINET data','Name','Plotter Busy');
    LINET=linetExtract(linet_fn,t1,t2,rc,x0,y0,z0);
    waitbar(0.2,wbh,'Finish loading LINET data','Name','Plotter Busy');
    
    %t=LINET(:,3);         % time
    %x=LINET(:,6)/1000;      % East
    %y=LINET(:,7)/1000;      % North
    data=LINET(:,8)/1000;      % Altitude  
    %r=LINET(:,2)/1000;    % Distance from a sensor
    
    %CG=[NaN NaN NaN NaN NaN NaN NaN NaN];
    %DLS=[NaN NaN NaN NaN NaN NaN NaN NaN];
    
    % Plot title
    tit='LINET';
    
    waitbar(0.30,wbh,'Finish loading LINET data','Name','Plotter Busy');
   
end

%% Exit if there is no data in the given time range
if (length(data)-sum(isnan(data)))==0
    delete(wbh)
    errordlg('No data in the given time range.',...
        'Histogram Error','modal')
    return
end




waitbar(0.4,wbh,'Finish loading  data..','Name','Plotter Busy');


n= histc(data,-100:1000:20000);

fg=figure;
set(fg,'visible','off')

waitbar(0.6,wbh,'Creating the Hitogram..','Name','Plotter Busy');

barh(n)
tit=sprintf('%s Histogram    \n%s-%s-%s  UT:%s - %s',...
    tit,fn(end-31:end-28),fn(end-26:end-25),fn(end-23:end-22),sec2hhmmss(t1),sec2hhmmss(t2));
title(tit)

waitbar(0.8,wbh,'Creating the Hitogram..','Name','Plotter Busy');

xlabel('Number of LDAR points')
ylabel('Altitude (km)')
box on
grid on

waitbar(1.0,wbh,'Creating the Hitogram..','Name','Plotter Busy');
delete(wbh)

set(fg,'visible','on')