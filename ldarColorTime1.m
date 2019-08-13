%% IDAR radiation sources color coded by time by Sumedhe Karunarathne
% Generate color coded data according to time for IDAR data in given file
% name 'fn' in given time range 't1' & 't2' and inside the given radius
% 'rc'. Have ability to draw graphs seperately or as sub plots.

function ldarColorTime1(fnn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,tl,tu,...
    rcc,x0,y0,z0,seperatee,graphsOnn,sen_set,g)

%clc

ldar_fn=fnn;                                % File name for idar 

t1=tl;                                % Lower time bound
t2=tu;                                % Upper time bound

rc=rcc;                               % Bound radus 

mz = 3; % marker size

seperate=seperatee;                             
% if seperate =1, this will draw 4 seperate graph windaws instead of 
% subplot plotting
                                        
graphsOn=graphsOnn;
% If seperate ==1, then only seleted graphs will plot by graphOn
% respectively
%       (1) altitude vs time
%       (2) altitude vs east
%       (3) north vs east
%       (4) north vs altitude
%       (5) altitude vs north
% example: to turn on (1) and (3) graphs only
% graphOn=[1,0,1,0,0]


absant_fn={};

fn=ldar_fn;

%sen_set

if sen_set.ldarOn==1 || sen_set.cglssOn == 1
    if exist(ldar_fn,'file')==0
        absant_fn=[absant_fn fn];
        ldar_fn='';
    end
else
    ldar_fn='';
end

if sen_set.nldnOn==1
    if exist(nldn_fn,'file')==0
        absant_fn=[absant_fn nldn_fn];
        nldn_fn='';
    end
else
     nldn_fn='';
end

if sen_set.linetOn==1
    if exist(linet_fn,'file')==0
        absant_fn=[absant_fn linet_fn];
        linet_fn='';
    end
else
     linet_fn='';
end

if sen_set.pbfaOn==1
    if exist(pbfa_fn,'file')==0
        absant_fn=[absant_fn pbfa_fn];
        pbfa_fn='';
    end
else
    pbfa_fn='';
end

if sen_set.pbfaOOn==1
    if exist(pbfaO_fn,'file')==0
        absant_fn=[absant_fn pbfaO_fn];
        pbfaO_fn='';
    end
else
    pbfaO_fn='';
end


if sen_set.nldn2On==1
    if exist(nldn2_fn,'file')==0
        absant_fn=[absant_fn nldn2_fn];
        nldn2_fn='';
    end
else
    nldn2_fn='';
end

if isempty(absant_fn)==0
    answer=questdlg(absant_fn','Filese Missing','OK','Exit','OK');
    if strcmp(answer,'Exit')
        return
    end
end


wbh=waitbar(0.0,'Loading LDAR data','Name','Plotter Busy');

t=[];
x=[];
y=[];
z=[];
r=[];
tit='';



if isempty(ldar_fn)==0
        
     waitbar(0.05,wbh,'Loading LDAR data','Name','Plotter Busy');
    
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,t1,t2,rc,x0,y0,z0,0);
    
    waitbar(0.15,wbh,'Finish loading LDAR data','Name','Plotter Busy');
    
    if sen_set.cglssOn
         t=CG(:,10);         % time
         x=CG(:,6)/1000;      % East
         y=CG(:,7)/1000;      % North
         z=CG(:,8)/1000;      % Altitude
         r=CG(:,11)/1000;    % Distance from a sensor
    end
    
    if sen_set.ldarOn
        t=[t;DLS(:,10)];         % time
        x=[x;DLS(:,6)/1000];      % East
        y=[y;DLS(:,7)/1000];      % North
        z=[z;DLS(:,8)/1000];      % Altitude
        r=[r;DLS(:,11)/1000];    % Distance from a sensor
    end
    
    % Plot title
    tit=[tit ' LDAR'];
    
    waitbar(0.08,wbh,'Loading LDAR data','Name','Plotter Busy');
else
    CG=[NaN NaN NaN NaN NaN NaN NaN NaN];
    DLS=[NaN NaN NaN NaN NaN NaN NaN NaN];
    waitbar(0.08,wbh,'Loading LDAR data','Name','Plotter Busy');
end


if isempty(linet_fn)==0
    waitbar(0.10,wbh,'Loading LINET data','Name','Plotter Busy');
    LINET=linetExtract2(linet_fn,t1,t2,rc,x0,y0,z0);
    waitbar(0.12,wbh,'Finish loading LINET data','Name','Plotter Busy');
    
    t=[t;LINET(:,3)];           % time
    x=[x;LINET(:,6)/1000];      % East
    y=[y;LINET(:,7)/1000];      % North
    z=[z;LINET(:,8)/1000];      % Altitude  
    r=[r;LINET(:,2)/1000];      % Distance from a sensor
    
    
    % Plot title
    tit=[tit ' LINET'];
    
    waitbar(0.15,wbh,'Loading LINET data','Name','Plotter Busy');
else
        LINET=[NaN NaN NaN NaN NaN NaN NaN NaN];
        waitbar(0.15,wbh,'Loading LINET data','Name','Plotter Busy');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if isempty(nldn_fn)==0
%     waitbar(0.151,wbh,'Loading LINET2 data','Name','Plotter Busy');
%     NLDN=linetExtract(nldn_fn,t1,t2,rc,x0,y0,z0);
%     waitbar(0.155,wbh,'Finish loading LINET2 data','Name','Plotter Busy');
%     
%     t=[t;NLDN(:,3)];           % time
%     x=[x;NLDN(:,6)/1000];      % East
%     y=[y;NLDN(:,7)/1000];      % North
%     z=[z;NLDN(:,8)/1000];      % Altitude  
%     r=[r;NLDN(:,2)/1000];      % Distance from a sensor
%     
%     
%     % Plot title
%     tit=[tit ' LINET2'];
%     
%     waitbar(0.15,wbh,'Loading LINET data','Name','Plotter Busy');
% else
%         LINET=[NaN NaN NaN NaN NaN NaN NaN NaN];
%         waitbar(0.15,wbh,'Loading LINET data','Name','Plotter Busy');
% end

if isempty(nldn_fn)==0
    waitbar(0.16,wbh,'Loading PBFA-A data','Name','Plotter Busy');
    PBFAA=pbfaExtract(nldn_fn,t1,t2,rc,x0,y0,z0);
    waitbar(0.18,wbh,'Loading PBFA data','Name','Plotter Busy');
    
    t=[t;PBFAA(:,6)];           % time
    x=[x;PBFAA(:,3)/1000];      % East
    y=[y;PBFAA(:,4)/1000];      % North
    z=[z;PBFAA(:,5)/1000];      % Altitude  
    r=[r;PBFAA(:,1)/1000];      % Distance from a sensor
    
   
    % Plot title
    tit=[tit ' PBFA-A'];
    
    waitbar(0.20,wbh,'Loading PBFA data','Name','Plotter Busy');
else
    PBFAA=[NaN NaN NaN NaN NaN NaN];
    waitbar(0.20,wbh,'Loading PBFA data','Name','Plotter Busy');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(pbfa_fn)==0
    waitbar(0.16,wbh,'Loading PBFA data','Name','Plotter Busy');
    PBFA=pbfaExtract(pbfa_fn,t1,t2,rc,x0,y0,z0);
    waitbar(0.18,wbh,'Loading PBFA data','Name','Plotter Busy');
    
    t=[t;PBFA(:,6)];           % time
    x=[x;PBFA(:,3)/1000];      % East
    y=[y;PBFA(:,4)/1000];      % North
    z=[z;PBFA(:,5)/1000];      % Altitude  
    r=[r;PBFA(:,1)/1000];      % Distance from a sensor
    
   
    % Plot title
    tit=[tit ' PBFA'];
    
    waitbar(0.20,wbh,'Loading PBFA data','Name','Plotter Busy');
else
    PBFA=[NaN NaN NaN NaN NaN NaN];
    waitbar(0.20,wbh,'Loading PBFA data','Name','Plotter Busy');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(pbfaO_fn)==0
    waitbar(0.16,wbh,'Loading PBFA data','Name','Plotter Busy');
    PBFAO=pbfaExtract(pbfaO_fn,t1,t2,rc,x0,y0,z0);
    waitbar(0.18,wbh,'Loading PBFA data','Name','Plotter Busy');
    
    t=[t;PBFAO(:,6)];           % time
    x=[x;PBFAO(:,3)/1000];      % East
    y=[y;PBFAO(:,4)/1000];      % North
    z=[z;PBFAO(:,5)/1000];      % Altitude  
    r=[r;PBFAO(:,1)/1000];      % Distance from a sensor
    
   
    % Plot title
    tit=[tit ' PBFA-O'];
    
    waitbar(0.20,wbh,'Loading PBFA data','Name','Plotter Busy');
else
    PBFAO=[NaN NaN NaN NaN NaN NaN];
    waitbar(0.20,wbh,'Loading PBFA data','Name','Plotter Busy');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(nldn2_fn)==0
    waitbar(0.16,wbh,'Loading PBFA data','Name','Plotter Busy');
    [NLDN2c NLDN2g] = nldnExtract(nldn2_fn,0,0,t1,t2,...
        rc,x0,y0,z0,0);
    waitbar(0.18,wbh,'Loading PBFA data','Name','Plotter Busy');
    
   
    t=[t;NLDN2c(:,8);NLDN2g(:,8)];           % time
    x=[x;NLDN2c(:,2)/1000;NLDN2g(:,2)/1000];      % East
    y=[y;NLDN2c(:,3)/1000;NLDN2g(:,3)/1000];      % North
    z=[z;NLDN2c(:,4)/1000;NLDN2g(:,4)/1000];      % Altitude  
    r=[r;NLDN2c(:,9)/1000;NLDN2g(:,9)/1000];      % Distance from a sensor
    
   
    % Plot title
    tit=[tit ' NLDN2'];
    
    waitbar(0.20,wbh,'Loading NLDN2 data','Name','Plotter Busy');
else
    NLDN2c=NaN(1,9); 
    NLDN2g=NaN(1,9);
    waitbar(0.20,wbh,'Loading NLDN2 data','Name','Plotter Busy');
end

% Exit if there is no data in the given time range
if (length(t)-sum(isnan(t)))==0
    delete(wbh)
    errordlg('No data in the given time range.',...
        'LDAR 3D Error','modal')
    return
end



   % Plot title
    tit=sprintf('%s Color Coded by Time    \n%s-%s-%s     UT: %s - %s',...
        tit,fn(end-31:end-28),fn(end-26:end-25),fn(end-23:end-22),sec2hhmmss(t1),sec2hhmmss(t2));





%% Exit if there is no data in the given time range
if (length(t)-sum(isnan(t)))==0
    delete(wbh)
    errordlg('No data in the given time range.',...
        'LDAR CT Error','modal')
    return
end





% Loading additional Land Marks given in 'add_points.txt'
if sen_set.lmOn==1
    [lmID,lmx,lmy,lmz]=textread('add_points.txt','%s %f %f %f','headerlines',8);
    % Convert to Km
    lmx=lmx/1000;
    lmy=lmy/1000;
    lmz=lmz/1000;
else
    lmx=[];
    lmy=[];
    lmz=[];
    lmID='';
end

waitbar(0.4,wbh,'Plotting the graph','Name','Plotter Busy');




fg=figure;
set(fg,'visible','off')


map=colormap;
%miv=min(t);
%mav=max(t);

miv = g.t1;
mav = g.t2;

clrstep = (mav-miv)/size(map,1);

marker='o';

%% Plotting z-t graph
if seperate==1 && graphsOn(1)==1
    
    % Legend control
    lg_ldar = 0;
    lg_cglss = 0;
    lg_pbfa = 0;
    lg_pbfaA = 0;
    lg_pbfaO = 0;
    lg_nldn2c = 0;
    lg_nldn2g = 0;
    lg_linet = 0;
    lg = {};
    
    waitbar(0.5,wbh,'Plotting the graph','Name','Plotter Busy');
    %figure
    hold on
    for nc=1:size(map,1)
        iv = find(t>=miv+(nc-1)*clrstep & t<=miv+nc*clrstep);
        
        for i1 = 1:length(iv)
            %Is it a LDAR2
            if  sum(x(iv(i1)) ==  DLS(:,6)/1000) >= 1 ...
                    && sum(y(iv(i1)) == DLS(:,7)/1000) >= 1
                %disp('DLS')
                hLine = plot(t(iv(i1)),z(iv(i1)),'o','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_ldar
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'LDAR2'];
                    lg_ldar = 1;
                end
                
                % Is it a CGLSS
            elseif sum(x(iv(i1)) ==  CG(:,6)/1000) >= 1 ...
                    && sum(y(iv(i1)) == CG(:,7)/1000) >= 1
                %disp('CG')
                hLine  = plot(t(iv(i1)),z(iv(i1)),'s','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_cglss
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'CGLSS'];
                    lg_cglss = 1;
                end
                
                % Is it a PBFA
            elseif sum(x(iv(i1)) ==  PBFA(:,3)/1000) >= 1 ...
                    && sum(y(iv(i1)) == PBFA(:,4)/1000) >= 1
                hLine = plot(t(iv(i1)),z(iv(i1)),'d','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfa
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'PBFA'];
                    lg_pbfa = 1;
                end
                
                %disp('PBFA')
                %Is it PBFA Auto?
            elseif sum(x(iv(i1)) ==  PBFAA(:,3)/1000) >= 1 ...
                    && sum(y(iv(i1)) == PBFAA(:,4)/1000) >= 1
                
                hLine = plot(t(iv(i1)),z(iv(i1)),'p','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfaA
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'PBFA-A'];
                    lg_pbfaA = 1;
                end
                
                %disp('PBFA')
                %Is it PBFA Old?
            elseif sum(x(iv(i1)) ==  PBFAO(:,3)/1000) >= 1 ...
                    && sum(y(iv(i1)) == PBFAO(:,4)/1000) >= 1
                hLine = plot(t(iv(i1)),z(iv(i1)),'v','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfaO
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'PBFA-O'];
                    lg_pbfaO = 1;
                end
                
                %disp('PBFA')   
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  NLDN2c(:,2)/1000) >= 1 ...
                    && sum(y(iv(i1)) == NLDN2c(:,3)/1000) >= 1
                hLine = plot(t(iv(i1)),z(iv(i1)),'<','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_nldn2c
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'NLDN2-C'];
                    lg_nldn2c = 1;
                end
                
                %disp('PBFA')
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  NLDN2g(:,2)/1000) >= 1 ...
                    && sum(y(iv(i1)) == NLDN2g(:,3)/1000) >= 1
                hLine = plot(t(iv(i1)),z(iv(i1)),'>','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_nldn2g
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'NLDN2-G'];
                    lg_nldn2g = 1;
                end
                
                %disp('PBFA'
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  LINET(:,6)/1000) >= 1 ...
                    && sum(y(iv(i1)) == LINET(:,7)/1000) >= 1
                hLine = plot(t(iv(i1)),z(iv(i1)),'>','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_linet
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'LINET'];
                    lg_linet = 1;
                end
            else
                disp('Non')
                plot(t(iv(i1)),z(iv(i1)),'p','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1)
            end
            
            
        end
        %plot(t(iv),z(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
    legend(lg)
    hold off
    
    waitbar(0.6,wbh,'Plotting the graph','Name','Plotter Busy');
    
    title(tit)
    ylabel('Altitude (km)')
    xlabel('Time (s)')
    box on

    % Re-format the colorbar
    h=colorbar;
    
    waitbar(0.7,wbh,'Plotting the graph','Name','Plotter Busy');
    
    %set(h,'ylim',[1 length(map)]);
    yal=linspace(1,length(map),10);
    set(h,'ytick',yal);
    % Create the yticklabels
    ytl=linspace(miv,mav,10);
    s=char(10,4);
    
    waitbar(0.8,wbh,'Plotting the graph','Name','Plotter Busy');
    
    for i=1:10
        if min(abs(ytl)) >= 0.001
            B=sprintf('%-4.3f',ytl(i));
        else
            B=sprintf('%-3.1E',ytl(i));
        end
        s(i,1:length(B))=B;
    end
    
    waitbar(0.9,wbh,'Plotting the graph','Name','Plotter Busy');
    
    set(h,'yticklabel',s);
    
    waitbar(1.0,wbh,'Done','Name','Plotter Busy');
    
    delete(wbh)
    
    set(fg,'visible','on')
        
   
elseif seperate ~=1
    waitbar(0.6,wbh,'Plotting the graph','Name','Plotter Busy');
    subplot(3,2,1:2)

    hold on
    for nc=1:size(map,1)
        iv = find(t>=miv+(nc-1)*clrstep & t<=miv+nc*clrstep) ;
        plot(t(iv),z(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
    hold off

    title(tit)
    ylabel('Altitude (km)')
    xlabel('Time (s)')
    box on

    % Re-format the colorbar
    h=colorbar;

    %set(h,'ylim',[1 length(map)]);
    yal=linspace(1,length(map),10);
    set(h,'ytick',yal);
    % Create the yticklabels
    ytl=linspace(miv,mav,10);
    s=char(10,4);
    for i=1:10
        if min(abs(ytl)) >= 0.001
            B=sprintf('%-4.3f',ytl(i));
        else
            B=sprintf('%-3.1E',ytl(i));
        end
        s(i,1:length(B))=B;
    end
    
    set(h,'yticklabel',s);


end


%% Plotting z-x graph

if seperate==1 && graphsOn(2)==1
   
    waitbar(0.5,wbh,'Plotting the graph','Name','Plotter Busy');
    
    %figure
    xlabel('East (km)')
    ylabel('Altitude (km)')
    title(tit)
    box on
    
    % Legend control
    lg_ldar = 0;
    lg_cglss = 0;
    lg_pbfa = 0;
    lg_pbfaA = 0;
    lg_pbfaO = 0;
    lg_nldn2c = 0;
    lg_nldn2g = 0;
    lg_linet = 0;
    lg = {};
    
    waitbar(0.6,wbh,'Plotting the graph','Name','Plotter Busy');
    
    hold on
    for nc=1:size(map,1)
        iv = find(t>=miv+(nc-1)*clrstep & t<=miv+nc*clrstep) ;
        
                for i1 = 1:length(iv)
            %Is it a LDAR2
            if  sum(x(iv(i1)) ==  DLS(:,6)/1000) >= 1 ...
                    && sum(y(iv(i1)) == DLS(:,7)/1000) >= 1
                %disp('DLS')
                hLine = plot(x(iv(i1)),z(iv(i1)),'o','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_ldar
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'LDAR2'];
                    lg_ldar = 1;
                end
                
                % Is it a CGLSS
            elseif sum(x(iv(i1)) ==  CG(:,6)/1000) >= 1 ...
                    && sum(y(iv(i1)) == CG(:,7)/1000) >= 1
                %disp('CG')
                hLine  = plot(x(iv(i1)),z(iv(i1)),'s','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_cglss
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'CGLSS'];
                    lg_cglss = 1;
                end
                
                % Is it a PBFA
            elseif sum(x(iv(i1)) ==  PBFA(:,3)/1000) >= 1 ...
                    && sum(y(iv(i1)) == PBFA(:,4)/1000) >= 1
                hLine = plot(x(iv(i1)),z(iv(i1)),'d','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfa
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'PBFA'];
                    lg_pbfa = 1;
                end
                
                %disp('PBFA')
                %is it PBFA-A?
            elseif sum(x(iv(i1)) ==  PBFAA(:,3)/1000) >= 1 ...
                    && sum(y(iv(i1)) == PBFAA(:,4)/1000) >= 1
                hLine = plot(x(iv(i1)),z(iv(i1)),'p','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfaA
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'PBFA-A'];
                    lg_pbfaA = 1;
                end
                
                %disp('PBFA')
                                %is it PBFA-A?
            elseif sum(x(iv(i1)) ==  PBFAO(:,3)/1000) >= 1 ...
                    && sum(y(iv(i1)) == PBFAO(:,4)/1000) >= 1
                hLine = plot(x(iv(i1)),z(iv(i1)),'v','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfaO
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'PBFA-O'];
                    lg_pbfaO = 1;
                end
                
                %disp('PBFA')
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  NLDN2c(:,2)/1000) >= 1 ...
                    && sum(y(iv(i1)) == NLDN2c(:,3)/1000) >= 1
                hLine = plot(x(iv(i1)),z(iv(i1)),'<','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_nldn2c
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'NLDN2-C'];
                    lg_nldn2c = 1;
                end
                
                %disp('PBFA')
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  NLDN2g(:,2)/1000) >= 1 ...
                    && sum(y(iv(i1)) == NLDN2g(:,3)/1000) >= 1
                hLine = plot(x(iv(i1)),z(iv(i1)),'>','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_nldn2g
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'NLDN2-G'];
                    lg_nldn2g = 1;
                end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  LINET(:,6)/1000) >= 1 ...
                    && sum(y(iv(i1)) == LINET(:,7)/1000) >= 1
                hLine = plot(x(iv(i1)),z(iv(i1)),'>','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_linet
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'LINET'];
                    lg_linet = 1;
                end
                
            else
                %disp('Non')
                plot(x(iv(i1)),z(iv(i1)),'p','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1)
            end
            
            
        end
        %plot(x(iv),z(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
    legend(lg)
    waitbar(0.7,wbh,'Plotting the graph','Name','Plotter Busy');
    
%     plot(sen_set.x/1000,sen_set.z/1000,'rp','markerfacecolor','r')
%     text(sen_set.x/1000,sen_set.z/1000,sen_set.sen_IDs)
    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot(sen_set.x(ns)/1000,sen_set.z(ns)/1000,'kp','markerfacecolor','r')
            text(sen_set.x(ns)/1000,sen_set.z(ns)/1000,sen_set.sen_IDs{ns})
        end
    end
    
    % Plotting Additional land marks
    plot(lmx,lmz,'bp','markerfacecolor','b')
    text(lmx,lmz,lmID)
    
    hold off
    
    waitbar(0.8,wbh,'Plotting the graph','Name','Plotter Busy');
    
    % Re-format the colorbar
    h=colorbar;

    %set(h,'ylim',[1 length(map)]);
    yal=linspace(1,length(map),10);
    set(h,'ytick',yal);
    % Create the yticklabels
    ytl=linspace(miv,mav,10);
    s=char(10,4);
    
    waitbar(0.9,wbh,'Plotting the graph','Name','Plotter Busy');
    
    for i=1:10
        if min(abs(ytl)) >= 0.001
            B=sprintf('%-4.3f',ytl(i));
        else
            B=sprintf('%-3.1E',ytl(i));
        end
        s(i,1:length(B))=B;
    end
    
    set(h,'yticklabel',s);
    
    waitbar(1.0,wbh,'Done','Name','Plotter Busy');
    delete(wbh)
    
    set(fg,'visible','on')
    
elseif seperate ~=1
    waitbar(0.7,wbh,'Plotting the graph','Name','Plotter Busy');
    subplot(3,2,3)

    hold on
    for nc=1:size(map,1)
        iv = find(t>=miv+(nc-1)*clrstep & t<=miv+nc*clrstep) ;
        plot(x(iv),z(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
    
%     plot(sen_set.x/1000,sen_set.z/1000,'rp','markerfacecolor','r')
%     text(sen_set.x/1000,sen_set.z/1000,sen_set.sen_IDs)

    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot(sen_set.x(ns)/1000,sen_set.z(ns)/1000,'kp','markerfacecolor','r')
            text(sen_set.x(ns)/1000,sen_set.z(ns)/1000,sen_set.sen_IDs{ns})
        end
    end

    % Plotting Additional land marks
    plot(lmx,lmz,'bp','markerfacecolor','b')
    text(lmx,lmz,lmID)
    
    hold off

    ylabel('Altitude (km)')
    box on

end

%subplot(3,2,4)

%% Plotting y-x graph

% Legend control
lg_ldar = 0;
lg_cglss = 0;
lg_pbfa = 0;
lg_pbfaA = 0;
lg_pbfaO = 0;
lg_nldn2c = 0;
lg_nldn2g = 0;
lg_linet = 0;
lg = {};

if seperate==1 && graphsOn(3)==1
    
    % Plot radar
    if sen_set.radarOn && ~isempty(sen_set.radarFn)
        radar_plot(sen_set,g)
    end    
    
    waitbar(0.5,wbh,'Plotting the graph','Name','Plotter Busy');
    %clc
    %figure
    hold on
    
    date = str2double(sprintf('%s%s%s',g.YYYY{:},g.MM{:},g.DD{:}));

    if sen_set.mapOn
        if date > 20110610 && date < 20110628
            xdd = [x; 0; 9.391];
            ydd = [y; 0; 3.089];
            xlim([min(xdd) max(xdd)]*1.3)
            ylim([min(ydd) max(ydd)]*1.3)
            rapidCity_map
        else
            xdd = [x; -60.1; 3.4];
            ydd = [y; -53.1; 49.3];
            xlim([min(xdd) max(xdd)]*1.3)
            ylim([min(ydd) max(ydd)]*1.3)
            florida_map
        end
    end
    
    for nc=1:size(map,1)
        iv = find(t>=miv+(nc-1)*clrstep & t<=miv+nc*clrstep);
        for i1 = 1:length(iv)
            %Is it a LDAR2
            if  sum(x(iv(i1)) ==  DLS(:,6)/1000) >= 1 ...                
                && sum(y(iv(i1)) == DLS(:,7)/1000) >= 1
                %disp('DLS')
                hLine = plot(x(iv(i1)),y(iv(i1)),'o','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_ldar
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off');                        
                else
                    lg = [lg 'LDAR2'];
                    lg_ldar = 1;
                end
                    
            % Is it a CGLSS    
            elseif sum(x(iv(i1)) ==  CG(:,6)/1000) >= 1 ...                
                && sum(y(iv(i1)) == CG(:,7)/1000) >= 1
                %disp('CG')
                hLine  = plot(x(iv(i1)),y(iv(i1)),'o','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_cglss
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                    
                else
                    lg = [lg 'CGLSS'];
                    lg_cglss = 1;
                end
                
            % Is it a PBFA
            elseif sum(x(iv(i1)) ==  PBFA(:,3)/1000) >= 1 ...                
                && sum(y(iv(i1)) == PBFA(:,4)/1000) >= 1
                hLine = plot(x(iv(i1)),y(iv(i1)),'o','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_pbfa
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                    
                else
                    lg = [lg 'PBFA'];
                    lg_pbfa = 1;
                end
                
                %disp('PBFA')
            % Is it a PBFA
            elseif sum(x(iv(i1)) ==  PBFAA(:,3)/1000) >= 1 ...                
                && sum(y(iv(i1)) == PBFAA(:,4)/1000) >= 1
                hLine = plot(x(iv(i1)),y(iv(i1)),'o','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_pbfaA
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                    
                else
                    lg = [lg 'PBFA-A'];
                    lg_pbfaA = 1;
                end
                
                %disp('PBFA')
            % Is it a PBFA Old?
            elseif sum(x(iv(i1)) ==  PBFAO(:,3)/1000) >= 1 ...                
                && sum(y(iv(i1)) == PBFAO(:,4)/1000) >= 1
                hLine = plot(x(iv(i1)),y(iv(i1)),'o','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_pbfaO
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                    
                else
                    lg = [lg 'PBFA-O'];
                    lg_pbfaO = 1;
                end
                
                %disp('PBFA')
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  NLDN2c(:,2)/1000) >= 1 ...
                    && sum(y(iv(i1)) == NLDN2c(:,3)/1000) >= 1
                hLine = plot(x(iv(i1)),y(iv(i1)),'*','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_nldn2c
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'NLDN2-C'];
                    lg_nldn2c = 1;
                end
                
                %disp('PBFA')
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  NLDN2g(:,2)/1000) >= 1 ...
                    && sum(y(iv(i1)) == NLDN2g(:,3)/1000) >= 1
                hLine = plot(x(iv(i1)),y(iv(i1)),'+','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_nldn2g
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'NLDN2-G'];
                    lg_nldn2g = 1;
                end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            elseif sum(x(iv(i1)) ==  LINET(:,6)/1000) >= 1 ...
                    && sum(y(iv(i1)) == LINET(:,7)/1000) >= 1
                hLine = plot(x(iv(i1)),y(iv(i1)),'+','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_linet
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'LINET'];
                    lg_linet = 1;
                end
                
            else
                %disp('Non')
                plot(x(iv(i1)),y(iv(i1)),'p','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1)
            end
                        
            
        end
        %plot(x(iv),y(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
 
    legend(lg)
    
    waitbar(0.6,wbh,'Plotting the graph','Name','Plotter Busy');
    
%     plot(sen_set.x/1000,sen_set.y/1000,'rp','markerfacecolor','r')
%     text(sen_set.x/1000,sen_set.y/1000,sen_set.sen_IDs) 


    
    if date > 20110610 && date < 20110628
        
        plot([0,9.391],[0,3.089],'kp','markerfacecolor','r','markersize',mz)
        text(0,0,'West')
        text(9.391,3.089,'East')
        
    else        
        for ns = 1:20
            if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
                plot(sen_set.x(ns)/1000,sen_set.y(ns)/1000,'kp','markerfacecolor','r','MarkerSize',mz)
                text(sen_set.x(ns)/1000,sen_set.y(ns)/1000,sen_set.sen_IDs{ns})
            end
        end
    end
    
    
    
    
    % Plotting Additional land marks
    plot(lmx,lmy,'bp','markerfacecolor','b')
    text(lmx,lmy,lmID)
    
    hold off
    
    waitbar(0.7,wbh,'Plotting the graph','Name','Plotter Busy');
    
    ylabel('North (km)')
    xlabel('East (km)')
    title(tit)
    box on
    
    
    % Re-format the colorbar
    h=colorbar;
    
    waitbar(0.8,wbh,'Plotting the graph','Name','Plotter Busy');
    
    % Format colorbar
    format_colorbar(h,miv,mav)
    
    delete(wbh)
    
   
    set(fg,'visible','on')
    
elseif seperate ~=1
    
    waitbar(0.8,wbh,'Plotting the graph','Name','Plotter Busy');
    
    subplot(3,2,5)
    hold on
    for nc=1:size(map,1)
        iv = find(t>=miv+(nc-1)*clrstep & t<=miv+nc*clrstep) ;
        plot(x(iv),y(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
    
%     plot(sen_set.x/1000,sen_set.y/1000,'rp','markerfacecolor','r')
%     text(sen_set.x/1000,sen_set.y/1000,sen_set.sen_IDs)

    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot(sen_set.x(ns)/1000,sen_set.y(ns)/1000,'kp','markerfacecolor','r')
            text(sen_set.x(ns)/1000,sen_set.y(ns)/1000,sen_set.sen_IDs{ns})
        end
    end
    
    % Plotting Additional land marks
    plot(lmx,lmy,'bp','markerfacecolor','b')
    text(lmx,lmy,lmID)
    
    hold off

    ylabel('North (km)')
    xlabel('East (km)')
    box on
end



%% plotting y-z graph

if seperate==1 && graphsOn(4)==1
    
    waitbar(0.5,wbh,'Plotting the graph','Name','Plotter Busy');
    
    % Legend control
    lg_ldar = 0;
    lg_cglss = 0;
    lg_pbfa = 0;
    lg_pbfaA = 0;
    lg_pbfaO = 0;
    lg_nldn2c = 0;
    lg_nldn2g = 0;
    lg_linet = 0;
    lg = {};
    
    %figure
    
    hold on
    for nc=1:size(map,1)
        iv = find(t>=miv+(nc-1)*clrstep & t<=miv+nc*clrstep) ;
        
                for i1 = 1:length(iv)
            %Is it a LDAR2
            if  sum(x(iv(i1)) ==  DLS(:,6)/1000) >= 1 ...                
                && sum(y(iv(i1)) == DLS(:,7)/1000) >= 1
                %disp('DLS')
                hLine = plot(z(iv(i1)),y(iv(i1)),'o','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_ldar
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off');                        
                else
                    lg = [lg 'LDAR2'];
                    lg_ldar = 1;
                end
                    
            % Is it a CGLSS    
            elseif sum(x(iv(i1)) ==  CG(:,6)/1000) >= 1 ...                
                && sum(y(iv(i1)) == CG(:,7)/1000) >= 1
                %disp('CG')
                hLine  = plot(z(iv(i1)),y(iv(i1)),'s','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_cglss
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                    
                else
                    lg = [lg 'CGLSS'];
                    lg_cglss = 1;
                end
                
            % Is it a PBFA
            elseif sum(x(iv(i1)) ==  PBFA(:,3)/1000) >= 1 ...                
                && sum(y(iv(i1)) == PBFA(:,4)/1000) >= 1
                hLine = plot(z(iv(i1)),y(iv(i1)),'d','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfa
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                    
                else
                    lg = [lg 'PBFA'];
                    lg_pbfa = 1;
                end
                
                %disp('PBFA')
            elseif sum(x(iv(i1)) ==  PBFAA(:,3)/1000) >= 1 ...                
                && sum(y(iv(i1)) == PBFAA(:,4)/1000) >= 1
                hLine = plot(z(iv(i1)),y(iv(i1)),'p','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfaA
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                    
                else
                    lg = [lg 'PBFA-A'];
                    lg_pbfaA = 1;
                end
                
                %disp('PBFA')
            elseif sum(x(iv(i1)) ==  PBFAO(:,3)/1000) >= 1 ...                
                && sum(y(iv(i1)) == PBFAO(:,4)/1000) >= 1
                hLine = plot(z(iv(i1)),y(iv(i1)),'v','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfaO
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                    
                else
                    lg = [lg 'PBFA-O'];
                    lg_pbfaO = 1;
                end
                
                %disp('PBFA') 
                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  NLDN2c(:,2)/1000) >= 1 ...
                    && sum(y(iv(i1)) == NLDN2c(:,3)/1000) >= 1
                hLine = plot(z(iv(i1)),y(iv(i1)),'<','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_nldn2c
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'NLDN2-C'];
                    lg_nldn2c = 1;
                end
                
                %disp('PBFA')
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  NLDN2g(:,2)/1000) >= 1 ...
                    && sum(y(iv(i1)) == NLDN2g(:,3)/1000) >= 1
                hLine = plot(z(iv(i1)),y(iv(i1)),'>','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_nldn2g
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'NLDN2-G'];
                    lg_nldn2g = 1;
                end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            elseif sum(x(iv(i1)) ==  LINET(:,6)/1000) >= 1 ...
                    && sum(y(iv(i1)) == LINET(:,7)/1000) >= 1
                hLine = plot(z(iv(i1)),y(iv(i1)),'>','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_linet
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'LINET'];
                    lg_linet = 1;
                end
            else
                %disp('Non')
                plot(z(iv(i1)),y(iv(i1)),'p','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1)
            end
                        
            
        end
        %plot(z(iv),y(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
    
    legend(lg)
    
    waitbar(0.6,wbh,'Plotting the graph','Name','Plotter Busy');
    
%     plot(sen_set.z/1000,sen_set.y/1000,'rp','markerfacecolor','r')
%     text(sen_set.z/1000,sen_set.y/1000,sen_set.sen_IDs)

    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot(sen_set.z(ns)/1000,sen_set.y(ns)/1000,'kp','markerfacecolor','r')
            text(sen_set.z(ns)/1000,sen_set.y(ns)/1000,sen_set.sen_IDs{ns})
        end
    end
    
    % Plotting Additional land marks
    plot(lmz,lmy,'bp','markerfacecolor','b')
    text(lmz,lmy,lmID)
    
    waitbar(0.7,wbh,'Plotting the graph','Name','Plotter Busy');
    
    hold off
    
    ylabel('North (km)')
    xlabel('Altitude (km)')
    title(tit)
    box on
    
    
    % Re-format the colorbar
    h=colorbar;

    waitbar(0.8,wbh,'Plotting the graph','Name','Plotter Busy');
    
    %set(h,'ylim',[1 length(map)]);
    yal=linspace(1,length(map),10);
    set(h,'ytick',yal);
    % Create the yticklabels
    ytl=linspace(miv,mav,10);
    s=char(10,4);
    
    waitbar(0.9,wbh,'Plotting the graph','Name','Plotter Busy');
    
    for i=1:10
        if min(abs(ytl)) >= 0.001
            B=sprintf('%-4.3f',ytl(i));
        else
            B=sprintf('%-3.1E',ytl(i));
        end
        s(i,1:length(B))=B;
    end
    
    set(h,'yticklabel',s);
    
    waitbar(1.0,wbh,'Done','Name','Plotter Busy');
    
    delete(wbh)
    
    set(fg,'visible','on')
    
elseif seperate ~=1
    waitbar(0.85,wbh,'Plotting the graph','Name','Plotter Busy');
    
    subplot(3,2,6)
    hold on
    for nc=1:size(map,1)
        iv = find(t>=miv+(nc-1)*clrstep & t<=miv+nc*clrstep) ;
        plot(z(iv),y(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
    
%     plot(sen_set.z/1000,sen_set.y/1000,'rp','markerfacecolor','r')
%     text(sen_set.z/1000,sen_set.y/1000,sen_set.sen_IDs)
    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot(sen_set.z(ns)/1000,sen_set.y(ns)/1000,'kp','markerfacecolor','r')
            text(sen_set.z(ns)/1000,sen_set.y(ns)/1000,sen_set.sen_IDs{ns})
        end
    end

    
    % Plotting Additional land marks
    plot(lmz,lmy,'bp','markerfacecolor','b')
    text(lmz,lmy,lmID)
    
    hold off

    xlabel('Altitude (km)')
    box on
    

    
end



%% plotting z-y graph

if seperate==1 && graphsOn(5)==1
    %figure
    
    waitbar(0.5,wbh,'Plotting the graph','Name','Plotter Busy');
    
    % Legend control
    lg_ldar = 0;
    lg_cglss = 0;
    lg_pbfa = 0;
    lg_pbfaA = 0;
    lg_pbfaO = 0;
    lg_nldn2c =0;
    lgn_nldn2g = 0;
    lg_linet = 0;
    lg = {};
    
    hold on
    for nc=1:size(map,1)
        iv = find(t>=miv+(nc-1)*clrstep & t<=miv+nc*clrstep) ;
                for i1 = 1:length(iv)
            %Is it a LDAR2
            if  sum(x(iv(i1)) ==  DLS(:,6)/1000) >= 1 ...                
                && sum(y(iv(i1)) == DLS(:,7)/1000) >= 1
                %disp('DLS')
                hLine = plot(y(iv(i1)),z(iv(i1)),'o','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_ldar
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off');                        
                else
                    lg = [lg 'LDAR2'];
                    lg_ldar = 1;
                end
                    
            % Is it a CGLSS    
            elseif sum(x(iv(i1)) ==  CG(:,6)/1000) >= 1 ...                
                && sum(y(iv(i1)) == CG(:,7)/1000) >= 1
                %disp('CG')
                hLine  = plot(y(iv(i1)),z(iv(i1)),'s','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
                if lg_cglss
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                    
                else
                    lg = [lg 'CGLSS'];
                    lg_cglss = 1;
                end
                
            % Is it a PBFA
            elseif sum(x(iv(i1)) ==  PBFA(:,3)/1000) >= 1 ...                
                && sum(y(iv(i1)) == PBFA(:,4)/1000) >= 1
                hLine = plot(y(iv(i1)),z(iv(i1)),'d','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfa
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                    
                else
                    lg = [lg 'PBFA'];
                    lg_pbfa = 1;
                end
                
                %disp('PBFA')
            % Is it a PBFA-A?
            elseif sum(x(iv(i1)) ==  PBFAA(:,3)/1000) >= 1 ...                
                && sum(y(iv(i1)) == PBFAA(:,4)/1000) >= 1
                hLine = plot(y(iv(i1)),z(iv(i1)),'p','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfaA
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                    
                else
                    lg = [lg 'PBFA-A'];
                    lg_pbfaA = 1;
                end
                
                %disp('PBFA')  
                
                            % Is it a PBFA-A?
            elseif sum(x(iv(i1)) ==  PBFAO(:,3)/1000) >= 1 ...                
                && sum(y(iv(i1)) == PBFAO(:,4)/1000) >= 1
                hLine = plot(y(iv(i1)),z(iv(i1)),'v','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfaO
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                    
                else
                    lg = [lg 'PBFA-O'];
                    lg_pbfaO = 1;
                end
                
                %disp('PBFA') 
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  NLDN2c(:,2)/1000) >= 1 ...
                    && sum(y(iv(i1)) == NLDN2c(:,3)/1000) >= 1
                hLine = plot(y(iv(i1)),z(iv(i1)),'<','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_nldn2c
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'NLDN2-C'];
                    lg_nldn2c = 1;
                end
                
                %disp('PBFA')
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  NLDN2g(:,2)/1000) >= 1 ...
                    && sum(y(iv(i1)) == NLDN2g(:,3)/1000) >= 1
                hLine = plot(y(iv(i1)),z(iv(i1)),'>','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_nldn2g
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'NLDN2-G'];
                    lg_nldn2g = 1;
                end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif sum(x(iv(i1)) ==  LINET(:,6)/1000) >= 1 ...
                    && sum(y(iv(i1)) == LINET(:,7)/1000) >= 1
                hLine = plot(y(iv(i1)),z(iv(i1)),'>','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_linet
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off')
                else
                    lg = [lg 'LINET'];
                    lg_linet = 1;
                end
           
            else
                %disp('Non')
                plot(y(iv(i1)),z(iv(i1)),'p','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1)
            end
                        
            
        end
        % plot(y(iv),z(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
    
    legend(lg)
    waitbar(0.6,wbh,'Plotting the graph','Name','Plotter Busy');
    
%     text(sen_set.y/1000,sen_set.z/1000,sen_set.sen_IDs)
%     plot(sen_set.y/1000,sen_set.z/1000,'rp','markerfacecolor','r')
%     
    
    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot(sen_set.y(ns)/1000,sen_set.z(ns)/1000,'kp','markerfacecolor','r')
            text(sen_set.y(ns)/1000,sen_set.z(ns)/1000,sen_set.sen_IDs{ns})
        end
    end
    
    % Plotting Additional land marks
    plot(lmy,lmz,'bp','markerfacecolor','b')
    text(lmy,lmz,lmID)
    
    waitbar(0.7,wbh,'Plotting the graph','Name','Plotter Busy');
    
    hold off
    
    ylabel('Altitude (km)')
    xlabel('North (km)')
    title(tit)
    box on
    
    waitbar(0.8,wbh,'Plotting the graph','Name','Plotter Busy');
    
    % Re-format the colorbar
    h=colorbar;

    %set(h,'ylim',[1 length(map)]);
    yal=linspace(1,length(map),10);
    set(h,'ytick',yal);
    % Create the yticklabels
    ytl=linspace(miv,mav,10);
    s=char(10,4);
    
    waitbar(0.9,wbh,'Plotting the graph','Name','Plotter Busy');
    
    for i=1:10
        if min(abs(ytl)) >= 0.001
            B=sprintf('%-4.3f',ytl(i));
        else
            B=sprintf('%-3.1E',ytl(i));
        end
        s(i,1:length(B))=B;
    end
    
    
    set(h,'yticklabel',s);
    
    waitbar(1.0,wbh,'Done','Name','Plotter Busy');
    delete(wbh)
    
    set(fg,'visible','on')
    
    
elseif seperate ~=1
    
    waitbar(0.95,wbh,'Plotting the graph','Name','Plotter Busy');
    
    subplot(3,2,4)
    hold on
    for nc=1:size(map,1)
        iv = find(t>=miv+(nc-1)*clrstep & t<=miv+nc*clrstep) ;
        plot(y(iv),z(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
    
%     plot(sen_set.y/1000,sen_set.z/1000,'rp','markerfacecolor','r')
%     text(sen_set.y/1000,sen_set.z/1000,sen_set.sen_IDs)

    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot(sen_set.y(ns)/1000,sen_set.z(ns)/1000,'kp','markerfacecolor','r')
            text(sen_set.y(ns)/1000,sen_set.z(ns)/1000,sen_set.sen_IDs{ns})
        end
    end
    
    % Plotting Additional land marks
    plot(lmy,lmz,'bp','markerfacecolor','b')
    text(lmy,lmz,lmID)
    
    hold off

    xlabel('North (km)')
    box on
    
    waitbar(1.0,wbh,'Done','Name','Plotter Busy');
    
    delete(wbh)
    
    set(fg,'visible','on')
    
    
end

if graphsOn(3) == 1
    plotter_tools(12)
end

if seperate == 1
    f = uimenu('Label','Plotter');
    uimenu(f,'Label','Box XY','Callback','plotter_tools(12)','Accelerator','E');
    uimenu(f,'Label','Fill Image','Callback','plotter_tools(13)','Accelerator','F');
    uimenu(f,'Label','Find Distance','Callback','plotter_tools(14)','Accelerator','D');
    
    if sen_set.radarOn
        uimenu(f,'Label','Vertical Radar Plane','Callback','plotter_tools(27)');
        uimenu(f,'Label','Radar altitude curser','Callback','radarTools(1)');
    else
        uimenu(f,'Label','Vertical Radar Plane','Callback','plotter_tools(27)','enable','off');
        uimenu(f,'Label','Radar altitude curser','Callback','radarTools(1)','enable','off');
    end
end


%% Changing the figure size
% Change the size of the figure
% if seperate~=1
%     units=get(fg,'units');
%     set(fg,'units','normalized','outerposition',[0.4 0.1 0.2 0.8]);
%     set(fg,'units',units);
% end


function format_colorbar(h,miv,mav)

% 5 ytick marks should be good.


cyl = get(h,'YLim');

set(h,'YTick',cyl(1):(cyl(2)-cyl(1))/10:cyl(2))

ytl=linspace(0,mav-miv,11);

s=char(11,4);

for i=1:11
    %if min(abs(ytl)) >= 0.001
        B=sprintf('%-4.3f',ytl(i));
    %else
    %    B=sprintf('%-3.1E',ytl(i));
    %end
    s(i,1:length(B))=B;
end


set(h,'yticklabel',s);
ylabel(h,sprintf('Seconds after %0.6f UT',miv))
