function LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,t1,t2,rc,x0,y0,z0,graph_type,sen_set,g)
%% LDAR-3D will plot LDAR points in space
% 
% syntax : 'LDAR_3D'  will plot points with color coded by time
%          'LDAR_3D(graph_number)' will plot the specific graph defined by
%          graph number.
%               if graph_number=1 ; Color coded by time
%               if graph_number=2 ; Color coded by type (CGLSS or DLS)
%               if graph_number=3 ; color coded by altitde
%               if graph_number=4 ; color coded by distance from the sensor
%                               5 ; color coded by flash type altitude vs time
%                               6 ; color coded by flash type North vs East
% Ex. 'LDAR_3D(3)'  will plot graph with color coded by altitude

%% Critical Inputs 
% fn='ldar2_20092280200.txt';             % File name for idar 
% 
% t1=7794; %7232;                         % Lower time bound
% t2=7795;
%rc=50000;                               % Bound radus 


%% Change default graph tipe here
% if nargin<1
%     graph_type=1;
%     % Default graph type : LDAR Points color coded by time
% else
%     graph_type=graph_number;
% end

%% Loading data

%% Check for missing files and exit if user prompts

absant_fn={};

fn=ldar_fn;

% Marker size
mz = 4;


%sen_set

if sen_set.ldarOn==1 || sen_set.cglssOn == 1
    if exist(ldar_fn,'file')==0
        absant_fn=[absant_fn fn];
        ldar_fn='';
    end
else
    ldar_fn='';
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


if sen_set.nldnOn==1
    if exist(nldn_fn,'file')==0
        absant_fn=[absant_fn nldn_fn];
        nldn_fn='';
    end
else
     nldn_fn='';
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
    
  
    % Plot title
    if sen_set.ldarOn
       tit=[tit ' LDAR'];
    else
        DLS=[NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
    end
    
    if sen_set.cglssOn
        tit=[tit ' CGLSS'];
    else
        CG=[NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
    end
    
    t=[CG(:,10);DLS(:,10)];         % time
    x=[CG(:,6);DLS(:,6)]/1000;      % East
    y=[CG(:,7);DLS(:,7)]/1000;      % North
    z=[CG(:,8);DLS(:,8)]/1000;      % Altitude
    r=[CG(:,11);DLS(:,11)]/1000;    % Distance from a sensor
    
    waitbar(0.08,wbh,'Loading LDAR data','Name','Plotter Busy');
else
    CG=[NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
    DLS=[NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
    waitbar(0.08,wbh,'Loading LDAR data','Name','Plotter Busy');
end


if isempty(linet_fn)==0
    waitbar(0.10,wbh,'Loading LINET data','Name','Plotter Busy');
    LINET=linetExtract(linet_fn,t1,t2,rc,x0,y0,z0);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if isempty(nldn_fn)==0
%     waitbar(0.151,wbh,'Loading LINET2 data','Name','Plotter Busy');
%     NLDN=linetExtract(nldn_fn,t1,t2,rc,x0,y0,z0);
%     waitbar(0.152,wbh,'Finish loading LINET2 data','Name','Plotter Busy');
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
%     waitbar(0.155,wbh,'Loading LINET2 data','Name','Plotter Busy');
% else
%         NLDN=[NaN NaN NaN NaN NaN NaN NaN NaN];
%         waitbar(0.155,wbh,'Loading LINET2 data','Name','Plotter Busy');
% end

if isempty(nldn_fn)==0
    waitbar(0.16,wbh,'Loading PBFA data','Name','Plotter Busy');
    PBFAA=pbfaExtract(nldn_fn,t1,t2,rc,x0,y0,z0);
    waitbar(0.18,wbh,'Loading PBFA data','Name','Plotter Busy');
    
    t=[t;PBFAA(:,6)];           % time
    x=[x;PBFAA(:,3)/1000];      % East
    y=[y;PBFAA(:,4)/1000];      % North
    z=[z;PBFAA(:,5)/1000];      % Altitude  
    r=[r;PBFAA(:,1)/1000];      % Distance from a sensor
    
   
    % Plot title
    tit=[tit ' PBFA-O'];
    
    waitbar(0.20,wbh,'Loading PBFA data','Name','Plotter Busy');
else
    PBFAA=[NaN NaN NaN NaN NaN NaN];
    waitbar(0.20,wbh,'Loading PBFA data','Name','Plotter Busy');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


if isempty(nldn2_fn)==0
    waitbar(0.16,wbh,'Loading NLDN2 data','Name','Plotter Busy');
    [NLDN2c NLDN2g] = nldnExtract(nldn2_fn,0,0,t1,t2,...
        rc,x0,y0,z0,0);
    waitbar(0.18,wbh,'Loading NLDN2 data','Name','Plotter Busy');
    
    t=[t;NLDN2c(:,8);NLDN2g(:,8)];           % time
    x=[x;NLDN2c(:,2)/1000;NLDN2g(:,2)/1000];      % East
    y=[y;NLDN2c(:,3)/1000;NLDN2g(:,3)/1000];      % North
    z=[z;NLDN2c(:,4)/1000;NLDN2g(:,4)/1000];      % Altitude  
    r=[r;NLDN2c(:,9)/1000;NLDN2g(:,9)/1000];      % Distance from a sensor
    
   
    % Plot title
    tit=[tit ' NLDN2'];
    
    waitbar(0.20,wbh,'Loading NLDN data','Name','Plotter Busy');
else
    PBFAO=[NaN NaN NaN NaN NaN NaN];
    waitbar(0.20,wbh,'Loading NLDN data','Name','Plotter Busy');
end


% Exit if there is no data in the given time range
if (length(t)-sum(isnan(t)))==0
    delete(wbh)
    errordlg('No data in the given time range.',...
        'LDAR 3D Error','modal')
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

if graph_type==1
    %% Plotting graph color coded by time
    waitbar(0.3,wbh,'Plotting the graph','Name','Plotter Busy');
%     t=[DLS(:,10);CG(:,10)];                       % time
%     x=[DLS(:,6);CG(:,6)]/1000;                    % East
%     y=[DLS(:,7);CG(:,7)]/1000;                    % North
%     z=[DLS(:,8);CG(:,8)]/1000;                   % Altitude
%     
    
    fg=figure;
    set(fg,'visible','off')
    
    map=colormap;
    miv=min(t);
    mav=max(t);
    
    waitbar(0.4,wbh,'Plotting the graph','Name','Plotter Busy');
    
    clrstep = (mav-miv)/size(map,1) ;
    
    waitbar(0.5,wbh,'Plotting the graph','Name','Plotter Busy');
    hold all
    
    lg_ldar = 0;
    lg_cglss = 0;
    lg_pbfa = 0;
    lg_pbfaA = 0;
    lg_pbfaO = 0;
    lg_nldn2c = 0;
    lg_nldn2g = 0;
    
    lg = {};
    
    
    for nc=1:size(map,1)
        iv = find(t>=miv+(nc-1)*clrstep & t<=miv+nc*clrstep) ;
        
        for i1 = 1:length(iv)
            %Is it a LDAR2
            if  sum(x(iv(i1)) ==  DLS(:,6)/1000) >= 1 ...                
                && sum(y(iv(i1)) == DLS(:,7)/1000) >= 1
                %disp('DLS')
                hLine = plot3(x(iv(i1)),y(iv(i1)),z(iv(i1)),'o','LineStyle','none','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
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
                hLine  = plot3(x(iv(i1)),y(iv(i1)),z(iv(i1)),'marker','s','LineStyle','none','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz);
                
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
                hLine = plot3(x(iv(i1)),y(iv(i1)),z(iv(i1)),'marker','d','LineStyle','none','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfa
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'PBFA'];
                    lg_pbfa = 1;
                end
                
                %disp('PBFA')
            % Is it a PBFA-Auto?
            elseif sum(x(iv(i1)) ==  PBFAA(:,3)/1000) >= 1 ...
                    && sum(y(iv(i1)) == PBFAA(:,4)/1000) >= 1
                hLine = plot3(x(iv(i1)),y(iv(i1)),z(iv(i1)),'marker','p','LineStyle','none','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfaA
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'PBFA-A'];
                    lg_pbfaA = 1;
                end
                
                %disp('PBFA')
            % Is it a PBFA-Old?
            elseif sum(x(iv(i1)) ==  PBFAO(:,3)/1000) >= 1 ...
                    && sum(y(iv(i1)) == PBFAO(:,4)/1000) >= 1
                hLine = plot3(x(iv(i1)),y(iv(i1)),z(iv(i1)),'marker','v','LineStyle','none','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_pbfaO
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'PBFA-O'];
                    lg_pbfaO = 1;
                end
                
                %disp('PBFA')
            
            %is it NLDNc
            elseif sum(x(iv(i1)) ==  NLDN2c(:,2)/1000) >= 1 ...
                    && sum(y(iv(i1)) == NLDN2c(:,3)/1000) >= 1
                hLine = plot3(x(iv(i1)),y(iv(i1)),z(iv(i1)),'marker','>','LineStyle','none','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_nldn2c
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'NLDN2-C'];
                    lg_nldn2c = 1;
                end
                
                %disp('NLDN2-C')
            
            %is it NLDNg
            elseif sum(x(iv(i1)) ==  NLDN2g(:,2)/1000) >= 1 ...
                    && sum(y(iv(i1)) == NLDN2g(:,3)/1000) >= 1
                hLine = plot3(x(iv(i1)),y(iv(i1)),z(iv(i1)),'marker','<','LineStyle','none','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1);
                
                if lg_nldn2g
                    % Exclude from legend
                    set(get(get(hLine,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');
                else
                    lg = [lg 'NLDN2-G'];
                    lg_nldn2g = 1;
                end
                
                %disp('NLDN2-C')
            else
                %disp('Non')

                plot3(x(iv(i1)),y(iv(i1)),z(iv(i1)),'marker','p','color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz+1)
            end
        end     
        
        
        %plot3(x(iv),y(iv),z(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
    legend(lg)
    waitbar(0.6,wbh,'Plotting the graph','Name','Plotter Busy');
    %plot3(3.380,4.410,0.002,'rp','markerfacecolor','r')
    %text(3.5,4.50,0.002,'Sensor')
   
    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot3(sen_set.x(ns)/1000,sen_set.y(ns)/1000,sen_set.z(ns)/1000,'rp','markerfacecolor','r')
            text(sen_set.x(ns)/1000,sen_set.y(ns)/1000,sen_set.z(ns)/1000,sen_set.sen_IDs{ns})
        end
    end
    
    % Plotting Additional land marks
    plot3(lmx,lmy,lmz,'bp','markerfacecolor','b')
    text(lmx,lmy,lmz,lmID)
    
    waitbar(0.7,wbh,'Plotting the graph','Name','Plotter Busy');
    hold off
    
    % Adjusting altitude limit
    if  z<12
        zlim([0 12])
    end
    
    xlabel('East (km)')
    ylabel('North (km) ')
    zlabel('Altitude (km)')
    waitbar(0.8,wbh,'Plotting the graph','Name','Plotter Busy');
    tit=sprintf('%s Color Coded by Time    \n%s-%s-%s     UT: %s - %s',...
        tit,fn(end-31:end-28),fn(end-26:end-25),fn(end-23:end-22),...
        sec2hhmmss(t1),sec2hhmmss(t2));

    title(tit)
    %grid on
    box on
    
    % Re-format the colorbar
    h=colorbar;
    waitbar(0.9,wbh,'Plotting the graph','Name','Plotter Busy');
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
    
    waitbar(1.0,wbh,'Done','Name','Plotter Busy');
    
    set(h,'yticklabel',s);
    set(get(h,'YLabel'),'String','Time (s)')    
    view(3)
    delete(wbh)
    set(fg,'visible','on')
    
    
  
elseif graph_type==2
    %% Plotting graph Color coded by flash type 
    waitbar(0.3,wbh,'Plotting the graph','Name','Plotter Busy');
    fg=figure;
    set(fg,'visible','off')
    
    plot3(DLS(:,6)/1000,DLS(:,7)/1000,DLS(:,8)/1000,'ko','markerfacecolor','k','MarkerSize',mz)
    hold on
    plot3(CG(:,6)/1000,CG(:,7)/1000,CG(:,8)/1000,'go','markerfacecolor','g','MarkerSize',mz)
    
    plot3(LINET(:,6)/1000,LINET(:,7)/1000,LINET(:,8)/1000,'ro','markerfacecolor','r','MarkerSize',mz)
    
    %plot3(NLDN(:,6)/1000,NLDN(:,7)/1000,NLDN(:,8)/1000,'bo','markerfacecolor','b','MarkerSize',mz)
    plot3(PBFAA(:,3)/1000,PBFAA(:,4)/1000,PBFAA(:,5)/1000,'bp','markerfacecolor','b','MarkerSize',mz+1)
    
    plot3(PBFA(:,3)/1000,PBFA(:,4)/1000,PBFA(:,5)/1000,'co','markerfacecolor','c','MarkerSize',mz)
    
    plot3(PBFAO(:,3)/1000,PBFAO(:,4)/1000,PBFAO(:,5)/1000,'mv','markerfacecolor','m','MarkerSize',mz+1)
    
    plot3(NLDN2c(:,2)/1000,NLDN2c(:,3)/1000,NLDN2c(:,4)/1000,'g>','markerfacecolor','g','MarkerSize',mz+1)
    plot3(NLDN2g(:,2)/1000,NLDN2g(:,3)/1000,NLDN2g(:,4)/1000,'r<','markerfacecolor','r','MarkerSize',mz+1)
    
    waitbar(0.4,wbh,'Plotting the graph','Name','Plotter Busy');
    
%     plot3(sen_set.x/1000,sen_set.y/1000,sen_set.z/1000,'rp','markerfacecolor','r')
%     text(sen_set.x/1000,sen_set.y/1000,sen_set.z/1000,sen_set.sen_IDs)

    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot3(sen_set.x(ns)/1000,sen_set.y(ns)/1000,sen_set.z(ns)/1000,'rp','markerfacecolor','r')
            text(sen_set.x(ns)/1000,sen_set.y(ns)/1000,sen_set.z(ns)/1000,sen_set.sen_IDs{ns})
        end
    end
      
    % Plotting Additional land marks
    plot3(lmx,lmy,lmz,'bp','markerfacecolor','b')
    text(lmx,lmy,lmz,lmID)
    
    waitbar(0.5,wbh,'Plotting the graph','Name','Plotter Busy');
    hold off
    
    lg={};
    
    if isnan(DLS(1,1))==0
        lg=[lg 'LDAR2'];
    end
    
    if isnan(CG(1,1))==0 
        lg=[lg 'CGLSS'];
    end
        
    if isnan(LINET(1,1))==0
        lg=[lg 'LINET'];
    end
    
%     if isnan(NLDN(1,1))==0
%         lg=[lg 'LINET2'];
%     end
    if isnan(PBFAA(1,1))==0
        lg=[lg 'PBFA-A'];
    end    

    if isnan(PBFA(1,1))==0
        lg=[lg 'PBFA'];
    end
    
    if isnan(PBFAO(1,1))==0
        lg=[lg 'PBFA-O'];
    end
    
    if isnan(NLDN2c(1,1))==0
        lg=[lg 'NLDN2-C'];
    end
    
    if isnan(NLDN2g(1,1))==0
        lg=[lg 'NLDN2-G'];
    end
    
    lg=[lg 'Sensors'];
    
    legend(lg)
    
    waitbar(0.6,wbh,'Plotting the graph','Name','Plotter Busy');
    
    tit=sprintf('%s Color Coded by Type    \n%s-%s-%s     UT: %s - %s',...
        tit,fn(end-31:end-28),fn(end-26:end-25),fn(end-23:end-22),...
        sec2hhmmss(t1),sec2hhmmss(t2));
    
    waitbar(0.8,wbh,'Plotting the graph','Name','Plotter Busy');
    
    title(tit)
    
    xlabel('East (km)')
    ylabel('North (km) ')
    zlabel('Altitude (km)')

    %grid on
    box on    
    waitbar(1.0,wbh,'Plotting the graph','Plotter Busy');
    
    delete(wbh)
    set(fg,'visible','on')
    

elseif graph_type==3
    %% Plotting the graph color coded by altitude
    
    waitbar(0.3,wbh,'Plotting the graph','Plotter Busy');
    
      
    
    fg=figure;
    set(fg,'visible','off')
    
    waitbar(0.5,wbh,'Plotting the graph','Plotter Busy');
    
    map=colormap;
    miv=min(z);
    mav=max(z);
    
    waitbar(0.4,wbh,'Plotting the graph','Plotter Busy');
    
    clrstep = (mav-miv)/size(map,1) ;
    
    marker='o';
    
    waitbar(0.6,wbh,'Plotting the graph','Plotter Busy');
    
    hold on
    for nc=1:size(map,1)
        iv = find(z>=miv+(nc-1)*clrstep & z<=miv+nc*clrstep) ;
        plot3(x(iv),y(iv),z(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
    
    waitbar(0.7,wbh,'Plotting the graph','Plotter Busy');
    
%     plot3(sen_set.x/1000,sen_set.y/1000,sen_set.z/1000,'rp','markerfacecolor','r')
%     text(sen_set.x/1000,sen_set.y/1000,sen_set.z/1000,sen_set.sen_IDs)
    
    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot3(sen_set.x(ns)/1000,sen_set.y(ns)/1000,sen_set.z(ns)/1000,'rp','markerfacecolor','r')
            text(sen_set.x(ns)/1000,sen_set.y(ns)/1000,sen_set.z(ns)/1000,sen_set.sen_IDs{ns})
        end
    end

    % Plotting Additional land marks
    plot3(lmx,lmy,lmz,'bp','markerfacecolor','b')
    text(lmx,lmy,lmz,lmID)
    
    waitbar(0.8,wbh,'Plotting the graph','Plotter Busy');
    hold off
    
    % Addusting altitude limit
    if  z<12
        zlim([0 12])
    end
    
    xlabel('East (km)')
    ylabel('North (km) ')
    zlabel('Altitude (km)')
    
    waitbar(0.85,wbh,'Plotting the graph','Plotter Busy');
    
    tit=sprintf('%s Color Coded by Altitude    \n%s-%s-%s     UT: %s - %s',...
        tit,fn(end-31:end-28),fn(end-26:end-25),fn(end-23:end-22),...
        sec2hhmmss(t1),sec2hhmmss(t2));

    title(tit)
    %grid on
    box on
    waitbar(0.9,wbh,'Plotting the graph','Plotter Busy');
        
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
            B=sprintf('%-4.2f',ytl(i));
        else
            B=sprintf('%-3.2E',ytl(i));
        end
        s(i,1:length(B))=B;
    end
    
    waitbar(0.95,wbh,'Plotting the graph','Plotter Busy');
    
    set(h,'yticklabel',s)
    set(get(h,'YLabel'),'String','Altitude (km)')

    view(3)
    waitbar(1.0,wbh,'Done','Plotter Busy');
   
    delete(wbh)
    set(fg,'visible','on')
    
  

elseif graph_type==4
    %% Plotting the graph color coded by Distance from the sensor
    
    waitbar(0.3,wbh,'Plotting the graph','Plotter Busy');
    
    
    
    
    fg=figure;
    set(fg,'visible','off')
    
    waitbar(0.5,wbh,'Plotting the graph','Plotter Busy');
    
    map=colormap;
    miv=min(r);
    mav=max(r);
    
    clrstep = (mav-miv)/size(map,1) ;
    
    waitbar(0.4,wbh,'Plotting the graph','Plotter Busy');
    
    marker='o';
    
    hold on
    for nc=1:size(map,1)
        iv = find(r>=miv+(nc-1)*clrstep & r<=miv+nc*clrstep) ;
        plot3(x(iv),y(iv),z(iv),marker,'color',map(nc,:),'markerfacecolor',map(nc,:),'MarkerSize',mz)
    end
    
    waitbar(0.6,wbh,'Plotting the graph','Plotter Busy');
    
    % Plotting the sensor point
%     plot3(sen_set.x/1000,sen_set.y/1000,sen_set.z/1000,'rp','markerfacecolor','r')
%     text(sen_set.x/1000,sen_set.y/1000,sen_set.z/1000,sen_set.sen_IDs)

    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot3(sen_set.x(ns)/1000,sen_set.y(ns)/1000,sen_set.z(ns)/1000,'rp','markerfacecolor','r','MarkerSize',mz)
            text(sen_set.x(ns)/1000,sen_set.y(ns)/1000,sen_set.z(ns)/1000,sen_set.sen_IDs{ns})
        end
    end
    
    % Plotting Additional land marks
    plot3(lmx,lmy,lmz,'bp','markerfacecolor','b')
    text(lmx,lmy,lmz,lmID)
    
    waitbar(0.7,wbh,'Plotting the graph','Plotter Busy');
    
    hold off
    
    % Addusting altitude limit
    if  z<12
        zlim([0 12])
    end
    
    xlabel('East (km)')
    ylabel('North (km) ')
    zlabel('Altitude (km)')
    
    tit=sprintf('%s Color Coded the distance from FM14   \n%s-%s-%s     UT: %s - %s',...
        tit,fn(end-31:end-28),fn(end-26:end-25),fn(end-23:end-22),...
        sec2hhmmss(t1),sec2hhmmss(t2));
    
    waitbar(0.8,wbh,'Plotting the graph','Plotter Busy');
    title(tit)
    %grid on
    box on
    
    waitbar(0.9,wbh,'Plotting the graph','Plotter Busy');
    
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
            B=sprintf('%-4.2f',ytl(i));
        else
            B=sprintf('%-3.2E',ytl(i));
        end
        s(i,1:length(B))=B;
    end
    set(h,'yticklabel',s)
    set(get(h,'YLabel'),'String','Distance from the sensor (km)')

    view(3)
    waitbar(1.0,wbh,'Done','Plotter Busy');
    
    delete(wbh)
    set(fg,'visible','on')

elseif graph_type==5
    
    %% Plotting graph Color coded by flash type Altitude vs. time
    waitbar(0.3,wbh,'Plotting the graph','Name','Plotter Busy');
    fg=figure;
    set(fg,'visible','off')
    
    plot(DLS(:,10),DLS(:,8)/1000,'ko','markerfacecolor','k','MarkerSize',mz)
    hold on
    plot(CG(:,10),CG(:,8)/1000,'go','markerfacecolor','g','MarkerSize',mz)
    
    plot(LINET(:,3),LINET(:,8)/1000,'ro','markerfacecolor','r','MarkerSize',mz)
    
    %plot(NLDN(:,3),NLDN(:,8)/1000,'bo','markerfacecolor','b','MarkerSize',mz)
    plot(PBFAA(:,6),PBFAA(:,5)/1000,'bp','markerfacecolor','b','MarkerSize',mz+1)
    
    plot(PBFA(:,6),PBFA(:,5)/1000,'co','markerfacecolor','c','MarkerSize',mz)
    
    plot(PBFAO(:,6),PBFAO(:,5)/1000,'mv','markerfacecolor','m','MarkerSize',mz+1)
    
    waitbar(0.4,wbh,'Plotting the graph','Name','Plotter Busy');
    
    
    waitbar(0.5,wbh,'Plotting the graph','Name','Plotter Busy');
    hold off
    
    lg={};
    
    if isnan(DLS(1,1))==0
        lg=[lg 'IC'];
    end
    
    if isnan(CG(1,1))==0 
        lg=[lg 'CG'];
    end
        
    if isnan(LINET(1,1))==0
        lg=[lg 'LINET'];
    end
    
%     if isnan(NLDN(1,1))==0
%         lg=[lg 'LINET2'];
%     end
    if isnan(PBFAA(1,1))==0
        lg=[lg 'PBFA-A'];
    end

    if isnan(PBFA(1,1))==0
        lg=[lg 'PBFA'];
    end
    
    if isnan(PBFAO(1,1))==0
        lg=[lg 'PBFA-O'];
    end
    
    legend(lg)
    
    waitbar(0.6,wbh,'Plotting the graph','Name','Plotter Busy');
    
    tit=sprintf('%s Color Coded by Type    \n%s-%s-%s     UT: %s - %s',...
        tit,fn(end-31:end-28),fn(end-26:end-25),fn(end-23:end-22),...
        sec2hhmmss(t1),sec2hhmmss(t2));
    
    waitbar(0.8,wbh,'Plotting the graph','Name','Plotter Busy');
    
    title(tit)
    
    xlabel('Time (s)')
    ylabel('Altitude(km) ')
    
    %grid on
    box on    
    waitbar(1.0,wbh,'Plotting the graph','Plotter Busy');
    
    delete(wbh)
    set(fg,'visible','on')

elseif graph_type==6
    
        
    %% Plotting graph Color coded by flash type North vs. East
    waitbar(0.3,wbh,'Plotting the graph','Name','Plotter Busy');
    fg=figure;
    hold all;
    set(fg,'visible','off')
    
    % Plot radar
    if sen_set.radarOn && ~isempty(sen_set.radarFn)
        radar_plot(sen_set,g)
    end
    
    plot(DLS(:,6)/1000,DLS(:,7)/1000,'ko','markerfacecolor','k','MarkerSize',mz)
    
    plot(CG(:,6)/1000,CG(:,7)/1000,'go','markerfacecolor','g','MarkerSize',mz)
    
    plot(LINET(:,6)/1000,LINET(:,7)/1000,'ro','markerfacecolor','r','MarkerSize',mz)
    
    %plot(NLDN(:,6)/1000,NLDN(:,7)/1000,'bo','markerfacecolor','b','MarkerSize',mz)
    plot(PBFAA(:,3)/1000,PBFAA(:,4)/1000,'bp','markerfacecolor','b','MarkerSize',mz+1)
    
    plot(PBFA(:,3)/1000,PBFA(:,4)/1000,'co','markerfacecolor','c','MarkerSize',mz)
    
    plot(PBFAO(:,3)/1000,PBFAO(:,4)/1000,'mv','markerfacecolor','m','MarkerSize',mz+1)
    
    waitbar(0.4,wbh,'Plotting the graph','Name','Plotter Busy');
    
    % Plotting the sensor point
%     plot(sen_set.x/1000,sen_set.y/1000,'rp','markerfacecolor','r')
%     text(sen_set.x/1000,sen_set.y/1000,sen_set.sen_IDs)
    for ns = 1:20
        if sen_set.x(ns)~=0 || sen_set.y(ns)~=0 || sen_set.z(ns)~=0
            plot(sen_set.x(ns)/1000,sen_set.y(ns)/1000,'rp','markerfacecolor','r')
            text(sen_set.x(ns)/1000,sen_set.y(ns)/1000,sen_set.sen_IDs{ns})
        end
    end

    
    % Plotting Additional land marks
    plot(lmx,lmy,'bp','markerfacecolor','b')
    text(lmx,lmy,lmID)
    
    
    
    waitbar(0.5,wbh,'Plotting the graph','Name','Plotter Busy');
    hold off
    
    lg={};
    
    if isnan(DLS(1,1))==0
        lg=[lg 'LDAR2'];
    end
    
    if isnan(CG(1,1))==0 
        lg=[lg 'CGLSS'];
    end
        
    if isnan(LINET(1,1))==0
        lg=[lg 'LINET'];
    end
    
%     if isnan(NLDN(1,1))==0
%         lg=[lg 'LINET2'];
%     end
    
    if isnan(PBFAA(1,1))==0
        lg=[lg 'PBFA-A'];
    end

    if isnan(PBFA(1,1))==0
        lg=[lg 'PBFA'];
    end
    
    if isnan(PBFAO(1,1))==0
        lg=[lg 'PBFA-O'];
    end
    
    
    legend(lg)
    
    waitbar(0.6,wbh,'Plotting the graph','Name','Plotter Busy');
    
    tit=sprintf('%s Color Coded by Type    \n%s-%s-%s     UT: %s - %s',...
        tit,fn(end-31:end-28),fn(end-26:end-25),fn(end-23:end-22),...
        sec2hhmmss(t1),sec2hhmmss(t2));
    
    waitbar(0.8,wbh,'Plotting the graph','Name','Plotter Busy');
    
    title(tit)
    
    xlabel('East (km)')
    ylabel('North (km) ')
    
    %grid on
    box on 
    grid on
    waitbar(1.0,wbh,'Plotting the graph','Plotter Busy');
    
    if sen_set.mapOn
        xdd = [x; -60.1; 3.4];
        ydd = [y; -53.1; 49.3];
        
        xlim([min(xdd) max(xdd)]*1.3)
        ylim([min(ydd) max(ydd)]*1.3)
        
        florida_map
    end
    
    delete(wbh)
    
    set(fg,'visible','on')
end

if graph_type ~= 5
    plotter_tools(12)
end

f = uimenu('Label','Plotter');
uimenu(f,'Label','Box XY','Callback','plotter_tools(12)','Accelerator','E');
uimenu(f,'Label','Fill Image','Callback','plotter_tools(13)','Accelerator','F');
uimenu(f,'Label','Find Distance','Callback','plotter_tools(14)','Accelerator','D');



