function NBP_analasys_post_process2

%% file name for data
%fn = 'E:\Sumedhe\Documents\NBP\20110814-NBP_info2.xlsx';
fn = 'C:\Users\sumedhe\Desktop\Nadee\Nadee-official-sumpc\NBP- vertical radar scan\2011-08-14\20110814-NBP_info4.xlsx';

% base folder for saving plots
 bf = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\Results\';

% Save plots?
savePlots = 0;

sns_x = [ -0.7524, 0.3372, -0.3149, -1.1658, -0.2701, -2.7640, -6.0091, 0.1825, -5.7394, -2.0637]*10;
sns_y = [  1.6555, 0.4446, -0.6838, -1.7020,  0.2631,  4.9254, -3.3983, -5.3008, 1.1923, 0.1569]*10;
sns_n = { 'K02', 'K14',    'K24',    'BCC',    'K17',    'EDW',    'STC',    'FLT',    'OVD',    'FFI'};

try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

sen_set = h.sen_set;
g = h.g;

sheet = 1;
xlRange = 'A2:AD306';

[ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
% 
pbfat = ndata(:,10);
pbfax = ndata(:,11)/1000;
pbfay = ndata(:,12)/1000;
pbfaz = ndata(:,13)/1000;
type = ndata(:,17);
% figure
% 
% % Generate xy radar plot and plot PBFA points according to their types
% %hist(pbfat,[0:300:84600])
% 
% 
% startT =80687.7413463;
% endT = 80879.6580636;
% 
%  t1 = 80687.7413463;  t2 = 80879.6580636; %1
% % t1 = 80907;    t2 = 81163; %2
% % t1 = 81163;    t2 = 81419.5; %3
% % t1 = 81419.5;  t2 = 81676.5; %4
% % t1 = 81676.5;   t2 = 81933; %5
% % t1 = 81933;     t2 = 82189.5; %6
% 
%  lol = sum(pbfat < t1) + 1;
%  ul = sum(pbfat < t2);
% 
% inds1 = find(pbfat > t1);
% inds2 = find(pbfat < t2);
% 
% inds = intersect(inds1,inds2);
% 
% rdfns = [];
% times = [];
% 
% for i=1:24
%     %bf = sprintf('H:/data/netCDF/MEL/2011/08/14/%2.2i/',i);
% %     bf = sprintf('C:/data/2011-08-05 -- 2011-08-16/data/netCDF/MEL/2011/08/14/%2.2i/',i);
%     bf = sprintf('//SADAQ7/data/newRaid2/data/netCDF/MEL/2011/08/14/%2.2i/',i);
%     files = ls(bf);
% 
%     [L L2] = size(files);
% 
%     for j = 3:L
%         rdfns = [rdfns; [bf files(j,:)]];
%         times = [times; str2num(files(j,end-16:end-15))*3600+...
%             str2num(files(j,end-14:end-13))*60 + ...
%             str2num(files(j,end-12:end-11))];
% 
%     end
% 
% end

% 
% % start times and end times
% dt = (times(1:end-1)+times(2:end))/2;
% 
% stT = [times(1); dt];
% enT = [dt; times(end)];
% 
% 
% % plot radar data for each file
% L = length(times);
% 
% % vericle radar base data [index x1 y1 x2 y2]
% %vrbd(362,:) = [14.16 -9.98 29.92 -4.36];
% 
% %insert x,y values of line
%      vrbd(368,:) = [13.25 -28.92 36.35 -21.08  ];
% 
% %for i = 1:L
% % 352-368
% for i = 363
%     figure
%     hold all
%     sen_set.radarFn = rdfns(i,:);
%     
%     %radar elevation angle
%     sen_set.radEleAngInd = 12;
% 
%     radar_plot(sen_set,g)
% 
% 
%     % Plot NBP data
%    ind1 = find(pbfat > stT(i));
%    ind2 = find(pbfat < enT(i));
% 
% 
%     inds0 = intersect(ind1,ind2);
% 
%     LNBP = length(inds0);
% 
%     % Type1
%     typeInds = find(type == 1);
%     inds = intersect(typeInds,inds0);
%     plot(pbfax(inds),pbfay(inds),'ko','markerfacecolor','m');
%     text(pbfax(inds)+0.5,pbfay(inds),num2str(inds),'FontWeight','bold');
% 
% 
% 
%     % Type2
%     typeInds = find(type == 2);
%     inds = intersect(typeInds,inds0);
%     plot(pbfax(inds),pbfay(inds),'ko','markerfacecolor','c');
%     text(pbfax(inds)+0.5,pbfay(inds),num2str(inds),'FontWeight','bold')
% 
% 
%     % Type3
%     typeInds = find(type == 3);
%     inds = intersect(typeInds,inds0);
%     plot(pbfax(inds),pbfay(inds),'ko','markerfacecolor','k');
%     text(pbfax(inds)+0.5,pbfay(inds),num2str(inds),'FontWeight','bold')
% 
% 
% 
%     % setup
%     daspect([1 1 1])
%     xlabel('East (km)')
%     ylabel('North (km)')
% 
%     plot(sns_x,sns_y,'kp','markerfacecolor','k')
%     florida_map
%     box on
%     saveFn = [rdfns(i,end-29:end) ' nNBPs = ' num2str(LNBP)]
%     title(saveFn ,'interpreter','none')
% 
%     xlim([0 40])
%     ylim([-30 10])
% 
%     tools2fig
%     f = uimenu('Label','Radar');
%     uimenu(f,'Label','Vertical Radar Plane','Callback','plotter_tools(27)');
%     uimenu(f,'Label','Radar altitude curser','Callback','radarTools(1)');
%     uimenu(f,'Label','Radar Contour','Callback','radar_contour');
% 
%     % Plot vertical radar plane line
%     plot([vrbd(i,1),vrbd(i,3)],[vrbd(i,2),vrbd(i,4)],'k', 'LineWidth',2)
% 
%     %saveas(gcf,['C:\Users\sumedhe\Desktop\Nadee\Nadee-official-sumpc\NBP- vertical radar scan\2011-08-14\radSeq\' saveFn '.fig']);
%     %saveas(gcf,['C:\Users\sumedhe\Desktop\Nadee\Nadee-official-sumpc\NBP- vertical radar scan\2011-08-14\radSeq\' saveFn '.png']);
%     %delete(gcf)
% 
%     %% Plotting Verticle radar plane
%     plot_vert_radar_plane2(rdfns(i,:),vrbd(i,1),vrbd(i,2),vrbd(i,3),vrbd(i,4),0)
%     figH = gcf;
%     radarTools(4) % for smoothed radar data
%     %radarTools(2) % for real 2D radar data
%     delete(figH)
% 
% 
%     % Type1
%     typeInds = find(type == 1);
%     inds = intersect(typeInds,inds0);
%     r = convert2D(vrbd(i,:),pbfax(inds),pbfay(inds));
%     plot(r,pbfaz(inds),'ko','markerfacecolor','m');
%     text(r+0.5,pbfaz(inds),num2str(inds),'FontWeight','bold')
% 
%     % Type2
%     typeInds = find(type == 2);
%     inds = intersect(typeInds,inds0);
%     r = convert2D(vrbd(i,:),pbfax(inds),pbfay(inds));
%     plot(r,pbfaz(inds),'ko','markerfacecolor','c');
%     text(r+0.5,pbfaz(inds),num2str(inds),'FontWeight','bold')
% 
%     % Type3
%     typeInds = find(type == 3);
%     inds = intersect(typeInds,inds0);
%     r = convert2D(vrbd(i,:),pbfax(inds),pbfay(inds));
%     plot(r,pbfaz(inds),'ko','markerfacecolor','k');
%     text(r+0.5,pbfaz(inds),num2str(inds),'FontWeight','bold')
% 
%     % title
%     [rdn, rfn ,extn ] = fileparts(rdfns(i,:));
%     title(sprintf('%s%s\n%s -- %s UT', rfn,extn,sec2hhmmss(stT(i)),sec2hhmmss(enT(i))),...
%     'Interpreter','none')
% 
% end


%% Generate radar verticle scans to determine NBP types

% Get all radar filenames
rdfns = [];
times = [];

for i=19:20
%     bf = sprintf('C:/data/2011-08-05 -- 2011-08-16/data/netCDF/MEL/2011/08/14/%2.2i/',i);
    bf = sprintf('T:/2011-08-05 -- 2011-08-16/data/netCDF/MEL/2011/08/14/%2.2i/',i);
    %bf = sprintf('//SADAQ7/data/newRaid2/data/netCDF/MEL/2011/08/14/%2.2i/',i);
    files = ls(bf);
    
    [L L2] = size(files);
    
    for j = 3:L
        rdfns = [rdfns; [bf files(j,:)]];
        times = [times; str2num(files(j,end-16:end-15))*3600+...
            str2num(files(j,end-14:end-13))*60 + ...
            str2num(files(j,end-12:end-11))];
        
    end
    
end


% start times and end times
dt = (times(1:end-1)+times(2:end))/2;

stT = [times(1); dt];
enT = [dt; times(end)];


for i = 1:length(pbfat)
    pbfat(i)
    if ~isnan(pbfat(i))
        % identify radar filename
        ind1 = find(stT <= pbfat(i));
        ind2 = find(enT > pbfat(i));
        ind = intersect(ind1,ind2);
                
        % Generate horizontal plot
        sen_set.radarFn = rdfns(ind,:);
        sen_set.radEleAngInd = 5;
        sbfh = figure;
        hrdh = figure;
        radar_plot(sen_set,g)
        hold all
               
        % Plot PBFA point
        plot(pbfax(i),pbfay(i),'ko','markerfacecolor','k')
        text(pbfax(i)+2,pbfay(i),num2str(i))        
        haxh = gca;
        
        % Let's generate vertical radar figures
        angs = [0 22.5 45 67.5 90 112.5 135 157.5 ]*pi/180;
        
        for k = 1:length(angs)
            
            ang = angs(k);
            
            figure(hrdh)
            x1 = pbfax(i) - 20*cos(ang);
            x2 = pbfax(i) + 20*cos(ang);
            y1 = pbfay(i) - 20*sin(ang);
            y2 = pbfay(i) + 20*sin(ang);
           
            try
                plot_vert_radar_plane2(rdfns(ind,:),x1,y1,x2,y2,0)
                figH = gcf;
                radarTools(4) % for smoothed radar data
                %radarTools(2) % for real 2D radar data
                delete(figH)
            end
            
            % Plot NBP point
            r = convert2D([x1 y1 x2 y2],pbfax(i),pbfay(i));
            plot(r,pbfaz(i),'ko','markerfacecolor','k');
            text(r+0.5,pbfaz(i),num2str(i),'FontWeight','bold')
            fgh = gcf;
            axh = gca;
            
            % Insert the figure into a subplot
            figure(sbfh)
            sbh = subplot(3,3,k);
            copyobj(get(axh,'children'),sbh)
            
            xlabel(get(get(axh,'xlabel'),'string'))
            ylabel(get(get(axh,'ylabel'),'string'))
            
            box on
            daspect([1 1 1])
            set(gcf,'renderer','painters')
            
            
            
            % cleanup
            delete(fgh)
        end
   
        
        insert last figure as the horizontal one
        figure(sbfh)
        sbh = subplot(3,3,9);
        
        copyobj(get(haxh,'children'),sbh)
        
        hold all
        setup_fig(rdfns(ind,:),[pbfax(i) - 40 pbfax(i) + 40] ,[pbfay(i) - 40,pbfay(i) + 40])
        
        box on
        daspect([1 1 1])
        set(gcf,'renderer','painters')
        
        % save as matlab figure
        fn = sprintf('C:/Users/sumedhe/Desktop/VerticleRadarNBP2/%3.3i',i);
        %saveas(gcf,[fn '.fig'])
        
        
        
        
        set(gca, 'Color', 'none')
        
        % Setup for save and then save as fig and png
        figure(sbfh)
        %set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 1280 1024]/200);
        set(gcf,'units','normalized','outerposition',[0 0 1 1])  
        %export_fig(fn, '-png' ,'-transparent')

        % cleanup
        %delete(sbfh)
        %delete(hrdh)
        return
    end
end

%% Convert NBP xyz to rz and plot 
% %(this program is for  generating several NBPS in a same vertical radr plane)
% 
% % % User inputs
%  xyz = [-10046.1 49812.5 13826 ]; % In meters (NBP point)
% %  %  after generating the intermediate figure before vertival radar, type >> guidata(gcf)
% x1= -14.3860; % In km (cordinates of line on wich vertical rada will be generated)
% x2= 0.5614;
% y1= 53.3333;
% y2= 41.8246;
% 
% % % Program
% xyz = xyz/1000;
% r = convert2D([x1 y1 x2 y2],xyz(1),xyz(2));
% 
% figure(gcf)
% hold all
% plot(r,xyz(3),'dk')

%% Plot LDAR2 in a time range
% 
% % This program will plot ldar2 points in given time range on already plotted
% % veritical scan. You have to choose the timinng range from Plotter 2
% % first.
% 
% % User inputs
% % AA' for NBP cluster(7)
% % x1= -39.0695;
% % x2= -24.6951;
% % y1= 44.8199;
% % y2= 29.4038;
% % AA' for zNBP quadreplet (old line)
% % x1= -24.5568;
% % x2= 5.2841;
% %y1= 56.5568;
% % y2= 42.3750;
% 
% % % AA' for zNBP quadreplet (new line)
% %  x1= -14.3860;
% % x2= 0.5614;
% % y1= 53.3333;
% % y2= 41.8246;
% % % 
% % % 
% % AA through NBP# 206,212 and IC closer to them
%  x1= 15.0766;
%  x2= -10.8705;
%  y1= -14.3154;
%  y2= 2.6066;
%  
%  
%  LDAR_color = 'r'; %[51 204 255]/255;
%  CGLSS_color = 'r';
%  PBFA_color = 'r';
%  markerSize = 3;
% % % 
% % % % Get plotter 2 data
% try h=guidata(findall(0,'Tag','plotter2'));
% catch; disp('Run plotter2 first! Kaniye'); return;
% end
% % 
% g = h.g;
% settings = h.sen_set;
% loc_dir = settings.loc_dir;
% 
% if g.mm < 30
%     ext=0;
% else
%     ext=30;
% end
% % 
% dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
%     -datenum(str2double(g.YYYY),0,0));
% 
% ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
%     loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);
% 
% % Time shift for LDAR and Linet
% if settings.ldar_tshiftOn==1
%     sn=settings.ldar_tshift_sn;
%     x=settings.x;
%     y=settings.y;
%     z=settings.z;
%     x0=x(sn)-settings.x0;
%     y0=y(sn)-settings.y0;
%     z0=z(sn)-settings.z0;
% else
%     x0=0;
%     y0=0;
%     z0=0;
% end
% 
% [CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
%         x0,y0,z0,0);
%     
% pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
%         loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
%         g.MM{:},g.DD{:},g.hh,ext);
%     
% PBFA=pbfaExtract(pbfa_fn,g.t1,g.t2,str2double(settings.ldar_r),...
%         x0,y0,z0);
% %     
% % % Plot data
%  figure(gcf)
% % %clf
% hold all
% r = convert2D([x1,y1,x2,y2],DLS(:,6)/1000,DLS(:,7)/1000);
% plot(r,DLS(:,8)/1000,'o','markerfacecolor',LDAR_color,...
%     'color',LDAR_color,'markersize',markerSize)
% 
% r = convert2D([x1,y1,x2,y2],CG(:,6)/1000,CG(:,7)/1000);
% plot(r,CG(:,8)/1000,'s','markerfacecolor',CGLSS_color,...
%     'color',CGLSS_color,'markersize',markerSize)
% 
% r = convert2D([x1,y1,x2,y2],PBFA(:,3)/1000,PBFA(:,4)/1000);
% plot(r,PBFA(:,5)/1000,'p','markerfacecolor',PBFA_color,...
%     'color',PBFA_color,'markersize',markerSize)


%% Plot horizontal xy plot for given time range 
% Requiered - (1) Time range will be used from plotter2
%             (2) Plot the xy figure LDAR/CGLSS plot
%             (3) Make the figure "current" by clicking on it.

% LDAR_color = 'm'; %[51 204 255]/255;
% CGLSS_color = 'm';
% PBFA_color = 'm';
% markerSize = 1;
% % 
% % % Get plotter 2 data
% try h=guidata(findall(0,'Tag','plotter2'));
% catch; disp('Run plotter2 first! Kaniye'); return;
% end
% 
% g = h.g;
% settings = h.sen_set;
% loc_dir = settings.loc_dir;
% 
% if g.mm < 30
%     ext=0;
% else
%     ext=30;
% end
% 
% dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
%     -datenum(str2double(g.YYYY),0,0));
% 
% ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
%     loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);
% 
% % Time shift for LDAR and Linet
% if settings.ldar_tshiftOn==1
%     sn=settings.ldar_tshift_sn;
%     x=settings.x;
%     y=settings.y;
%     z=settings.z;
%     x0=x(sn)-settings.x0;
%     y0=y(sn)-settings.y0;
%     z0=z(sn)-settings.z0;
% else
%     x0=0;
%     y0=0;
%     z0=0;
% end
% 
% [CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
%         x0,y0,z0,0);
%     
% pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
%         loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
%         g.MM{:},g.DD{:},g.hh,ext);
%     
% PBFA=pbfaExtract(pbfa_fn,g.t1,g.t2,str2double(settings.ldar_r),...
%         x0,y0,z0);
%     
%     
% figure(gcf); hold all;
% plot(DLS(:,6)/1000,DLS(:,7)/1000,'o','color',LDAR_color,'markerfacecolor',LDAR_color)
% plot(CG(:,6)/1000,CG(:,7)/1000,'s','color',CGLSS_color,'markerfacecolor',CGLSS_color)
% plot(PBFA(:,3)/1000,PBFA(:,4)/1000,'d','color',PBFA_color,'markerfacecolor',PBFA_color)


% %% Get info about doublets triplets etc.
% % Written on 2015/06/12
% clc
% fileName = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\NBP- vertical radar scan\2011-08-14\20110814-NBP_info24.xlsx';
% sheet = 2;
% xlRange = 'A2:AD139';
% [ndata, txtdata, data] = xlsread(fileName, sheet, xlRange);
% 
% index = ndata(:,1);
% pbfat = ndata(:,12);
% pbfax = ndata(:,13);
% pbfay = ndata(:,14);
% pbfaz = ndata(:,15);
% 
% sheet = 3;
% xlRange = 'A2:AD31';
% [ndata, txtdata, data] = xlsread(fileName, sheet, xlRange);
% 
% index = [index; ndata(:,1)];
% pbfat = [pbfat; ndata(:,10)];
% pbfax = [pbfax; ndata(:,11)];
% pbfay = [pbfay; ndata(:,12)];
% pbfaz = [pbfaz; ndata(:,13)];
% 
% sheet = 4;
% xlRange = 'A2:AD8';
% [ndata, txtdata, data] = xlsread(fileName, sheet, xlRange);
% 
% index = [index; ndata(:,1)];
% pbfat = [pbfat; ndata(:,10)];
% pbfax = [pbfax; ndata(:,11)];
% pbfay = [pbfay; ndata(:,12)];
% pbfaz = [pbfaz; ndata(:,13)];
% 
% 
% % Doublets
% dbI = [22	31
% 45	49
% 51	56
% 57	62
% 90	92
% 116	121
% 117	119
% 137	139
% 148	180
% 149	290
% 160	164
% 163	170
% 174	177
% 194	196
% 204	234
% 206	212
% 220	221
% 244	245
% 250	251
% 258	259
% 270	273
% 272	275
% ];
% 
% dts = [];
% dxs = [];
% dys = [];
% dzs = [];
% 
% 
% for i = 1:length(dbI)
%     ind1 = find(index == dbI(i,1));
%     ind2 = find(index == dbI(i,2));
%     
%     dts = [dts abs(pbfat(ind1) - pbfat(ind2))];
%     dxs = [dxs abs(pbfax(ind1) - pbfax(ind2))];
%     dys = [dys abs(pbfay(ind1) - pbfay(ind2))];
%     dz = abs(pbfaz(ind1) - pbfaz(ind2));
%     dzs = [dzs dz];
%     
%     if dz > 1500
%         fprintf('%i \t%i\n',dbI(i,1),dbI(i,2))
%     end
%     
% end
% 
% dxs = dxs/1000;
% dys = dys/1000;
% dzs = dzs/1000;
% fprintf('                 \t       \t min \t max \t mean \t std \t median \n')
% fprintf('Doublets (N = %i)\t dt(s)\t %0.1f\t %0.1f \t %0.1f \t %0.1f \t %0.1f \n',...
%     length(dbI),min(dts),max(dts),mean(dts),std(dts),median(dts))
% fprintf('                 \t dx(km)\t %0.2f\t %0.2f \t %0.2f \t %0.2f \t %0.2f \n',...
%     min(dxs),max(dxs),mean(dxs),std(dxs),median(dxs))
% fprintf('                 \t dy(km)\t %0.2f\t %0.2f \t %0.2f \t %0.2f \t %0.2f \n',...
%     min(dys),max(dys),mean(dys),std(dys),median(dys))
% fprintf('                 \t dz(km)\t %0.2f\t %0.2f \t %0.2f \t %0.2f \t %0.2f \n',...
%     min(dzs),max(dzs),mean(dzs),std(dzs),median(dzs))
% 
% 
% % Triplets
% trI = [  13	17 19
%     113	115 116
% 222	226	228
% ];
%      
% dts = [];
% dxs = [];
% dys = [];
% dzs = [];
% 
% for i = 1:length(trI)
%     inds(1) = find(index == trI(i,1));
%     inds(2) = find(index == trI(i,2));
%     inds(3) = find(index == trI(i,3));
%     
%        
%     for j = 2:3
%         dts = [dts abs(pbfat(inds(1)) - pbfat(inds(j)))];
%         dxs = [dxs abs(pbfax(inds(1)) - pbfax(inds(j)))];
%         dys = [dys abs(pbfay(inds(1)) - pbfay(inds(j)))];
%         dz = abs(pbfaz(inds(1)) - pbfaz(inds(j)));
%         dzs = [dzs dz];
%         
%         if dz > 1600
%             fprintf('%i\t%i\t%i\n',trI(i,1),trI(i,2),trI(i,3))
%         end
%     end
% end
% 
% dxs = dxs/1000;
% dys = dys/1000;
% dzs = dzs/1000;
% 
% fprintf('Triplets (N = %i)\t dt(s)\t %0.1f\t %0.1f \t %0.1f \t %0.1f \t %0.1f \n',...
%     length(trI),min(dts),max(dts),mean(dts),std(dts),median(dts))
% fprintf('                 \t dx(km)\t %0.2f\t %0.2f \t %0.2f \t %0.2f \t %0.2f \n',...
%     min(dxs),max(dxs),mean(dxs),std(dxs),median(dxs))
% fprintf('                 \t dy(km)\t %0.2f\t %0.2f \t %0.2f \t %0.2f \t %0.2f \n',...
%     min(dys),max(dys),mean(dys),std(dys),median(dys))
% fprintf('                 \t dz(km)\t %0.2f\t %0.2f \t %0.2f \t %0.2f \t %0.2f \n',...
%     min(dzs),max(dzs),mean(dzs),std(dzs),median(dzs))
% 
% 
% % Qdplets
% qdI = [9	11  57  62
%     20	21  30 	53
% ];
% 
% dts = [];
% dxs = [];
% dys = [];
% dzs = [];
% 
% for i = 1:2
%     inds(1) = find(index == qdI(i,1));
%     inds(2) = find(index == qdI(i,2));
%     inds(3) = find(index == qdI(i,3));
%     inds(4) = find(index == qdI(i,4));
%     
%   for j = 2:4
%         dts = [dts abs(pbfat(inds(1)) - pbfat(inds(j)))];
%         dxs = [dxs abs(pbfax(inds(1)) - pbfax(inds(j)))];
%         dys = [dys abs(pbfay(inds(1)) - pbfay(inds(j)))];
%         dz = abs(pbfaz(inds(1)) - pbfaz(inds(j)));
%         dzs = [dzs dz];
%         
%         if dz > 1600
%             fprintf('%i\t%i\t%i\n',trI(i,1),trI(i,2),trI(i,3))
%         end
%     end
% end
% dxs = dxs/1000;
% dys = dys/1000;
% dzs = dzs/1000;
% fprintf('Quadraplets (N = %i)\t dt(s)\t %0.1f\t %0.1f \t %0.1f \t %0.1f \t %0.1f \n',...
%     2,min(dts),max(dts),mean(dts),std(dts),median(dts))
% fprintf('                 \t dx(km)\t %0.2f\t %0.2f \t %0.2f \t %0.2f \t %0.2f \n',...
%     min(dxs),max(dxs),mean(dxs),std(dxs),median(dxs))
% fprintf('                 \t dy(km)\t %0.2f\t %0.2f \t %0.2f \t %0.2f \t %0.2f \n',...
%     min(dys),max(dys),mean(dys),std(dys),median(dys))
% fprintf('                 \t dz(km)\t %0.2f\t %0.2f \t %0.2f \t %0.2f \t %0.2f \n',...
%     min(dzs),max(dzs),mean(dzs),std(dzs),median(dzs))
% 
% 
% % Sedtuplet
% sxI = [124	125	127 128 129 132];
% dts = nan(1,5);
% dxs = nan(1,5);
% dys = nan(1,5);
% dzs = nan(1,5);
% 
% ind1 =  find(index == sxI(1));
% 
% for i = 1:5
%     ind2 = find(index == sxI(i+1));
%     dts(i) = abs(pbfat(ind1) - pbfat(ind2));
%     dxs(i) = abs(pbfax(ind1) - pbfax(ind2));
%     dys(i) = abs(pbfay(ind1) - pbfay(ind2));
%     dzs(i) = abs(pbfaz(ind1) - pbfaz(ind2));
% end
%     
% dxs = dxs/1000;
% dys = dys/1000;
% dzs = dzs/1000;
% fprintf('Sextuplet (N = %i)\t dt(s)\t %0.1f\t %0.1f \t %0.1f \t %0.1f \t %0.1f \n',...
%     1,min(dts),max(dts),mean(dts),std(dts),median(dts))
% fprintf('                 \t dx(km)\t %0.2f\t %0.2f \t %0.2f \t %0.2f \t %0.2f \n',...
%     min(dxs),max(dxs),mean(dxs),std(dxs),median(dxs))
% fprintf('                 \t dy(km)\t %0.2f\t %0.2f \t %0.2f \t %0.2f \t %0.2f \n',...
%     min(dys),max(dys),mean(dys),std(dys),median(dys))
% fprintf('                 \t dz(km)\t %0.2f\t %0.2f \t %0.2f \t %0.2f \t %0.2f \n',...
%     min(dzs),max(dzs),mean(dzs),std(dzs),median(dzs))
% 


%% Plotting the map
% map_file = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\NBP- vertical radar scan\LatexDocument\pics\Map\base_map.jpg';
% I = imread(map_file);
% 
% figure
% imagesc(I)
% hold all
% LL1 = [ -82.402481,27.704536];
% xy1 = [251 622];
% 
% LL2 = [ -80.60811,28.57854];
% xy2 = [578, 441];
% 
% xc = -86.0:0.5:-78.6;
% nx = xy1(1) + (xc - LL1(1))*(xy2(1) - xy1(1))/(LL2(1)-LL1(1));
% 
% set(gca,'xtick',nx,'xticklabel',xc)
% 
% grid on
% daspect([1 1 1])


%% Plot vertical radar data
% This program is for generating the same line(AA') in different radar PPI , so
% tht vertical radar plane generated will be the same for each time period

% % % User input
   a = guidata(gcf);
  fn = a.fileName;
% % cordinates ofvertical radar line (AA') obtained using >>guidata(gcf)
% % AA' cordinates for: NBP cluster 124,125,126,127,128,129,132
% x1= -39.0695;
% x2= -24.6951;
% y1= 44.8199;
% y2= 29.4038;
% 
%  % AA' cordinates for: NBP quadreplet 20,21,30,53
% x1= -14.3860;
% x2= 0.5614;
% y1= 53.3333;
% y2= 41.8246;


% cordinates for over the flash in '2011-07-22-Originationof flash project'
 x1= 16.1728;
 x2= 43.4568;
 y1= 38.3333;
 y2= 5.3704;


%AA through NBP# 206,212 and IC closer to them
%  x1= 15.0766;
%  x2= -10.8705;
%  y1= -14.3154;
%  y2= 2.6066;
% % 
 plot_locations = 1;
% 
 hold all
 plot([x1 x2],[y1 y2],'ks-')
% 
% 
% % Plot data
 figure
plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)

%% Plot vertical radar data (data extract from xy figure
% Written on 2015-04-01 to regenerate plots
%  (1) Run plotter2 and update the date
%  (2) Choose the radar angle (if you need 2d plot, for some reason)
%  (3) Choose user input correctly
%  (4) If auto = 1, just open the previous Distance - altitude plot and run
%  this program
%  (5) If auto = 0, draw the 2D plot from plotter2, indicate x1,x2,theta
%  values. and run this.

% 
% % User inputs
% plot_locations = 1;
% auto = 1;
% x1 = -16.8;
% y1 = 0;
% theta = 90;
% D = 9; % Plot distance in  km
% save_figure = 1; % Make this 1 if you want to save figures automatically.
% save_dn = 'C:\Users\Sumedhe\Desktop\NadeeRadar\';
% 
% 
% if auto
%     % Get plotter2 data
%     try h=guidata(findall(0,'Tag','plotter2'));
%     catch; disp('Run plotter2 first!'); return;
%     end
%     
%     g = h.g;
%     
%     % Get the title of the plot
%     %figh = gcf;
%     axh = gca;
%     tit = get(get(axh,'title'),'String');
%     
%     year = tit(2,1:4);
%     
%     if strcmp(year,'2011')
%         stit = tit(2,11:end);
%         % Let's update the time and return
%         n1 = strfind(stit,'UT:');
%         n2 = strfind(stit,'-');        
%         t1 = hhmmss2sec(stit(n1+3:n2-1),0);
%         t2 = hhmmss2sec(stit(n2+1:end),0);
%         g.t1 = t1;
%         g.t2 = t2;
%         g2plotter(g)     
%         return
%         
%     else
%         
%         %rarar_fn = tit(1,1:end-11);
%         n = strfind(tit(2,:),'--');
%         t1 = hhmmss2sec(tit(2,1:n-1),0);
%         n1 = strfind(tit(2,:),'UT');
%         t2 = hhmmss2sec(tit(2,n+2:n1-1),0);
%         g.t1 = t1;
%         g.t2 = t2;
%         g2plotter(g)
%         h=guidata(findall(0,'Tag','plotter2'));
%         g = h.g;        
%     end
% 
%     
%     % Flash number
%     fileN = get(gcf,'FileName');
%     
%     if isempty(fileN)
%         flashN = inputdlg('Enter flash number:');
%         
%         if isempty(flashN) || strcmp(flashN{:},'')
%             disp('no flash number entered')
%             return
%         else
%             flashN = str2double(flashN);
%         end
%         
%     else        
%         % flashN = fileN(end-17:end-16); % flash no less than 10
%         flashN = fileN(end-18:end-16); % Greater than 10
%         flashN(strfind(flashN,'-'))='';
%         flashN = str2double(flashN);
%     end
%     
% 
%     % Get x1,x2,y1,y2
%     xlab = get(get(axh,'xlabel'),'String');
%     n1 = strfind(xlab,'(');
%     n3 = strfind(xlab,')');
%     n2 = n1(2)-1 + strfind(xlab(n1(2):n3(2)),',');
%     n4 = strfind(xlab,'=');
%     theta = str2double(xlab(n4+1:end-1));
%     x1 = str2double(xlab(n1(2)+1:n2-1));
%     y1 = str2double(xlab(n2+1:n3(2)-1));
%     xL = xlim(axh);
%     x2 = x1 + xL(2)*cosd(theta);
%     y2 = y1 + xL(2)*sind(theta);
%     
%     h.g = g;
%     
%     % Plot the xy figure
%     % Set the radar file name
%     [~ , radarFullFiles] = getRadarFiles(h);
%     sen_set = h.sen_set;
%     sen_set.radarOn = 1;
%     sen_set.radarFn = radarFullFiles{1};
%     h.sen_set = sen_set;
%     ldar_plot_execute(h);
%     hold all
%     plot([x1,x2],[y1,y2],'ks-','linewidth',1)
%     
%     
%     % COrrect xlim and ylim
%     xlim([min([x1-10,x1+10,x2-10,x2+10]),max([x1-10,x1+10,x2-10,x2+10])])
%     ylim([min([y1-10,y1+10,y2-10,y2+10]),max([y1-10,y1+10,y2-10,y2+10])])
%     
%     if save_figure
%         fn = sprintf('%s%s-%s-%s-flash%3.3i-RADAR-xy',save_dn,g.YYYY{:},g.MM{:},g.DD{:},flashN);
%         saveas(gcf,[fn '.fig'])
%         saveas(gcf,[fn '.png'])
%     end
%         
%     
% else
%     x2 = x1 + D*cosd(theta);
%     y2 = y1 + D*sind(theta);
% end
% 
% % Plot vertical radar plane
% a = guidata(gcf);
% fn = a.fileName;
% plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations);
% figH = gcf;
% 
% if save_figure
%     fn = sprintf('%s%s-%s-%s-flash%3.3i-RADAR-3D',save_dn,g.YYYY{:},g.MM{:},g.DD{:},flashN)
%     saveas(gcf,[fn '.fig'])
%     saveas(gcf,[fn '.png'])
% end
%     
%     % Plot smoothed 2d radar plane
% figure(figH); radarTools(4) % For 2D smooth curve
% if save_figure
%     fn = sprintf('%s%s-%s-%s-flash%3.3i-RADAR-VS2D-smooth',save_dn,g.YYYY{:},g.MM{:},g.DD{:},flashN);
%     saveas(gcf,[fn '.fig'])
%     saveas(gcf,[fn '.png'])
% end
% 
% figure(figH); radarTools(2) % For 2D real curve
% if save_figure
%     fn = sprintf('%s%s-%s-%s-flash%3.3i-RADAR-VS2D-real',save_dn,g.YYYY{:},g.MM{:},g.DD{:},flashN);
%     saveas(gcf,[fn '.fig'])
%     saveas(gcf,[fn '.png'])
% end
% 
% 
% 







%% Do not comment below this line
 function r = convert2D(vrbd,x,y)

x1 = vrbd(1);
x2 = vrbd(3);
y1 = vrbd(2);
y2 = vrbd(4);

% slope and intercept of the line (need later)
m = (y2 - y1)/(x2 - x1);
b = (x2*y1 - x1*y2)/(x2 - x1);

% Snap the location coordinates to line
if m == inf
    y = y;
    x = x1;
else
    x = (m*y+x-m*b)/(m*m + 1);
    y = (m*m*y+m*x+b)/(m*m + 1);
end

% Make them 2D
r = sqrt((x - x1).^2 + (y - y1).^2);

function setup_fig(rfn,xl,yl)


sns_x = [ -0.7524, 0.3372, -0.3149, -1.1658, -0.2701, -2.7640, -6.0091, 0.1825, -5.7394, -2.0637]*10;
sns_y = [  1.6555, 0.4446, -0.6838, -1.7020,  0.2631,  4.9254, -3.3983, -5.3008, 1.1923, 0.1569]*10;
sns_n = { 'K02', 'K14',    'K24',    'BCC',    'K17',    'EDW',    'STC',    'FLT',    'OVD',    'FFI'};


% setup

xlabel('East (km)')
ylabel('North (km)')

plot(sns_x,sns_y,'kp','markerfacecolor','k')
text(sns_x,sns_y,sns_n)
florida_map
box on

[rdn, rfn ,extn ] = fileparts(rfn);
title(sprintf('%s%s', rfn,extn),...
    'Interpreter','none')


xlim(xl)
ylim(yl)


daspect([1 1 1])



function ldar_plot_execute(handles)

if strcmp(handles.sen_set.autoBaseFolder,'on')
    handles = set_auto_b_folder(handles);
end

g = handles.g;
sen_set = handles.sen_set;

% set the base folder


if g.mm < 30
    ext=0;
else
    ext=30;
end

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);


linet_fn=sprintf('%s/LINET/%s/%s/%s/linet_%s%s%s_%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

% nldn_fn=sprintf('%s/LINET2/%s/%s/%s/linet_%s%s%s_%2.2i%2.2i.txt',...
%     sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
%     g.MM{:},g.DD{:},g.hh,ext);

nldn_fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

%Old PBFA fn
pbfaO_fn=sprintf('%s/PBFA_old/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

nldn2_fn=sprintf('%s/NLDN2/%s/%s/NLDN2_%s%s%s.txt',...
        sen_set.loc_dir,g.YYYY{:},g.MM{:},g.YYYY{:},...
        g.MM{:},g.DD{:});
    
% Graph selection from the user
num=sen_set.ldar_graph_type;

% Do we need ldar time shift to be included?
if sen_set.ldar_tshiftOn==1
    x0=sen_set.x(sen_set.ldar_tshift_sn)-sen_set.x0;
    y0=sen_set.y(sen_set.ldar_tshift_sn)-sen_set.y0;
    z0=sen_set.z(sen_set.ldar_tshift_sn)-sen_set.z0;
else
    x0=0;
    y0=0;
    z0=0;
end

switch num
    case 1
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,sen_set,g)
    case 2
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,2,sen_set,g)
    case 3
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,3,sen_set,g)
    case 4
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,4,sen_set,g)
    case 5
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[1,0,0,0,0],sen_set,g);
    case 6
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[0,1,0,0,0],sen_set);
    case 7
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[0,0,1,0,0],sen_set,g);
    case 8
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[0,0,0,1,0],sen_set,g);
    case 9
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[0,0,0,0,1],sen_set,g);
    case 10
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,0,[0,1,0,0,0],sen_set,g);
    case 11
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,5,sen_set,g)
    case 12
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,6,sen_set,g)
    case 13
        LDAR_Hist1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,sen_set)
    case 14
        LDAR_frequency(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,sen_set)
        
    otherwise
        disp xc
        
end

function [radarFiles , radarFullFiles] = getRadarFiles(handles)

g = handles.g;
sen_set = handles.sen_set;

% Radar file folder
rdn = sprintf('%s/netCDF/%s/%s/%s/%s/%2.2i/',...
    sen_set.base_dir,sen_set.radarStations{sen_set.radarStationID},...
    g.YYYY{:},g.MM{:},g.DD{:},g.hh);

files = dir(rdn);

L = length(files);
times = nan(1,L-2);

radarFiles = {};
radarFullFiles = {};

if L > 2    
    for i = 3:L
        fn =  files(i).name;
        times(i-2) = str2double(fn(14:15))*3600+str2double(fn(16:17))*60+str2double(fn(18:19));
        radarFiles = [radarFiles; fn];
        radarFullFiles = [radarFullFiles; [rdn fn]];
    end
end

% if within the first 10 mins of the hour, let's load files from the
% previous hour
if g.mm < 10
    dNum1 = datenum(str2double(g.YYYY{:}),str2double(g.MM{:}),str2double(g.DD{:}),g.hh,0,0);
    
    % substract an hour
    dNum2 = dNum1 - 1/24;
    
    dvec = datevec(dNum2);
    
    % Radar file folder
    rdn = sprintf('%s/netCDF/%4.4i/%2.2i/%2.2i/%2.2i/',...
        sen_set.base_dir,dvec(1),dvec(2),dvec(3),dvec(4));
    
    files = dir(rdn);
    
    L = length(files);
    
    % We just need the last file of this folder
    if L > 2
        fn =  files(L).name;
        times =[times, str2double(fn(14:15))*3600+str2double(fn(16:17))*60+str2double(fn(18:19))];
        radarFiles = [radarFiles; fn];
        radarFullFiles = [radarFullFiles; [rdn fn]];
    end
end

% if within the first 10 mins of the hour, let's load files from the
% previous hour
if g.mm > 50
    dNum1 = datenum(str2double(g.YYYY{:}),str2double(g.MM{:}),str2double(g.DD{:}),g.hh,0,0);
    
    % add an hour
    dNum2 = dNum1 + 1/24;
    
    dvec = datevec(dNum2);
    
    % Radar file folder
    rdn = sprintf('%s/netCDF/%4.4i/%2.2i/%2.2i/%2.2i/',...
        sen_set.base_dir,dvec(1),dvec(2),dvec(3),dvec(4));
    
    files = dir(rdn);
    
    L = length(files);
    
    % We just need the first file of this folder
    if L > 2
        fn =  files(3).name;
        times =[times, str2double(fn(14:15))*3600+str2double(fn(16:17))*60+str2double(fn(18:19))];
        radarFiles = [radarFiles; fn];
        radarFullFiles = [radarFullFiles; [rdn fn]];
    end
end

if ~isempty(times)
    dts = abs(times - g.t1);
    [dts I] = sort(dts);
    radarFiles = radarFiles(I);
    radarFullFiles = radarFullFiles(I);
else
    radarFiles = '';
    radarFullFiles = '';
end