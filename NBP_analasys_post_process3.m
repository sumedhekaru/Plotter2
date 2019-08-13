function NBP_analasys_post_process3

clc
%% User inputs

% file name for data
fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814-NBP_info2.xlsx';
%
% % base folder for saving plots
% bf = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\Results\';
%
% % Save plots?
savePlots = 0;

sns_x = [ -0.7524, 0.3372, -0.3149, -1.1658, -0.2701, -2.7640, -6.0091, 0.1825, -5.7394, -2.0637]*10;
sns_y = [  1.6555, 0.4446, -0.6838, -1.7020,  0.2631,  4.9254, -3.3983, -5.3008, 1.1923, 0.1569]*10;
sns_n = { 'K02', 'K14',    'K24',    'BCC',    'K17',    'EDW',    'STC',    'FLT',    'OVD',    'FFI'};

%%
sheet = 1;
xlRange = 'A3:AD307';

[ndata, txtdata, data] = xlsread(fn, sheet, xlRange);

NBPinds = ndata(:,1);
pbfat = ndata(:,10);
pbfax = ndata(:,11);
pbfay = ndata(:,12);
pbfaz = ndata(:,13);




%% Open the file with new types (6 types as of Nov 05,2014)

%   0 - go to next (Can't determine the type) (type 0)
%   1 - Catogory 1 - Clean bipolar pulses
%   2 - Catogory 2 - Bipolar w extra pulses after overshoot peak
%   3 - Catogory 3 - extra pulses befor and afer overshoot
%   4 - Catogory 4 - Bipolar with extra pulse between first peak and zero cross
%   5 - Catogory 5 - extra pulses before leading peak
%   6 - Catogory 6 - Strainge type of NBP

%fID = fopen('C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\Types-20141030\20110814-type_info.txt','r');
%typeD = textscan(fID,'%f %f');
%type = typeD{2};
%fclose(fID);


%% Load NBP times to plotter2
% loadInd = 248;
% pbfat(loadInd)
% if ~isnan(pbfat(loadInd))
%     figh = figure;
%     xlim([pbfat(loadInd)-500e-6 pbfat(loadInd)+500e-6])
%     plotter_tools(3);
%     delete(figh)
%     %peak_modifier
% else
%     disp('There is no NBP for this point')
% end
% return


%% Current vs time (MTLL)
% t = 0e-6:0.01e-6:120e-6;
% 
% t1 = 10e-6;
% t2 = 50e-6;
% t3 = 75e-6;
% t2p = 20e-6;
% alpha =15/t2;
% k = (t2-t1)/t1;
% A = 1;
% 
% ii = t;
% 
% for i = 1:length(t)
%     if t(i) <= t1
%         ii(i) = A*exp(-(alpha*(t(i)-t1))^2);
%     elseif t(i) <= t1 + t2p
%         ii(i) = A*exp(-(alpha*(t(i)-t1)/k)^2);
%     elseif t(i) <= (t1 + t2p + t3)
%         ii(i) = -A*exp(-(alpha*(t2p/k))^2)*(t(i)-t1-t2p-t3)/t3;
%     else
%         ii(i) = 0;
%     end
% end
% 
% figure
% hold all
% plot(t,ii,'linewidth',2)
% ylim([0 1.1]);
% xlim([0 120e-6]);
% box on
% 
% plot([0 t1*2],[1 1],'k--')
% plot([t1 t1],[0  1],'k--')
% plot([t1 t1]+t2p,[0  2*A*exp(-(alpha*(t2p/k))^2)],'k--')
% plot([t1 t1]+t2p+t3,[0  2*A*exp(-(alpha*(t2p/k))^2)],'k--')
% set(gca,'xtick',[t1, t1+t2p ,t1+t2p+t3],'ytick',[1],'xticklabel',{'t_1','t''_2','t_3'},'yticklabel','A')
% tools2fig
%     
% 
% return

%% NBP generate vertically stacked figures time shifted to PBFA
% %Usage : NBP_plot_all_verticle(plot_range); where plot_range = 1:101
% 
% index = 286;
% 
% % Spacing between plots. (could be one value, could be set of values)
% aa.spacing = [30 30 50 50 40 30 40];
% aa.spacing = 10;
% %aa.spacing = [10 10 10 10 10 30 30]
% 
% % pre and post times (in microseconds)
% aa.t1t2 = [500 10000];
% %aa.t1t2 = [150000 150000];
% %aa.t1t2 = [25 125];`
% %aa.t1t2 = [200 200];
% 
% % Time shifts in us
% aa.dts = [0 0 0 0 0 0 0 0 ];
% 
% % x location of text on each plot
% aa.textx = 8000;
% aa.textx = sum(aa.t1t2)*0.95;
% 
% 
% NBP_plot_all_verticle(index:index,aa);
% tools2fig
% return


%% NBP total time info excel reader

% xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\TotalTime\NBP_total_durations-2015-07-29.xlsx';
% 
% sheet = 1;
% xlRange = 'A2:F67';
% 
% [ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);
% 
% indx = ndata(:,1);
% D = ndata(:,3);
% dE = ndata(:,4);
% dt = ndata(:,5);
% P = ndata(:,6);
% 
% % Total number of samples
% fprintf('Total number of samples = %i\n',length(unique(indx)));
% 
% % Total number of data obtained
% ind1 = find(~isnan(dE));
% ind2 = indx(ind1);
% ind3 = unique(ind2);
% length(ind3);
% fprintf('Total number of NBPs obtained for = %i \n',length(ind3));
% 
% 
% % Total number data points
% fprintf('Total number of data points = %i\n',length(ind1))
% 
% 
% 
% % Histograms
% figure
% 
% subplot(2,2,4)
% hist(dt,2.5:5:35)
% xlabel('Total time dt (ms)')
% ylabel('Count')
% legend([sprintf('(%0.1f',nanmean(dt)) '\pm' sprintf('%0.1f) ms',nanstd(dt))])
% set(gca,'xtick',0:5:35)
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% 
% 
% %figure
% subplot(2,2,2)
% hist(dE,-55:10:0)
% xlabel('Total-electrostatic change \DeltaE_2 (ms)')
% ylabel('Count')
% legend([sprintf('(%0.1f',nanmean(dE)) '\pm' sprintf('%0.1f) V/m',nanstd(dE))],'location','northwest')
% set(gca,'xtick',-60:10:0)
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% 
% 
% % figure
% % hist(D,1.25:2.5:15)
% % xlabel('Total time dt (ms)')
% % ylabel('Count')
% 
% fprintf('total E Max = %0.1f\n',nanmax(dE))
% fprintf('total E Min = %0.1f\n',nanmin(dE))
% 
% fprintf('Total t Max = %0.1f\n',nanmax(dt))
% fprintf('Total t Min = %0.1f\n',nanmin(dt))


%% Statistics of fast electrostatic data

% xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\Fast_E_static_data\Fast_E_static_data.xlsx';
% 
% sheet = 1;
% xlRange = 'A2:E82';
% 
% [ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);
% 
% indx = ndata(:,1);
% D = ndata(:,3);
% dE = ndata(:,4);
% dE2 = ndata(:,5);
% %dt = ndata(:,5);
% %P = ndata(:,6);
% 
% % Total number of samples
% fprintf('Total number of samples = %i\n',length(unique(indx)));
% 
% % Total number of data obtained
% ind1 = find(~isnan(dE));
% ind2 = indx(ind1);
% ind3 = unique(ind2);
% length(ind3);
% fprintf('Total number of NBPs obtained for = %i \n',length(ind3));
% 
% 
% % Total number data points
% fprintf('Total number of data points = %i\n',length(ind1))
% 
% %figure
% subplot(2,2,1)
% hist(dE,-18.75:2.5:0)
% xlabel('Fast-electrostatic change \DeltaE_1 (ms)')
% ylabel('Count')
% legend([sprintf('(%0.1f',nanmean(dE)) '\pm' sprintf('%0.1f) V/m',nanstd(dE))],'location','northwest')
% set(gca,'xtick',-20:2.5:0)
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% 
% % Statistics of fast electro static changes
% fprintf('Fast E Max = %0.1f\n',nanmax(dE))
% fprintf('Fast E Min = %0.1f\n',nanmin(dE))
% 
% 
% % Slow electrostatic change
% dE3 = dE2 - dE;
% 
% % Statistics of fast electro static changes
% fprintf('Slow E Max = %0.1f\n',nanmax(dE3))
% fprintf('Slow E Min = %0.1f\n',nanmin(dE3))
% 
% 
% 
% %figure
% subplot(2,2,3)
% hist(dE3,-42.5:5:0)
% xlabel('Slow-electrostatic change \DeltaE_3 (ms)')
% ylabel('Count')
% legend([sprintf('(%0.1f',nanmean(dE3)) '\pm' sprintf('%0.1f) V/m',nanstd(dE3))],'location','northwest')
% set(gca,'xtick',-45:5:0)
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% 
% 
% 
% return

%% Scatter plots of dE1, dE2, dE3 with R
% xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\TotalTime\NBP_total_durations-2015-07-29.xlsx';
% 
% sheet = 1;
% xlRange = 'A2:F67';
% 
% [ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);
% 
% indx1 = ndata(:,1);
% D1 = ndata(:,3);
% dE = ndata(:,4);
% dt = ndata(:,5);
% P1 = ndata(:,6);
% 
% 
% xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\Fast_E_static_data\Fast_E_static_data_b.xlsx';
% 
% sheet = 1;
% xlRange = 'A2:E82';
% 
% [ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);
% 
% indx2 = ndata(:,1);
% D2 = ndata(:,3);
% dE_fast = ndata(:,4);
% dE_total = ndata(:,5);
% dE_slow = dE_total - dE_fast;
% 
% 
% % Get altitudes
% xlFn = 'C:\Users\sumedhe\Desktop\NBP-timing-JGR-2015\20110814-NBP_info3.xlsx';
% sheet = 1;
% xlRange = 'A2:M306';
% 
% ndataz = xlsread(xlFn, sheet, xlRange);
% indzs = ndataz(:,1);
% zs = ndataz(:,13)/1000;
% 
% 
% Z1 = zeros(size(D1));
% 
% for i = 1:length(D1)
%     ind = find(indzs == indx1(i));
%     Z1(i) = zs(ind);
% end
% 
% Z2 = zeros(size(D2));
% 
% for i = 1:length(D2)
%     ind = find(indzs == indx2(i));
%     Z2(i) = zs(ind);
% end
% 
% R1 = sqrt(D1.^2 + Z1.^2);
% R2 = sqrt(D2.^2 + Z2.^2);
% 
% figure
% subplot(2,2,1)
% plot(R2,dE_fast,'ro')
% box on
% xlabel('Range (km)')
% ylabel('\Delta E_{fast} (V/m)')
% 
% subplot(2,2,2)
% plot(R2,dE_slow,'ro')
% box on
% xlabel('Range (km)')
% ylabel('\Delta E_{slow} (V/m)')
%    
% 
% subplot(2,2,3)
% plot(R2,dE_total,'ro')
% box on
% xlabel('Range (km)')
% ylabel('\Delta E_{total} (V/m)')
% 
% 
% subplot(2,2,4)
% plot(R1,dt,'ro')
% box on
% xlabel('Range (km)')
% ylabel('\Delta t (ms)')
% 
% 
% tools2fig
%    
% 
% return



%% Correlation figures of dE1, dE2, dE3, D
% xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\TotalTime\NBP_total_durations-2015-07-29.xlsx';
% 
% sheet = 1;
% xlRange = 'A2:F67';
% 
% [ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);
% 
% indx1 = ndata(:,1);
% D1 = ndata(:,3);
% dE = ndata(:,4);
% dt = ndata(:,5);
% P1 = ndata(:,6);
% 
% % Get altitudes
% xlFn = 'C:\Users\sumedhe\Desktop\NBP-timing-JGR-2015\20110814-NBP_info3.xlsx';
% sheet = 1;
% xlRange = 'A2:M306';
% 
% ndataz = xlsread(xlFn, sheet, xlRange);
% indzs = ndataz(:,1);
% zs = ndataz(:,13)/1000;
% 
% R1 = D1;
% 
% for i = 1:length(R1)
%     ind = indzs == indx1(i);
%     R1(i) = sqrt(D1(i)^2 + zs(ind)^2);
% end
% 
% 
% 
% 
% xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\Fast_E_static_data\Fast_E_static_data.xlsx';
% 
% sheet = 1;
% xlRange = 'A2:E82';
% 
% [ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);
% 
% indx2 = ndata(:,1);
% D2 = ndata(:,3);
% dE1 = ndata(:,4);
% dE2 = ndata(:,5);
% dE3 = dE2 - dE1;
% 
% R2 = D2;
% 
% for i = 1:length(R2)
%     ind = indzs == indx2(i);
%     R2(i) = sqrt(D2(i)^2 + zs(ind)^2);
% end
% 
% 
% 
% figure
% tools2fig
% subplot(2,2,1)
% inds = find(~isnan(dE1));
% y = dE1(inds);
% x = R2(inds);
% plot(x,y,'ro')
% xlabel('Range (km)')
% ylabel('\DeltaE_{fast} (V/m)')
% ylim([-60 5])
% 
% %figure
% subplot(2,2,2)
% inds = find(~isnan(dE3));
% y = dE3(inds);
% x = R2(inds);
% plot(x,y,'ro')
% xlabel('Range (km)')
% ylabel('\DeltaE_{slow} (V/m)')
% ylim([-60 5])
% 
% subplot(2,2,3)
% inds = find(~isnan(dE2))
% y = dE2(inds);
% x = R2(inds);
% plot(x,y,'ro')
% xlabel('Range (km)')
% ylabel('\DeltaE_{total} (V/m)')
% ylim([-60 5])
% 
% subplot(2,2,4)
% inds = find(~isnan(dt))
% y = dt(inds);
% x = R2(inds);
% plot(x,y,'ro')
% xlabel('Range (km)')
% ylabel('\Deltat (ms)')
% 
% 
% return
% 
% close all
% 
% figure
% tools2fig
% subplot(3,3,1)
% inds = find(~isnan(dE1));
% y = dE1(inds);
% x = D2(inds);
% plot(x,y,'ro','markerfacecolor','r')
% [fittedX,fittedY,P,S] = best_fit(x,y);
% xlabel('Range (km)')
% ylabel('Fast-electrostatic change \DeltaE_1 (V/m)')
% hold all
% plot(fittedX,fittedY,'LineWidth',1);
% %legend({'data', sprintf('Pearson %0.2f\nSpearman %0.2f',P,S)},'location','southeast')
% title(sprintf('r = %0.2f    r_s = %0.2f',P,S))
% 
% 
% %figure
% subplot(3,3,2)
% inds = find(~isnan(dE2));
% y = dE2(inds);
% x = D2(inds);
% plot(x,y,'ro','markerfacecolor','r')
% [fittedX,fittedY,P,S] = best_fit(x,y);
% xlabel('Range (km)')
% ylabel('Total-electrostatic change \DeltaE_2 (V/m)')
% hold all
% plot(fittedX,fittedY,'LineWidth',1);
% %legend({'data', sprintf('Pearson %0.2f\nSpearman %0.2f',P,S)},'location','southeast')
% title(sprintf('r = %0.2f    r_s = %0.2f',P,S))
% 
% %figure
% subplot(3,3,3)
% inds = find(~isnan(dE3));
% y = dE3(inds);
% x = D2(inds);
% plot(x,y,'ro','markerfacecolor','r')
% [fittedX,fittedY,P,S] = best_fit(x,y);
% xlabel('Range (km)')
% ylabel('Slow-electrostatic change \DeltaE_3 (V/m)')
% hold all
% plot(fittedX,fittedY,'LineWidth',1);
% %legend({'data', sprintf('Pearson %0.2f\nSpearman %0.2f',P,S)},'location','southeast')
% title(sprintf('r = %0.2f    r_s = %0.2f',P,S))
% 
% %figure
% subplot(3,3,4)
% inds = find(~isnan(dt));
% y = dt(inds);
% x = D2(inds);
% plot(x,y,'ro','markerfacecolor','r')
% [fittedX,fittedY,P,S] = best_fit(x,y);
% xlabel('Range (km)')
% ylabel('Total-time dt (ms)')
% hold all
% plot(fittedX,fittedY,'LineWidth',1);
% %legend({'data', sprintf('Pearson %0.2f\nSpearman %0.2f',P,S)},'location','southeast')
% title(sprintf('r = %0.2f    r_s = %0.2f',P,S))
% 
% 
% subplot(3,3,5)
% inds = find(~isnan(dt));
% y = dt(inds);
% x = dE(inds);
% plot(x,y,'ro','markerfacecolor','r')
% [fittedX,fittedY,P,S] = best_fit(x,y);
% xlabel('Total-electrostatic change \DeltaE_2(V/m)')
% ylabel('Total-time dt (ms)')
% hold all
% plot(fittedX,fittedY,'LineWidth',1);
% %legend({'data', sprintf('Pearson %0.2f\nSpearman %0.2f',P,S)},'location','southeast')
% title(sprintf('r = %0.2f    r_s = %0.2f',P,S))
% 
% 
% %figure
% subplot(3,3,7)
% inds = find(~isnan(dE2));
% y = dE2(inds);
% x = dE1(inds);
% plot(x,y,'ro','markerfacecolor','r')
% [fittedX,fittedY,P,S] = best_fit(x,y);
% xlabel('Fast-electrostatic change \DeltaE_1(V/m)')
% ylabel('Total-electrostatic change \DeltaE_2(V/m)')
% hold all
% plot(fittedX,fittedY,'LineWidth',1);
% %legend({'data', sprintf('Pearson %0.2f\nSpearman %0.2f',P,S)},'location','southeast')
% title(sprintf('r = %0.2f    r_s = %0.2f',P,S))
% 
% 
% 
% %figure
% subplot(3,3,8)
% inds = find(~isnan(dE3));
% y = dE3(inds);
% x = dE1(inds);
% plot(x,y,'ro','markerfacecolor','r')
% [fittedX,fittedY,P,S] = best_fit(x,y);
% xlabel('Fast-electrostatic change \DeltaE_1(V/m)')
% ylabel('Slow-electrostatic change \DeltaE_2(V/m)')
% hold all
% plot(fittedX,fittedY,'LineWidth',1);
% %legend({'data', sprintf('Pearson %0.2f\nSpearman %0.2f',P,S)},'location','southeast')
% title(sprintf('r = %0.2f    r_s = %0.2f',P,S))
% 
% return
% 
% figure


%% Plot 3-D charge moments on RADAR plot

% indexes of NBPs which has 3d charge moments
inds = [92 94 206 158 195 177 174 196 192 194];

% Read 3-D charge moment data
xlFn = 'C:\Users\sumedhe\Desktop\NBP-timing-JGR-2015\data\Fast_charge_moments.xlsx';
sheet = 1;
xlRange = 'A3:I81';
[ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);

% get locatiosn of NBPs




for i = 1:length(inds)
    load_data2plotter(inds(i))
    % get plotter2 data
    h=guidata(findall(0,'Tag','plotter2'));
    sen_set = h.sen_set;
    g = h.g;
    
    % charge moments
    ind = find(ndata(:,1) == inds(i));
    ind = ind(1);
    px = ndata(ind,2);
    py = ndata(ind,3);
    pz = ndata(ind,4);
    pr = sqrt(px^2+py^2);
    
    % charge moment horizontal angle
    fprintf('Px = %0.1f    Py = %0.1f     Pz = %0.1f\n',px,py,pz)
    p_ang = atan2d(py,px);
           
    % Get locations of NBPs
    ind = find(NBPinds == inds(i));
    x0 = pbfax(ind)/1000;
    y0 = pbfay(ind)/1000;
    z0 = pbfaz(ind)/1000;
    
    [radarFiles , radarFullFiles] = getRadarFiles(g,sen_set);
    
    % Genarate horizontal radar plot
    % Generate horizontal plot
    sen_set.radarFn = radarFullFiles{1};
    sen_set.radEleAngInd = 5;
    figure(10); cla; 
    radar_plot(sen_set,g)
    xlim([x0-40 x0+40])
    ylim([y0-40 y0+40])
    florida_map
    
    % Extra menu
    f = uimenu('Label','Plotter');
    uimenu(f,'Label','Vertical Radar Plane','Callback','plotter_tools(27)');
    uimenu(f,'Label','Radar altitude curser','Callback','radarTools(1)');
    
    % Plot NBP on the horizontal radar
    plot(x0,y0,'ko','markersize',5)    
    daspect([1 1 1])
    box on
    
    % North south plot
    dR = 20;
    
    angs = [p_ang p_ang+90];
    
    for j = 1:2
        ang = angs(j); % in degrees
        
        x1 = x0 - dR*cosd(ang);
        x2 = x0 + dR*cosd(ang);
        y1 = y0 - dR*sind(ang);
        y2 = y0 + dR*sind(ang);
        
        figure(10)
        plot([x1 x2],[y1 y2],'-sk','markersize',5)
        ax1 = gca;
        
        % Plot vertical radar data
        plot_vert_radar_plane2(radarFullFiles{1},x1,y1,x2,y2,0)
        f2 = gcf;
        %radarTools(4) % for smoothed radar data
        radarTools(2) % for real 2-d radar data
        f3 = gcf;
        delete(f2);
        
        r = convert2D([x1 y1 x2 y2],x0,y0);
        plot(r,z0,'ko');
        text(r+0.5,z0,num2str(inds(i)),'FontWeight','bold')
        
        if j == 1
            quiver(r,z0,pr,pz,0.005,'filled','color','k', 'MaxHeadSize',1)
            ax2 = gca;
            xl2 = get(get(gca,'xlabel'),'String');
        else
            quiver(r,z0,0,pz,0.005,'filled','color','k', 'MaxHeadSize',1)
            ax3 = gca;
            xl3 = get(get(gca,'xlabel'),'String');
        end
    end
    
    % make a final figure
    figure
    sbh = subplot(2,3,[1,4]);
    copyobj(get(ax1,'children'),sbh)
    daspect([1 1 1])
    set(gcf,'renderer','painters')
    xlim([x0-40 x0+40])
    ylim([y0-40 y0+40])
    title(radarFiles{1},'interpreter','none')
    box on
    
    
    sbh = subplot(2,3,2:3);
    copyobj(get(ax2,'children'),sbh)
    daspect([1 1 1])
    set(gcf,'renderer','painters')
    xlabel(xl2)
    ylabel('Altitude (km)')
    box on
    title(sprintf('Px = %0.1f    Py = %0.1f     Pz = %0.1f\n',px,py,pz))
    
    sbh = subplot(2,3,5:6);
    copyobj(get(ax3,'children'),sbh)    
    daspect([1 1 1])
    set(gcf,'renderer','painters')
    box on
    
    xlabel(xl3)
    ylabel('Altitude (km)')
    tools2fig
    
    % save figure
    bf = 'C:\Users\sumedhe\Desktop/NBP-timing-JGR-2015\data\charge_moments_on_radar\';
    svfn = sprintf('%s%3.3i-real',bf,inds(i));
    saveas(gcf,[svfn '.fig']);
    save_full_screen_png([svfn '.png'],150)
    close all
   
    
end


%% Estimate charge moments
% %xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\TotalTime\NBP_total_durations-2015-07-29.xlsx';
% %sheet = 1;
% %xlRange = 'A2:E65';
% 
% xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\Fast_E_static_data\Fast_E_static_data_b.xlsx';
% sheet = 1;
% xlRange = 'A2:E82';
% 
% 
% [ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);
% 
% indx = ndata(:,1);
% sIDs = ndata(:,2);
% D = ndata(:,3);
% dE = ndata(:,5);
% 
% for i = 1:length(indx) 
%     
%     ind = indx(i);
%     
%     if isnan(dE(i))
%         fprintf('%i %10.1f %10.2e\n',ind,NaN,NaN)
%     else
%         
%         
%         ind0 = find(ind == indx);
%         
%         dEms = [];
%         IDs = [];
%         
%         for i=1:length(ind0);
%             if ~isnan(dE(ind0(i)))
%                 dEms = [dEms dE(ind0(i))];
%                 IDs = [IDs sIDs(ind0(i))];
%             end
%         end
%         
%         x = pbfax(ind);
%         y = pbfay(ind);
%         z = pbfaz(ind);
%         
%         
%         
%         [P, px, py, pz, ki_sqrd] = point_dipole([x y z],IDs,dEms,0);
%         %fprintf('%i %10.1f %10.1f %10.1f %10.1f %10.2e\n',ind,P, px, py, pz, ki_sqrd)
%         
%         [P2, ki_sqrd2] = point_dipole_1([x y z],IDs,dEms,0);
%         %fprintf('%i %10.1f %10.2e\n',ind,P2,ki_sqrd2)
%         
%         fprintf('%i %10.1f %10.1f %10.1f %10.1f %10.2e %10.1f %10.2e\n',ind,px, py, pz, P, ki_sqrd, P2, ki_sqrd2)
%     end
%        
% end
% 
% return

%% Estimate charge moments for Fast change
% %xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\TotalTime\NBP_total_durations-2015-07-29.xlsx';
% %sheet = 1;
% %xlRange = 'A2:E65';
% 
% xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\Fast_E_static_data\Fast_E_static_data_b.xlsx';
% sheet = 1;
% xlRange = 'A2:E82';
% 
% 
% [ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);
% 
% indx = ndata(:,1);
% sIDs = ndata(:,2);
% D = ndata(:,3);
% dE = ndata(:,4);
% 
% for i = 1:length(indx) 
%     
%     ind = indx(i);
%     
%     if isnan(dE(i))
%         fprintf('%i %10.1f %10.2e\n',ind,NaN,NaN)
%     else
%         
%         
%         ind0 = find(ind == indx);
%         
%         dEms = [];
%         IDs = [];
%         
%         for i=1:length(ind0);
%             if ~isnan(dE(ind0(i)))
%                 dEms = [dEms dE(ind0(i))];
%                 IDs = [IDs sIDs(ind0(i))];
%             end
%         end
%         
%         x = pbfax(ind);
%         y = pbfay(ind);
%         z = pbfaz(ind);
%         
%         
%         
%         [P, px, py, pz, ki_sqrd] = point_dipole([x y z],IDs,dEms,0);
%         %fprintf('%i %10.1f %10.1f %10.1f %10.1f %10.2e\n',ind,P, px, py, pz, ki_sqrd)
%         
%         [P2, ki_sqrd2] = point_dipole_1([x y z],IDs,dEms,0);
%         %fprintf('%i %10.1f %10.2e\n',ind,P2,ki_sqrd2)
%         
%         fprintf('%i %10.1f %10.1f %10.1f %10.1f %10.2e %10.1f %10.2e\n',ind,px, py, pz, P, ki_sqrd, P2, ki_sqrd2)
%     end
%        
% end
% 
% return

%% Calculate total charge moment (fast + slow)
%xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\TotalTime\NBP_total_durations-2015-07-29.xlsx';
%sheet = 1;
%xlRange = 'A2:E65';

% xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\Fast_E_static_data\Fast_E_static_data_b.xlsx';
% sheet = 1;
% xlRange = 'A2:E82';
% 
% 
% [ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);
% 
% indx = ndata(:,1);
% sIDs = ndata(:,2);
% D = ndata(:,3);
% dE_fast = ndata(:,4);
% dE_slow = ndata(:,5);
% dE = dE_fast + dE_slow;
% 
% for i = 1:length(indx) 
%     
%     ind = indx(i);
%     
%     if isnan(dE(i))
%         fprintf('%i %10.1f %10.2e\n',ind,NaN,NaN)
%     else
%         
%         
%         ind0 = find(ind == indx);
%         
%         dEms = [];
%         IDs = [];
%         
%         for i=1:length(ind0);
%             if ~isnan(dE(ind0(i)))
%                 dEms = [dEms dE(ind0(i))];
%                 IDs = [IDs sIDs(ind0(i))];
%             end
%         end
%         
%         x = pbfax(ind);
%         y = pbfay(ind);
%         z = pbfaz(ind);
%         
%         
%         
%         [P, px, py, pz, ki_sqrd] = point_dipole([x y z],IDs,dEms,0);
%         %fprintf('%i %10.1f %10.1f %10.1f %10.1f %10.2e\n',ind,P, px, py, pz, ki_sqrd)
%         
%         [P2, ki_sqrd2] = point_dipole_1([x y z],IDs,dEms,0);
%         %fprintf('%i %10.1f %10.2e\n',ind,P2,ki_sqrd2)
%         
%         fprintf('%i %10.1f %10.1f %10.1f %10.1f %10.2e %10.1f %10.2e\n',ind,px, py, pz, P, ki_sqrd, P2, ki_sqrd2)
%     end
%        
% end
% 
% return


%% Generate Charge moment plots

% % Read 3-D charge moment data
% inds = [92 94 206 158 195 177 174 196 192 194];
% xlFn = 'C:\Users\sumedhe\Desktop\NBP-timing-JGR-2015\data\Fast_charge_moments.xlsx';
% sheet = 1;
% xlRange = 'A3:I81';
% [ndata0, txtdata0, data0] = xlsread(xlFn, sheet, xlRange);
% 
% pz3d = [];
% for i = 1:length(inds)
% 
%     
%     % charge moments
%     ind = find(ndata0(:,1) == inds(i));
%     ind = ind(1);
%     %px = ndata(ind,2);
%     %py = ndata(ind,3);
%     pz3d = [pz3d ndata0(ind,4)];
%     %pr = sqrt(px^2+py^2);
% 
% end
% 
% 
% % 1D charge moment data
% xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\Fast_E_static_data\Fast_E_static_data_b.xlsx';
% sheet = 1;
% xlRange = 'A2:I82';
% 
% [ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);
% 
% 
% indx = ndata(:,1);
% P_fast = abs(ndata(:,6));
% P_total = abs(ndata(:,8));
% P_slow = P_total - P_fast;
% 
% 
% figure
% subplot(1,3,1)
% inds = find(~isnan(P_fast));
% x = P_fast(inds)/1000;
% [x,I] = unique(x);
% fprintf('Minimum P_fast = %0.1f\n',min(x))
% fprintf('Maximum P_fast = %0.1f\n',max(x))
% fprintf('Average P_fast = %0.1f\n',mean(x))
% 
% [N1 x1] = hist(x,(50:100:750)/1000);
% [N2 x2] = hist(abs(pz3d)/1000,(50:100:750)/1000);
% 
% bar(x1,[N1;N2]',1.0,'hist');
% 
% set(gca,'xtick',(0:100:800)/1000)
% xlim([0 800]/1000)
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% %title([sprintf('(%0.2f',nanmean(x)) '\pm' sprintf('%0.2f) C m',nanstd(x))])
% xlabel('Fast-charge moment (P_{fast}) (C km)')
% ylabel('Count')
% ax = gca;
% ax.XTickLabelRotation = 90;
% legend('P_z (1-D)', 'P_z (3-D)')
% 
% subplot(1,3,2)
% inds = find(~isnan(P_slow));
% x = P_slow(inds)/1000;
% 
% 
% [x,I] = unique(x);
% 
% fprintf('Minimum P_slow = %0.1f\n',min(x))
% fprintf('Maximum P_slow = %0.1f\n',max(x))
% fprintf('Average P_slow = %0.1f\n',mean(x))
% 
% [N1 x1] = hist(x,(100:200:2600)/1000);
% [N2 x2] = hist([1200 440 800]/1000,(100:200:2600)/1000);
% bar(x1,[N1;N2]',1.0)
% 
% set(gca,'xtick',(0:200:2600)/1000)
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% xlim([0 2200]/1000)
% %title([sprintf('(%0.0f',nanmean(x)) '\pm' sprintf('%0.0f) C m',nanstd(x))])
% xlabel('Slow-charge moment (P_{slow}) (Ck m)')
% ylabel('Count')
% ax = gca;
% ax.XTickLabelRotation = 90;
% legend('P_z (1-D)', 'P_z (3-D)')
% 
% 
% inds = find(~isnan(P_slow));
% x = P_fast(inds)/1000;
% y = P_slow(inds)/1000;
% [x,I] = unique(x);
% y = y(I);
% 
% %3d values
% x3d = [0.65 0.19 0.42];
% y3d = [1.20 0.44 0.80];
% 
% subplot(1,3,3)
% hold all
% plot(x,y,'ko','markerfacecolor','b')
% plot(x3d,y3d,'ko','markerfacecolor','y')
% legend('P_z (1-D)','P_z (3-D)')
% [fittedX,fittedY,P,S] = best_fit(x,y);
% cc = corrcoef(x,y);
% xlabel('Fast-charge moment (P_{fast}) (C km)')
% ylabel('Slow-charge moment (P_{slow}) (C km)')
% hold all
% plot(fittedX,fittedY,'LineWidth',1);
% title(sprintf('Co. Coeff. = %0.2f',cc(2)))
% box on
% ax = gca;
% ax.XTickLabelRotation = 90;
% 
% tools2fig
% 
% return
% 


%% E-Static Vs R^3 test

% xlFn = 'C:\Users\sumedhe\Desktop\NBP-timing-JGR-2015\20110814-NBP_info3.xlsx';
% sheet = 1;
% xlRange = 'A2:U306';
% 
% [ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);
% 
% indx1 = ndata(:,1);
% Zs = ndata(:,13);
% 
% xlFn = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\Fast_E_static_data\Fast_E_static_data.xlsx';
% 
% sheet = 1;
% xlRange = 'A2:E82';
% 
% [ndata, txtdata, data] = xlsread(xlFn, sheet, xlRange);
% 
% indx2 = ndata(:,1);
% D = ndata(:,3);
% dE1 = ndata(:,4);
% dE2 = ndata(:,5);
% dE3 = dE2 - dE1;
% 
% inds = indx2(1);
% is = 1;
% 
% % Best fit exponants
% n1s = [];
% n2s = [];
% 
% for i = 2:length(indx2)
%     if indx2(i) == inds(1)
%         inds = [inds, indx2(i)];
%         is = [is, i];
%     else
%         if length(inds) >= 2
%             dE1s = dE1(is);
%             dE2s = dE2(is);
%             Ds = D(is);
%             
%             % find Z
%             ii = find(indx1 == indx2(i));            
%             z = Zs(ii)/1000; % Altitude in km            
%             Rs = sqrt(Ds.^2+z^2);
%             sins = (2 - 3*(Ds./Rs).^2);
%             
%             
%             
%             
%             % Find best fit dE*R^n = C equation
%             x = log(Rs);
%             y = log(abs(dE1s./sins));
%             p = polyfit(x,y,1);
%             n1 = -p(1);
%             c1 = exp(p(2));
%             
%             y = log(abs(dE2s./sins));
%             p = polyfit(x,y,1);
%             n2 = -p(1);
%             c2 = exp(p(2));
%             
%             n1s = [n1 n1s];
%             n2s = [n2 n2s];
%             
%             dE3s = dE2s-dE1s;
%            
%                                                 
%             % Let's print info
%             fprintf('\n==================== NBP = %i ================\n',indx2(i-1))
%             fprintf('E_fast\tE_total\t\tR\tdE1*f\tdE2*f\tdE3*f\n')
%             Ratio1 = dE1s.*(Rs.^3)./sins/1000;
%             Ratio2 = dE2s.*(Rs.^3)./sins/1000;
%             Ratio3 = dE3s.*(Rs.^3)./sins/1000;
%             
%             for j = 1:length(inds)
%                 fprintf('%6.2f\t%6.2f\t%5.1f\t\t%6.3f\t\t%6.3f\t\t%6.3f\n',...
%                     dE1s(j),dE2s(j),Rs(j),Ratio1(j),Ratio2(j),Ratio3(j))
%             end    
%             
%             fprintf('Best fit exponant for E_fast = %0.2f\n',n1)
%             fprintf('Best fit exponant for E_slow = %0.2f\n',n2)
%             fprintf('f = R*3/(2 - 3*(D/R)^2)\n')
%         end
%         inds = indx2(i);
%         is = i;
%     end
% end
% 
% figure
% subplot(1,3,1)
% hold all
% plot(n1s,'bo')
% plot(xlim,[0 0]+nanmean(n1s))
% box on
% title(sprintf('Fast E exponent = (%0.1f+-%0.f)',nanmean(n1s),nanstd(n1s)))
% xlabel('Sample #')
% ylabel('Exponent')
% 
% 
% subplot(1,3,2)
% hold all
% plot(n2s,'ro')
% plot(xlim,[0 0]+nanmean(n2s))
% box on
% title(sprintf('Slow E exponent = (%0.1f+-%0.f)',nanmean(n2s),nanstd(n2s)))
% xlabel('Sample #')
% ylabel('Exponent')
% 
% subplot(1,3,3)
% hold all
% plot(n1s,'bo')
% plot(n2s,'ro')
% plot(xlim,[0 0]+nanmean([n1s n2s]))
% box on
% title(sprintf('All E exponent = (%0.1f+-%0.f)',nanmean([n1s n2s]),nanstd([n1s n2s])))
% xlabel('Sample #')
% ylabel('Exponent')
% 
% 
% return





%% Generate histogram with closest sensor information

% % Horizontal distance to each sensor
% L = length(pbfat);
% Ds = nan(L,10);
% sen_set = open('sensor_setting.mat');
% 
% 
% % Get indexes of all NBPs (remove none NBPs)
% TypeFileID =  fopen('C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\Types-20141030\20110814-type_info.txt');
% typeData = textscan(TypeFileID,'%f %f');
% 
% % indexes of undertmined 
% Ind_1 = find(isnan(typeData{2})); % NaNs
% Ind_2 = find(typeData{2} == 0); % Type 0 (undetermined)
% Ind_3 = find(typeData{2} == 6); % Typ2 6 (removed from original study);
% 
% validNBPs = setdiff(1:305,[Ind_1; Ind_2 ; Ind_3]);
% 
% for i = 1:10
%     Ds(:,i) = sqrt((pbfax - sen_set.x(i) ).^2 + (pbfay - sen_set.y(i)).^2);
% end
% 
% % Convert to km
% Ds = Ds/1000;
% % figure
% % Dn = reshape(Ds,[1,3050]);
% % [ys, xs] = hist(Dn,.5:1:11);
% % ys(end) = ys(end)/1000;
% % 
% % bar(xs,ys,1);
% % xlabel('Horizontal distance (km)')
% % ylabel('Number of stations')
% % title('Sensors located different distances from NBPs')
% % tools2fig
% 
% 
% 
% % % Let's write a file
% % fID = fopen('C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\distance_to_sensors.txt','w');
% % size(Ds)
% % [mm, ind] = min(Ds(:));
% % [i,j] = ind2sub([305,10],ind);
% % 
% % while ~isnan(ind) || ~(i == 1 && j == 1)
% %     fprintf(fID,'%3.3i\t%2.2i\t%5.1f\n',i,j,Ds(i,j));
% %     fprintf('%3.3i\t%2.2i\t%5.1f\n',i,j,Ds(i,j));
% %        Ds(i,j) = NaN;
% %     [mm, ind] = min(Ds(:));
% %     [i,j] = ind2sub([305,10],ind);%     
% % end
% 
% 
% mDs = [];
% inds =[];
% is = [];
% for i = 1:305;
%     % is this a real NBP, then continue
%     if sum(validNBPs == i)
%         % find minimum distance and the sensor
%         [mD, ind] = min(Ds(i,:));
%         is = [is i];
%         mDs = [mDs mD];
%         inds = [inds ind];
%         % save in an array
%     end
% end
% 
% % Let's write the file;
% mDdata = [is; inds; mDs];
% mDdata = mDdata';
% mDdata = sortrows(mDdata,3);
% 
% fID = fopen('C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\data\closest_to_sensors.txt','w');
% 
% 
% for i = 1:length(mDdata)
%     if mDdata(i,2) > 3
%         mDdata(i,2) = mDdata(i,2) + 1;
%     end
%     
%     fprintf('%3.3i\t%2.2i\t%5.1f\n',mDdata(i,:))
%     fprintf(fID,'%3.3i\t%2.2i\t%5.1f\n',mDdata(i,:));
% 
% end
% 
% fclose(fID);
% return

%% NBP generate normalized peak info
%NBP_normalized_peak_vals

%% NBP types histogram
% Ns = [nnz(type == 1) nnz(type == 2)+nnz(type == 3) nnz(type == 4) nnz(type == 5)]
% L = sum(Ns)
% figure
% hold all
% bar(Ns)
% 
% centers = 1:4;
% K = numel(centers);
% pcts = Ns;
% for k = 1:K
%     text(centers(k),pcts(k),[num2str(round(pcts(k)/L*1000)/10) '% (' num2str(pcts(k)) ')' ],'HorizontalAlignment','center','VerticalAlignment','bottom')
% end
% box on
% ylim([0 160])
% set(gca,'xtick',[1 2 3 4],'xTickLabel',{'A','B','C','D'})
% ylabel('Number of NBPs')
% xlabel('Types')
% tools2fig

%% Check file names vs types
% files = dir('C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\Types-20141030\5\*.png')
% inds = [];
% for i = 1:length(files)
%     inds = [inds str2num(files(i).name(1:4))];
% end
% inds = inds';
% inds2 = find(type == 5);
% 
% c = setdiff(inds2, inds)

%% xy plot and Altitude Histogram

%xy plot

% figure
% subplot(1,2,1)
% hold all
% 
% D = sqrt(pbfax.^2 + pbfay.^2)/1000;
% 
% indx = find(type == 1); plot(pbfax(indx)/1000,pbfay(indx)/1000,'ro','markerfacecolor','r','markersize',3);
% ind_all = indx;
% indx = [find(type == 2); find(type == 3)]; plot(pbfax(indx)/1000,pbfay(indx)/1000,'ko','markerfacecolor','k','markersize',2);
% ind_all = [ind_all; indx];
% indx = find(type == 4); plot(pbfax(indx)/1000,pbfay(indx)/1000,'bo','markerfacecolor','b','markersize',2);
% ind_all = [ind_all; indx];
% indx = find(type == 5); plot(pbfax(indx)/1000,pbfay(indx)/1000,'go','markerfacecolor','g','markersize',2);
% ind_all = [ind_all; indx];
% plot(sns_x,sns_y,'mp','markerfacecolor','m','markersize',5)
% 
% 
% fprintf('Distance Max %0.1f    Min %0.1f    Mean %0.1f  STD %0.1f\n', ...
%     max(D(ind_all)),min(D(ind_all)),mean(D(ind_all)),std(D(ind_all)))
% 
% 
% % Find number of NBPs within groups
% D_all = D(ind_all);
% rangesD = [0 30 60 90 120 150];
% 
% for i = 1:length(rangesD)-1
%     
%     DrInd1 = find(D_all > rangesD(i));
%     DrInd2 = find(D_all <= rangesD(i+1));
%     DrInds = intersect(DrInd1,DrInd2);
%     length(DrInds)
% end
%    
%     
% 
% 
% ylim([-150 200])
% xlim([-130 60])
% florida_map
% text(sns_x,sns_y,sns_n)
% 
% legend('Type A','Type B','Type C','Type D','Sensors')
% 
% box on
% xlabel('East (km)')
% ylabel('North (km)')
% daspect([1 1 1])
% 
% % Indexs of all the types
% ind_type = [find(type==1); find(type==2); find(type==3); find(type==4); find(type==5)];
% 
% 
% pbfaz = pbfaz(ind_type);
% pbfax = pbfax(ind_type);
% pbfay = pbfay(ind_type);
% pbfat = pbfat(ind_type);
% 
% fprintf('Altitude Max %0.1f    Min %0.1f    Mean %0.1f  STD %0.1f\n', ...
%     max(pbfaz),min(pbfaz),mean(pbfaz),std(pbfaz))
% 
% % plot data
% subplot(1,2,2)
% hold all
% [counts,bins] = hist(pbfaz/1000,(500:1000:20000)/1000);
% barh(bins,counts,1.0)
% ylabel('Altitude (km)')
% xlabel('Number of NBPs')
% % if savePlots
% %     saveas(gcf,[bf 'Z-hist.fig'])
% % end
% 
% 
% 
% 
% % spline fit
% cs = csape(bins,counts);
% xx = linspace(0,20,101);
% yy = ppval(cs,xx);
% %;
% %set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
% 
% 
% tools2fig
% 
% fprintf('Altitude range: %0.1f -- %0.1f km\n',min(pbfaz),max(pbfaz))
% 
% % altitude data from smith et al
% data = [156.2	894.3
% 156.2	850.0
% 157.1	829.8
% 162.0	805.7
% 186.1	762.3
% 177.4	740.1
% 186.1	717.9
% 192.8	696.7
% 192.8	675.5
% 226.6	652.4
% 245.8	631.2
% 323.9	610.0
% 378.9	587.8
% 467.6	565.6
% 657.5	544.4
% 761.6	522.2
% 830.0	500.1
% 850.3	477.9
% 933.2	456.7
% 1025.7	435.5
% 1109.6	412.3
% 1233.9	390.2
% 1533.7	369.0
% 1613.8	345.8
% 1488.4	323.7
% 1176.1	303.4
% 970.8	281.2
% 684.5	260.0
% 473.3	237.9
% 343.2	215.7
% 289.2	193.5
% 223.7	172.3
% 207.3	150.1
% 188.0	128.0
% 191.9	105.8
% 170.6	84.6
% 169.7	62.4
% 166.8	40.2
% 162.0	19.0];
% 
% % axis data (for scalling)
% xv1 = [155 0];
% xv2 = [1491  1500];
% yv1 = [916 0];
% yv2 = [40 20];
% 
% % scale
% x = (xv2(2) - xv1(2))/(xv2(1) - xv1(1))*(data(:,1) - xv1(1));
% y = (yv2(2) - yv1(2))/(yv2(1) - yv1(1))*(data(:,2) - yv1(1));
% 
% % Normalize data to match our data
% max(counts)/max(x)
% x = x*max(counts)/max(x);
% 
% ylim([0 20])
% xlim([0 70])
% 
% %barh(bins,counts,1.0)
% %plot(counts,bins,'r-','LineWidth',2)
% plot(x,y,'k-','LineWidth',2)
% %plot(yy,xx,'-g','LineWidth',2);
% [mm ind] = max(yy);
% xx(ind)
% 
% box on
% legend('Current study', '\it{Smith et al. 2004}')
% %pH = plot(counts,bins,'ro');
% 
% 
% % Max pbfas
% for i = 1:5
%     [mm mind] = max(pbfaz);
%     fprintf('%0.1f\t%0.1f\n',sqrt((pbfax(mind)/1000)^2+(pbfay(mind)/1000)^2),pbfaz(mind)/1000)
%     pbfaz(mind) = 12000;
% end
% 
% % Man pbfas
% for i = 1:5
%     [mm mind] = min(pbfaz);
%     fprintf('%0.6f\t%0.1f\t%0.1f\n',pbfat(mind),sqrt((pbfax(mind)/1000)^2+(pbfay(mind)/1000)^2),pbfaz(mind)/1000)
%     pbfaz(mind) = 12000;
% end


%% NBP altitude types histogram

% figure
% hold all
% cnt = [];
% lg = {};
% 
% indx = find(type == 1);
% [counts,bins] = hist(pbfaz(indx)/1000,(500:1000:20000)/1000);
% cnt =  [cnt; counts];
% lg = [lg sprintf('A (%0.1f\\pm%0.1f) km',mean(pbfaz(indx)/1000),std(pbfaz(indx)/1000))];
% %barh(bins,counts,1.0)
% 
% indx = [find(type == 2); find(type == 3)];
% [counts,bins] = hist(pbfaz(indx)/1000,(500:1000:20000)/1000);
% cnt =  [cnt; counts];
% lg = [lg sprintf('B (%0.1f\\pm%0.1f) km',mean(pbfaz(indx)/1000),std(pbfaz(indx)/1000))];
% %barh(bins,counts,1.0)
% 
% indx = find(type == 4);
% [counts,bins] = hist(pbfaz(indx)/1000,(500:1000:20000)/1000);
% cnt =  [cnt; counts];
% lg = [lg sprintf('C (%0.1f\\pm%0.1f) km',mean(pbfaz(indx)/1000),std(pbfaz(indx)/1000))];
% %barh(bins,counts,1.0)
% 
% indx = find(type == 5); 
% [counts,bins] = hist(pbfaz(indx)/1000,(500:1000:20000)/1000);
% cnt =  [cnt; counts];
% lg = [lg sprintf('D (%0.1f\\pm%0.1f) km',mean(pbfaz(indx)/1000),std(pbfaz(indx)/1000))];
% %barh(bins,counts,1.0)
% barh(bins,cnt')
% ylim([8,20])
% box on
% 
% legend(lg)
% tools2fig
%  xlabel('Number of NBPs')
% ylabel('Altitude (km)')



%% PBFA loction error stats
% erFn = 'C:\Users\Sumedhe\Desktop\NBP\PBFA_errors\Combined_method_10km_res_large_test.mat';
% erD = open(erFn);
% 
% % distance
% r = sqrt(erD.x.^2+erD.y.^2);
% 
% ranges = [0 30 60 90 120 150]*1000;
% 
% for i = 2:length(ranges)
%     
%     ind1 = find(r>ranges(i-1));
%     ind2 = find(r<=ranges(i));
%     ind = intersect(ind1,ind2);
%     [I,J] = ind2sub(size(r),ind);
%     
%     dx = nanmean(nanmean(erD.dx3(I,J)));
%     dy = nanmean(nanmean(erD.dy3(I,J)));
%     dz = nanmean(nanmean(erD.dz3(I,J)));
%     dt = nanmean(nanmean(erD.dt3(I,J)));
%     
%     fprintf('<=%2.2i km -- dx =%3.0fm\tdy =%3.0fm\tdz =%4.0fm\tdt =%4.1fus\n',ranges(i),dx,dy,dz,dt*1e6)
%     
% end






%% New time comparison ( Isolated and non isolated)
% tfn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\Times\20110814-time_info-modified.txt';
% 
% fID = fopen(tfn);
% tdata = textscan(fID,'%f %f %f %f %f %f %f','headerlines',2);
% fclose(fID);
% 
% tdata = cell2mat(tdata);
% 
% figure
% tools2fig
% 
% subplot(2,3,1)
% hold all
% [y1, x1] = hist(tdata(:,2),0.5:2.5:21);
% xlabel('Rise time (\mus)')
% 
% %spline fit
% splinefitplot([0 x1],[0 y1],0:.1:21,'r')
% xlim([0 21])
% ylim([0 110])
% box on
% 
% 
% subplot(2,3,2)
% hold all
% [y2, x2] = hist(tdata(:,3),0.5:1.5:15);
% xlim([0 15])
% xlabel('10-90% Rise time (\mus)')
% 
% %spline fit
% splinefitplot([0 x2],[0 y2],0:.1:15,'r')
% xlim([0 15])
% ylim([0 110])
% box on
% 
% subplot(2,3,3)
% hold all
% [y3, x3] = hist(tdata(:,7),0.5:10);
% xlim([0 10])
% xlabel('FWHM (\mus)')
% 
% %spline fit
% splinefitplot([0 x3],[0 y3],0:.1:10,'r')
% xlim([0 10])
% ylim([0 100])
% box on
% 
% 
% subplot(2,3,4)
% hold all
% [y4, x4] = hist(tdata(:,4),0.5:4:35);
% xlim([0 35])
% xlabel('Zero Cross Time (\mus)')
% 
% 
% %spline fit
% splinefitplot([0 x4],[0 y4],0:.1:35,'r')
% xlim([0 35])
% ylim([0 90])
% box on
% 
% 
% subplot(2,3,5)
% hold all
% [y5, x5] = hist(tdata(:,5),0:4:32);
% %xlim([4.5 39])
% xlabel('Negative peak Time (\mus)')
% 
% %spline fit
% splinefitplot(x5,y5,0:1:32,'r')
% xlim([0 35])
% ylim([0 60])
% box on
% 
% 
% 
% subplot(2,3,6)
% hold all
% [y6, x6] = hist(tdata(:,6),20:10:120);
% %xlim([4.5 39])
% xlabel('Total duration (\mus)')
% 
% %spline fit
% splinefitplot([20 x6],[0 y6],20:.5:120,'r')
% xlim([20 120])
% ylim([0 60])
% box on
% 
% % Let's work on isolated ones
% iso = cell2mat(data(:,29));
% iso2 = cell2mat(data(:,22));
% % truly isolated
% ind = intersect(find(iso >= 0.66), find(abs(iso2) >= 0.66));
% 
% 
% subplot(2,3,1)
% [y1, x1] = hist(tdata(ind,2),0.5:2.5:21);
% splinefitplot([0 x1],[0 y1],0:.1:21,'g')
% 
% subplot(2,3,2)
% [y2, x2] = hist(tdata(ind,3),0.5:1.5:15);
% splinefitplot([0 x2],[0 y2],0:.1:15,'g')
% 
% subplot(2,3,3)
% [y3, x3] = hist(tdata(ind,7),0.5:10);
% splinefitplot([0 x3],[0 y3],0:.1:10,'g')
% 
% subplot(2,3,4)
% [y4, x4] = hist(tdata(ind,4),0.5:4:35);
% splinefitplot([0 x4],[0 y4],0:.1:35,'g')
% 
% subplot(2,3,5)
% [y5, x5] = hist(tdata(ind,5),0:4:32);
% splinefitplot(x5,y5,0:1:32,'g')
% 
% subplot(2,3,6)
% [y6, x6] = hist(tdata(ind,6),20:10:120);
% splinefitplot([20 x6],[0 y6],20:.5:120,'g')
% 
% % Let's work on NONE-isolated ones
% iso = cell2mat(data(:,29));
% iso2 = cell2mat(data(:,22));
% % truly isolated
% ind = intersect(find(iso >= 0.66), find(abs(iso2) >= 0.66));
% 
% ind = setxor(1:305,ind);
% 
% subplot(2,3,1)
% [y1, x1] = hist(tdata(ind,2),0.5:2.5:21);
% splinefitplot([0 x1],[0 y1],0:.1:21,'b')
% legend('All','Isolated','Non-isolated')
% 
% 
% subplot(2,3,2)
% [y2, x2] = hist(tdata(ind,3),0.5:1.5:15);
% splinefitplot([0 x2],[0 y2],0:.1:15,'b')
% legend('All','Isolated','Non-isolated')
% 
% 
% subplot(2,3,3)
% [y3, x3] = hist(tdata(ind,7),0.5:10);
% splinefitplot([0 x3],[0 y3],0:.1:10,'b')
% legend('All','Isolated','Non-isolated')
% 
% 
% subplot(2,3,4)
% [y4, x4] = hist(tdata(ind,4),0.5:4:35);
% splinefitplot([0 x4],[0 y4],0:.1:35,'b')
% legend('All','Isolated','Non-isolated')
% 
% 
% subplot(2,3,5)
% [y5, x5] = hist(tdata(ind,5),0:4:32);
% splinefitplot(x5,y5,0:1:32,'b')
% legend('All','Isolated','Non-isolated')
% 
% 
% subplot(2,3,6)
% [y6, x6] = hist(tdata(ind,6),20:10:120);
% splinefitplot([20 x6],[0 y6],20:.5:120,'b')
% legend('All','Isolated','Non-isolated')


%% Timing (new) of different types

% tfn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\Times\20110814-time_info-modified.txt';
% 
% fID = fopen(tfn);
% tdata = textscan(fID,'%f %f %f %f %f %f %f','headerlines',2);
% fclose(fID);
% 
% tdata = cell2mat(tdata);
% 
% 
% 
% 
% 
% figure
% tools2fig
% str1 = {};
% str2 = {};
% str3 = {};
% str4 = {};
% 
% typesChar = ['0','A','B','C','D','E','F','G','H'];
% 
% 
% types = type;
% rst = [];
% fwhmt = [];
% zct = [];
% 
% for i = 1:4
%     
%     
%     if i == 1
%         inds = find(types == 1);
%         col = 'r';
%     elseif i ==2
%         inds = [find(types == 2); find(types == 3)];
%         col = 'k';
%     elseif i == 3
%         inds = find(types == 4);
%         col = 'b';
%     elseif i== 4
%         inds = find(types == 5);
%         col = 'g';
%     end
%     
% %         if i == 1
% %             inds = [find(types == 1); find(types == 5); find(types == 7)];
% %             col = 'r';
% %         elseif i ==2
% %             inds = [find(types == 6); find(types == 8)];
% %             col = 'b';
% %         elseif i == 3
% %             inds = [find(types == 2); find(types == 4); find(types == 3)];
% %             col = 'g';
% %         elseif i== 4
% %             inds = 1:length(types);
% %             col = 'k';
% %         end
% %     
%     
%     
%     
% %     switch i
% %         case 0
% %             inds = find(types == 0);
% %             col = 'k';
% %         case 1
% %             inds = find(types == 1);
% %             col = 'r';
% %         case 2
% %             inds = find(types == 2);
% %             col = 'g';
% %         case 3
% %             inds = find(types == 3);
% %             col = 'b';
% %         case 4
% %             inds = find(types == 4);
% %             col = 'c';
% %         case 5
% %             inds = find(types == 5);
% %             col = 'm';
% %     end
% %    
%     %     % Rise time
%     %     subplot(2,3,1)
%     %     hold all
%     %     [y1, x1] = hist(tdata(inds,2),0.5:2.5:21);
%     %     xlabel('Rise time (\mus)')
%     %
%     %     %spline fit
%     %     splinefitplot([0 x1],[0 y1],0:.1:21,col)
%     %     xlim([0 21])
%     %     ylim([0 50])
%     %     box on
%     
%   
%     
%     % 10-90 rise time
%     subplot(2,2,1)
%     hold all
%     [y2, x2] = hist(tdata(inds,3),0.5:1:15);
%     rst = [rst; y2];
%     rstx = x2;
%     xlabel('10-90% Rise time (\mus)')
%     
%     % print info for largest 5 rise times
%     if i== 4
%         rsdata = tdata(inds,3);
%         for kk = 1:5            
%             [mmm, mmi] = max(rsdata);
%             fprintf('%i\t%0.1f\n',inds(mmi),mmm)
%             rsdata(mmi) = -inf;            
%         end
%     end
%     
%     
%     %spline fit
%     splinefitplot([0 x2],[0 y2],0:.1:15,col)
%     xlim([0 15])
%     %ylim([0 50])
%     box on
%     
%     mn = nanmean(tdata(inds,3));
%     stdd = nanstd(tdata(inds,3));
%     str1 = [str1 sprintf('%s (%0.1f\\pm%0.1f) \\mus',typesChar(i+1),mn,stdd)];
%     legend(str1);
%     
%     % FWHM
%     subplot(2,2,2)
%     hold all
%     [y3, x3] = hist(tdata(inds,7),0.5:10);
%     xlabel('FWHM (\mus)')
%     fwhmt = [fwhmt; y3];
%     fwhmx = x3;
%     
%     
%     %spline fit
%     splinefitplot([0 x3],[0 y3],0:.1:10,col)
%     xlim([0 10])
%     %ylim([0 45])
%     box on
%     
%     mn = nanmean(tdata(inds,7));
%     stdd = nanstd(tdata(inds,7));
%     str2 = [str2 sprintf('%s (%0.1f\\pm%0.1f) \\mus',typesChar(i+1),mn,stdd)];
%     legend(str2);
%     
%     % 0 cross time
%     subplot(2,2,3)
%     hold all
%     [y4, x4] = hist(tdata(inds,4),2.5:5:35);
%     xlabel('Zero Cross Time (\mus)')
%     zct = [zct; y4];
%     zctx = x4;
%     
%     %spline fit
%     splinefitplot([0 x4],[0 y4],0:.1:35,col)
%     xlim([0 35])
%     %ylim([0 40])
%     box on
%     
%     
%     mn = nanmean(tdata(inds,4));
%     stdd = nanstd(tdata(inds,4));
%     str3 = [str3 sprintf('%s (%0.1f\\pm%0.1f) \\mus',typesChar(i+1),mn,stdd)];
%     legend(str3);
%     
%     
%     % Negative peak time
%     subplot(2,2,4)
%     hold all
%     [y5, x5] = hist(tdata(inds,5),0:4:32);
%     xlabel('Negative peak Time (\mus)')
%     
%     %spline fit
%     splinefitplot(x5,y5,0:.1:32,col)
%     xlim([0 35])
%     %ylim([0 32])
%     box on
%     
%     mn = nanmean(tdata(inds,5));
%     stdd = nanstd(tdata(inds,5));
%     str4 = [str4 sprintf('%s (%0.1f\\pm%0.1f) \\mus',typesChar(i+1),mn,stdd)];
%     legend(str4);
%     
%     %     % Total duration
%     %     subplot(2,3,6)
%     %     hold all
%     %     [y6, x6] = hist(tdata(inds,6),20:10:120);
%     %     %xlim([4.5 39])
%     %     xlabel('Total duration (\mus)')
%     %
%     %     %spline fit
%     %     splinefitplot([20 x6],[0 y6],20:.5:120,col)
%     %     xlim([20 120])
%     %     ylim([0 25])
%     %     box on
%     
%     
% end
% 
% tools2fig
% 
% figure(101)
% 
% % Rise time bar plot
% subplot(2,3,4)
% hold all
% h = bar(rstx',rst','BarWidth',1);
% set(h(1),'facecolor','r')
% set(h(2),'facecolor','k')
% set(h(3),'facecolor','b')
% set(h(4),'facecolor','g')
% %plot(rstx',rst','linewidth',2)
% ylabel('Number of NBPs')
% xlabel('10-90% Time (\mus)')
% legend(str1)
% box on
% set(gca,'xtick',[0:1.0:15],'xticklabel',{0 ,'','','','', 5, '','','','', 10,'','','','', 15})
% 
% 
% 
% 
% 
% % FWHM time bar plot
% subplot(2,3,5)
% hold all
% h = bar(fwhmx',fwhmt',1);
% set(h(1),'facecolor','r')
% set(h(2),'facecolor','k')
% set(h(3),'facecolor','b')
% set(h(4),'facecolor','g')
% legend('Type A','Type B','Type C','Type D')
% xlabel('FWHM Time (\mus)')
% legend(str2)
% box on
% 
% % zero cross time time bar plot
% subplot(2,3,6)
% hold all
% h = bar(zctx',zct',1);
% set(h(1),'facecolor','r')
% set(h(2),'facecolor','k')
% set(h(3),'facecolor','b')
% set(h(4),'facecolor','g')
% legend('Type A','Type B','Type C','Type D')
% xlabel('Zero Cross Time (\mus)')
% legend(str3)
% box on
% set(gca,'xtick',[0:5:40],'xticklabel',{'', 5, '', 15, '', 25, '', 35, '',45})
% 
% 
% %subplot(2,2,1)
% %legend('A','B','C','D','E','F','G','H')
% str = {};
% % times in one plot
% % 10-90 rise time
% figure
% hold all
% inds = 1:length(types);
% [y2, x2] = hist(tdata(inds,3),0.5:1.5:15);
% xlabel('Time (\mus)')
% 
% 
% 
% 
% %spline fit
% splinefitplot([0 x2],[0 y2],0:.1:15,'r')
% xlim([0 15])
% %ylim([0 50])
% box on
% 
% mn = nanmean(tdata(inds,3));
% stdd = nanstd(tdata(inds,3));
% str = [str sprintf('%s (%0.1f\\pm%0.1f) \\mus','10-90 % rise time',mn,stdd)];
% legend(str);
% 
% % FWHM
% [y3, x3] = hist(tdata(inds,7),0.5:10);
% xlabel('FWHM (\mus)')
% 
% %spline fit
% splinefitplot([0 x3],[0 y3],0:.1:10,'g')
% xlim([0 10])
% %ylim([0 45])
% box on
% 
% mn = nanmean(tdata(inds,7));
% stdd = nanstd(tdata(inds,7));
% str = [str sprintf('%s (%0.1f\\pm%0.1f) \\mus','FWHM',mn,stdd)];
% legend(str);
% 
% % 0 cross time
% 
% [y4, x4] = hist(tdata(inds,4),0.5:4:35);
% xlabel('Zero Cross Time (\mus)')
% 
% 
% %spline fit
% splinefitplot([0 x4],[0 y4],0:.1:35,'b')
% xlim([0 35])
% %ylim([0 40])
% box on
% 
% 
% mn = nanmean(tdata(inds,4));
% stdd = nanstd(tdata(inds,4));
% str = [str sprintf('%s (%0.1f\\pm%0.1f) \\mus','Zero cross time',mn,stdd)];
% legend(str);
% 
% 
% 
% 
% % Negative peak time
% 
% % [y5, x5] = hist(tdata(inds,5),0:4:32);
% % xlabel('Negative peak Time (\mus)')
% % 
% % %spline fit
% % splinefitplot(x5,y5,0:.1:32,'k')
% % xlim([0 35])
% % %ylim([0 32])
% % box on
% 
% % mn = nanmean(tdata(inds,5));
% % stdd = nanstd(tdata(inds,5));
% % str = [str sprintf('%s (%0.1f\\pm%0.1f) \\mus','Neg peak time',mn,stdd)];
% % legend(str);
% 
% figure(101)
% subplot(2,3,1)
% hold all
% [y2, x2] = hist(tdata(inds,3),1.25:2.5:35);
% [y3, x3] = hist(tdata(inds,7),1.25:2.5:35);
% [y4, x4] = hist(tdata(inds,4),1.25:2.5:35);
% bar(x4',[y2;y3;y4]')
% set(gca,'xtick',0:2.5:40,'xticklabel',{0,'', 5, '', 10, '', 15, '', 20, '', 25, '', 30, '', 35, '', 40})
% 
% % [y2, x2] = hist(tdata(inds,3),0.5:1:15);
% % [y3, x3] = hist(tdata(inds,7),0.5:10);
% % [y4, x4] = hist(tdata(inds,4),0.5:5:35);
% % 
% % plot(x2,y2)
% % plot(x3,y3)
% % plot(x4,y4)
% 
% legend('Rise time','Fall time','Zero cross time')
% legend(str)
% box on
% ylabel('Number of NBPs')
% xlabel('Time (\mus)')
% 
% tools2fig

%% Normalized peak distribution

% %file name for NBP data
% fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814-NBP_info2.xlsx';
% sheet = 1;
% xlRange = 'A3:AD307';
% ndata = xlsread(fn, sheet, xlRange);
% 
% % directory name for power data
% %dn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814_plots\MeanPowers\';
% dn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\Noramlized_peaks\';
% 
% indx = ndata(:,1);
% 
% % To put in the title
% pulse_kinds = { ...
%     '0: Not claer'
%     'A: Clean bipolar'
%     'B: Pulses after neg overshoot'
%     'C: Pulses befora and after neg overshoot'
%     'D: Pulses before zero cross'
%     'E: Pulses in rising edge'};
% 
% types = type;
% 
% figure(100)
% subplot(2,3,3)
% hold all
% np = [];
% 
% clor = {'r','k','b','g'};
% xs = 0:0.1:45;
% TypeStr ={'A','B','C','D'};
% lg_str = {};
% 
% cts = [];
% 
% for type = 1:4
%     %inds = [find(types == 6); find(types == 8)] ;
%     switch type
%         case 1
%             inds = find(types == 1);
%         case 2
%             inds = [find(types == 2); find(types == 3)];
%         case 3
%             inds = find(types == 4);
%         case 4
%             inds = find(types == 5);
%     end
% 
%     np1 = [];
% 
%     
%     
% 
%     for i = 1:length(inds)
% 
%         % Power file name
%         %pfn = sprintf('%sMeanPowersMeanPower-%3.3i.mat',dn,indx(i));
%         pfn = sprintf('%s%3.3i-normalized.fig',dn,indx(inds(i)));
% 
%         if exist(pfn,'file')
%             fprintf('Working on %s\n',pfn)
%             % b = open(pfn);
%             fgh = openfig(pfn,'reuse','invisible');
%             d = guidata(fgh);
%             close(fgh);
%             np1 = [np1 nanmean(d.normPeaks)];
%             np = [np nanmean(d.normPeaks)];
% 
%         else
%             fprintf('File not found %s\n',pfn)
%         end
%     end
% 
%     [counts,bins] = hist(np1,2.5:5:42.5)
%     cts = [cts; counts];
%  
%     %splinefitplot([0 bins],[0 counts],xs,clor{type})
%     lg_str = [lg_str sprintf('%s (%0.1f\\pm%0.1f) V/m' ,TypeStr{type},mean(np1),std(np1))]
%     
% end 
% 
% h = bar(bins',cts',1);
% set(h(1),'facecolor','r')
% set(h(2),'facecolor','k')
% set(h(3),'facecolor','b')
% set(h(4),'facecolor','g')
% box on
% 
% 
% % [counts,bins] = hist(np,2.5:5:42.5);
% % splinefitplot([0 bins],[0 counts],xs,'m')
% lg_str = [lg_str sprintf('%s (%0.1f\\pm%0.1f) V/m' ,'All NBPs',mean(np),std(np))]
% 
% legend(lg_str)
% % mean(np)
% % std(np)
% 
% 
% 
% figure
% hold all
% counts = [0 27    84    49    20    16     4     1     1     0];
% bins = [0 2.5000    7.5000   12.5000   17.5000   22.5000   27.5000   32.5000   37.5000 42.5];
% %mean = 11.0152;
% %sdev = 6.2302;
% %plot(bins,counts,'ro')
% xs = 0:0.1:45;
% splinefitplot(bins,counts,xs,'r')
% ylim([0 90])
% xlabel('Normalized peak amplitude (V/m)')
% ylabel('Number of NBPs')
% 
% tools2fig



%% Power distributions according to the types

% % file name for NBP data
% fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814-NBP_info2.xlsx';
% sheet = 1;
% xlRange = 'A3:AD307';
% ndata = xlsread(fn, sheet, xlRange);
% 
% % directory name for power data
% %dn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814_plots\MeanPowers\';
% dn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\PowerCurves5MHz-3\';
% 
% indx = ndata(:,1);
% 
% % To put in the title
% pulse_kinds = { ...
%     '0: Not claer'
%     'A: Clean bipolar'
%     'B: Pulses after neg overshoot'
%     'C: Pulses befora and after neg overshoot'
%     'D: Pulses before zero cross'
%     'E: Pulses in rising edge'};
% 
% types = type;
% 
% figure
% clor = {'r','k','b','g'};
% 
% saveD.ind = [];
% saveD.totalMeanP = [];
% 
% peakP = [];
% for type = 1:4
%     %inds = [find(types == 6); find(types == 8)] ;
%     switch type
%         case 1
%             inds = find(types == 1);
%         case 2
%             inds = [find(types == 2); find(types == 3)];
%         case 3
%             inds = find(types == 4);
%         case 4
%             inds = find(types == 5);
%     end
% 
%     
% 
%     fs = [];
%     amps = [];
%     coudnt_get = 0;
% 
%     for i = 1:length(inds)
% 
%         % Power file name
%         %pfn = sprintf('%sMeanPowersMeanPower-%3.3i.mat',dn,indx(i));
%         pfn = sprintf('%s%3.3i-power_dist.fig',dn,indx(inds(i)));
% 
%         if exist(pfn,'file')
%             fprintf('Working on %s\n',pfn)
%             % b = open(pfn);
%             fgh = openfig(pfn,'reuse','invisible');
%             d = guidata(fgh);
%             close(fgh);
% 
%             % save data for future use
%             %f = b.f;
%             %amp = b.amp;
%             
%             % Save data for future us (one time business)
%             %saveD.ind = [saveD.ind dn,indx(inds(i)) ];
%             %saveD.totalMeanP = [saveD.totalMeanP d.totalMeanP];
%             
%             f = d.fx;
%             amp = d.dataMeanP;
%           
% %             [mm indf] = max(amp);
% %             df = (f(indf+1)+f(indf))/2 - (f(indf)+f(indf-1))/2;
% %             peakP = [peakP mm*df];
%             
%             
%             % Normalize
%             amp = amp/max(amp);
% 
%             try
%                 fs = [fs; log10(f')];
%                 amps = [amps; amp'];
%             catch
%                 coudnt_get = coudnt_get + 1;
%             end
%         else
%             fprintf('File not found %s\n',pfn)
%         end
%     end
% 
%     subplot(2,3,type)
%     hold all
%     mn = mean(amps);
%     er = std(amps);
%     
%    
%     fh = fill([fs(1,2:end) fliplr(fs(1,2:end))],[mn(2:end)+er(2:end) fliplr(mn(2:end))-fliplr(er(2:end))],[128/250 128/250 190/250],'edgecolor','none');
%     set(get(get(fh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
%     %alpha(ph,0.2)
%     plot(fs(1,:),mean(amps), clor{type},'LineWidth',2)
%     xlabel('log_{10}(f in Hz)')
%     ylabel('<Normalized power density>')
%     legend(sprintf('Avg curve for %i NBPs',length(inds)))
%     title(pulse_kinds{type+1})
%     box on
%     
%     % FIx y limit
%     ylim([0 1.2])
%     
%     subplot(2,3,6)
%     hold all
%     plot(fs(1,:),mean(amps),clor{type}, 'LineWidth',2)
%     
%     
% end 
% 
% % save([dn 'totalMeanP.mat'],'-Struct','saveD')
% 
% subplot(2,3,6)
% ylim([0 1.2])
% box on
% xlabel('log_{10}(f in Hz)')
% ylabel('<Normalized power density>')
% legend('A','B','C','D')
% 
% % figure
% % plot(peakP)
% % fprintf('Max = %.3e   Min %.3e Mean %0.3e STD %0.3e',...
% %     max(peakP),min(peakP),mean(peakP),std(peakP))
% 
% tools2fig

%% Total mean power distribution


% fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\PowerCurves5Mhz-3\totalMeanP.mat';
% 
% a = open(fn);
% avgPs = a.totalMeanP;
% 
% 
% 
% figure
% hold all
% x = log10(avgPs);
% 
% % Histogram
% [ys, xs] = hist(x,5.875:0.25:9.5);
% h = bar(xs,ys,1);
% %set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
% 
% 
% min(avgPs)
% max(avgPs)
% 
% 
% 
% box on
% xlabel('log_{10}(<Power in watts>)')
% ylabel('Number of NBPs')
% legend(['\mu = ' sprintf('%0.2e W',mean(avgPs)) '    \sigma = \pm' sprintf('%0.2e W',std(avgPs))],'Location','northwest')
% 
% figure
% hold all
% splinefitplot(xs,ys,6.5:0.01:8.4,'r')


%% Peak currents histogram
% figure
% hist(abs(ndata(:,15)),0:10:120)
% title('Peak Current Distributions')
% 
% tools2fig
% I = abs(ndata(:,15));
% types = type;
% 
% cts = [];
% allInds = [];
% lg = {};
% typestr = {'A','B','C','D'};
% for type = 1:4
%     %inds = [find(types == 6); find(types == 8)] ;
%     switch type
%         case 1
%             inds = find(types == 1);
%         case 2
%             inds = [find(types == 2); find(types == 3)];
%         case 3
%             inds = find(types == 4);
%         case 4
%             inds = find(types == 5);
%     end
%     allInds = [allInds; inds];
%     [counts, bins] = hist(I(inds),5:10:125);
%     cts = [cts; counts];
%     
%     lg = [lg sprintf('%s: (%0.1f\\pm%0.1f) kA',typestr{type},nanmean(I(inds)),std(I(inds)))];
% end
% 
% 
% figure
% subplot(1,2,1)
% h = bar(bins',cts');
% set(h(1),'facecolor','r')
% set(h(2),'facecolor','k')
% set(h(3),'facecolor','b')
% set(h(4),'facecolor','g')
% legend(lg)
% xlabel('Peak currents (kA)')
% ylabel('Number of NBPs')
% 
% subplot(1,2,2)
% hist(I(allInds),5:10:125);
% xlabel('Peak currents (kA)')
% legend(sprintf('All (%i) NBPs: (%0.1f\\pm%0.1f) kA',...
%     length(allInds),mean(I(allInds)),std(I(allInds))))
% 
% fprintf('***** Peak Currents *******')
% fprintf('Max = %.1f kA\n',nanmax(abs(I(allInds))))
% fprintf('Min = %.1f kA\n',nanmin(abs(I(allInds))))
% fprintf('Avg = %.1f kA\n',nanmean(abs(I(allInds))))
% fprintf('STD = %.1f kA\n',nanstd(abs(I(allInds))))
% fprintf('MOD = %.1f kA\n',mode(abs(I(allInds))))
% 
% 
% tools2fig


%% Isolation
% figure
% subplot(1,2,1)
% hold all
% %rectangle('Position',
% ph1 = patch([.66 11 11 .66],[-0.5 -0.5 0.66 0.66],[209,234,211]/255);
% ph2 = patch([-0.5 -0.5 0.66 0.66],[.66 11 11 .66],[203,215,232]/255);
% ph3 = patch([0.66 11 11 0.66],[.66 0.66 11 11],[218,176,176]/255);
% ph4 = patch([-.5 0.66 0.66 -.5],[-.5 -.5  .66 0.66],[100,176,176]/255);
% %alpha([ph1,ph2,ph3],0.4)
% set(get(get(ph1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(ph2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(ph3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(ph4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% 
% types = type;
% typestr = {'A','B','C','D'};
% 
% iso1 = cell2mat(data(:,29));
% iso2 = cell2mat(data(:,22));
% 
% allInds = [];
% 
% dtdt = [];
% 
% for type = 1:4
%     %inds = [find(types == 6); find(types == 8)] ;
%     switch type
%         case 1
%             inds = find(types == 1);
%             mc = 'r';
%             mz = 5;
%         case 2
%             inds = [find(types == 2); find(types == 3)];
%             mc = 'k';
%             mz = 4;
%         case 3
%             inds = find(types == 4);
%             mc = 'b';
%             mz = 4;
%         case 4
%             inds = find(types == 5);
%             mc = 'g';
%             mz = 4;
%     end
% 
%     plot(ndata(inds,29),abs(ndata(inds,22)),'o','markeredgecolor',mc,'markerfacecolor',mc,'markersize',mz)
%     
%     % convet nan to inf because if there is no time available most probably
%     % it should mean that no close by pulses.
%     iso1(isnan(iso1)) = inf;
%     iso2(isnan(iso2)) = inf;
% 
%     
%     fprintf('\nType %s: Within Type  All NBPs\n',typestr{type})
%     ind = intersect(find(iso1(inds) >= 0.66), find(abs(iso2(inds)) >= 0.66));
%     fprintf('Truly isolated = %0.1f%%\t %0.1f%%\n',length(ind)/length(inds)*100,length(ind)/226*100);
%     %dtdt(type,1) = length(ind)/length(inds)*100;
%     dtdt(type,1) = length(ind);
%     
%     ind = intersect(find(iso1(inds) > 0.66), find(abs(iso2(inds)) < 0.66));
%     fprintf('After IC       = %0.1f%%\t %0.1f%%\n',length(ind)/length(inds)*100,length(ind)/226*100);
%     %dtdt(type,2) = length(ind)/length(inds)*100;
%     dtdt(type,2) = length(ind);
%     
%     ind = intersect(find(iso1(inds) < 0.66), find(abs(iso2(inds)) > 0.66));
%     fprintf('Before IC      = %0.1f%%\t %0.1f%%\n',length(ind)/length(inds)*100,length(ind)/226*100);
%     %dtdt(type,3) = length(ind)/length(inds)*100;
%     dtdt(type,3) = length(ind);
%     
%     ind = intersect(find(iso1(inds) < 0.66), find(abs(iso2(inds)) < 0.66));
%     fprintf('Within IC       = %0.1f%%\t %0.1f%%\n',length(ind)/length(inds)*100,length(ind)/226*100);
%     %dtdt(type,4) = length(ind)/length(inds)*100;
%     dtdt(type,4) = length(ind);
%     
%     allInds = [allInds; inds];
%     
%     
%     
% end
% 
% 
% legend('Type A','Type B','Type C','Type D','location','east')
% 
% 
% 
% 
% box on
% axis on
% 
% tools2fig
% set(gca, 'LineWidth', 1)
% xlabel('Post pulse time (s)')
% ylabel('Pre pulse time (s)')
% axis equal
% xlim([-.5 11])
% ylim([-.5 11])
% 
% % Statistics All NBPs
% fprintf('\nAll NBPs: Within Type  All NBPs\n',typestr{type})
% inds = allInds;
% ind = intersect(find(iso1(inds) >= 0.66), find(abs(iso2(inds)) >= 0.66));
% fprintf('Truly isolated = %0.1f%%\t %0.1f%%\n',length(ind)/length(iso1(inds))*100,length(ind)/226*100);
% %dtdt(5,1) = length(ind)/length(inds)*100;
% dtdt(5,1) = length(ind);
% 
% ind = intersect(find(iso1(inds) > 0.66), find(abs(iso2(inds)) < 0.66));
% fprintf('After IC       = %0.1f%%\t %0.1f%%\n',length(ind)/length(iso1(inds))*100,length(ind)/226*100);
% %dtdt(5,2) = length(ind)/length(inds)*100;
% dtdt(5,2) = length(ind);
% 
% ind = intersect(find(iso1(inds) < 0.66), find(abs(iso2(inds)) > 0.66));
% fprintf('Before IC      = %0.1f%%\t %0.1f%%\n',length(ind)/length(iso1(inds))*100,length(ind)/226*100);
% %dtdt(5,3) = length(ind)/length(inds)*100;
% dtdt(5,3) = length(ind);
% 
% ind = intersect(find(iso1(inds) < 0.66), find(abs(iso2(inds)) < 0.66));
% fprintf('Within IC       = %0.1f%%\t %0.1f%%\n',length(ind)/length(iso1(inds))*100,length(ind)/226*100);
% %dtdt(5,4) = length(ind)/length(inds)*100;
% dtdt(5,4) = length(ind);
% 
% 
% subplot(1,2,2)
% hold all
% h = bar(dtdt);
% set(gca,'xtick',[1 2 3 4 5],'xticklabel',{'Type A','Type B','Type C','Type D','All'})
% box on
% ylabel('Number of NBPs')
% legend('Isolated','After IC','Before IC','Within IC','location','northwest')
% 
% ybuff=0.5;
% % Put %values
% totals = [2 151 30 43 226];
% for i=1:length(h)
%     XDATA=get(get(h(i),'Children'),'XData');
%     YDATA=get(get(h(i),'Children'),'YData');
%     for j=1:size(XDATA,2)        
%         x=XDATA(1,j)+(XDATA(3,j)-XDATA(1,j))/2;
%         y=YDATA(2,j)+ybuff;
%         t=[num2str(YDATA(2,j)/totals(j)*100,3) ,'%'];
%         text(x,y,t,'Color','k','HorizontalAlignment','left','Rotation',90)
%     end
% end
% 
% ylim([0 100])
% 
% % a b labels
% for ii = 1:2
%     subplot(1,2,ii)
%     yl = ylim;
%     xl = xlim;
%     text(xl(1)- range(xl)*0.20,yl(2),['(',char(ii+96),')'],'fontsize',15)
% end

%% DONT CHANGE ANYTHING BELOW %%

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

function splinefitplot(x,y,xs,col,norm)

%spline fit
cs = csape(x,y);
yy = ppval(cs,xs);

if nargin < 5 
    norm = 0;
end

if norm
    mm = max(yy);
    yy = yy/mm;
    y = y/mm;
end

plot(xs,yy,[col '-'],'linewidth',2)
H = plot(x(2:end),y(2:end),[col 'o']);
set(get(get(H,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
box on


function NBP_plot_all_verticle(plot_range,aa)

% This function is written to catogorize NBPs according to their types.


% Open the NBP file
data = xlsread('C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814-NBP_info2.xlsx');

% Starting index acoording to column 1 
% plotting range

% base foloder for saving plots
bf = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\Stacked_plots\';


% Start the program


% Load data
% get plotter2 data
h=guidata(findall(0,'Tag','plotter2'));
sen_set = h.sen_set;
g = h.g;

% Turn off all plots
g.lpgraphs = zeros(1,60);
g.chgraphs = zeros(1,60);


% Seconds are most probably 0
g.ss = 0;


for index = plot_range
    % to ingnore header we need to add 2 lines to raws
    t = data(index+2,10);
    
    % If t is NaN, there is no PBFA for this NBP and let's go to the next
    % one
    if ~isnan(t)
        
        % PBFA location
        x0 = data(index+2,11);
        y0 = data(index+2,12);
        z0 = data(index+2,13);
        
        % Time input for plotter
        g.hh = floor(t/3600);
        g.mm = floor(((t - 3600*g.hh)/60)/5)*5;
        g.t1 = t - 0.0006;
        g.t2 = t + 0.0006;
        
        % Sensors to use
        sns = [1 2 3 5 6 8 9 10];
        
        % Channels to use
        ch = 3;
        
        % Let's turn on the just the channels we need
        g.chgraphs(sns*3 - (3-ch)) = 1;
        
        % Generate figure
        gen_multi_figure(index,g,sen_set,t,x0,y0,z0,aa)
        
        
        fg = gcf;
              
        
        % save the figure
%         set(fg, 'PaperUnits', 'centimeters');
%         set(gcf, 'PaperPosition', [0 0 40 25])
%          saveas(fg,sprintf('%s/%3.3i-stacked_plots.fig',...
%              bf,index))
%          saveas(fg,sprintf('%s/%3.3i-stacked_plots.png',...
%              bf,index))
%          
%          delete(fg)
        
    end
    
    
end


disp('done')


function gen_multi_figure(ind,g,sen_set,t0,x0,y0,z0,aa)

% Generate file names
a = file_names_generator(g,sen_set);

% Check unavailable files if you want
% fprintf('%s\n',a.absant_fn{:})

% TIme shifts
tshift=sen_set.t_shift;

% Channel to load
ch = 3;

% Load ch data
fgh = figure('units','normalized','outerposition',[0.05 0.05 0.95 0.95],'visible','off');
hold all

R = sqrt((sen_set.x-x0).^2 + (sen_set.y - y0).^2 + (sen_set.z - z0).^2);
D = sqrt((sen_set.x-x0).^2 + (sen_set.y - y0).^2)/1000;

wbh = waitbar(0,'Please wait..','name','plotting');

vOffset = 0;

%[Dn, sortIn] = sort(D,'descend');
[Dn, sortIn] = sort(D);

wbcounter = 0;

colors = [
         0         0    1.0000
    1.0000         0         0
         0    1.0000         0
         0         0    0.1724
    1.0000    0.1034    0.7241
    1.0000    0.8276         0
         0    0.3448         0
    0.5172    0.5172    1.0000
    0.6207    0.3103    0.2759
         0    1.0000    0.7586
         0    0.5172    0.5862];
counter = 0;

voffsets = 0;

for ind1 = sortIn
    
    i = ind1*3 - (3-ch);
    
    wbcounter = wbcounter + 1;
    waitbar(wbcounter/60,wbh)
    
    if ~strcmp(a.ch_fn{i},'')
        % Arrival time
        sn = ceil(i/3);
        
        % Lets load 75us from both sides of arrival time
        %g.t1 = t0 + R(sn)/3.0e8 - 25e-6;
        %g.t2 = t0 + R(sn)/3.0e8 + 75e-6;
        g.t1 = t0 + R(sn)/3.0e8 - aa.t1t2(1)*1e-6;
        g.t2 = t0 + R(sn)/3.0e8 + aa.t1t2(2)*1e-6;
        
        
        [t, y, ch_freq_str] = FA_Extract1(a.ch_fn{i},a.h_fn{i},g.t1,g.t2,tshift(i),sen_set,i);
        
        t = (t - g.t1)*1e6;
        if ~isempty(t) && range(y)>0.001
            
            
%             filD = fdesign.highpass('Fst,Fp,Ast,Ap',0.1,0.2,40,1,4000);
%             Hd = design(filD,'butter');
%             %Then you just filter with
%             
%             y = filter(Hd,y);

            %y(isnan(y)) = 0;
            %y = detrend(y);           
            
            
            
            counter = counter + 1;
            
            y = y*g.factor(i);
            offset = nanmean(y(1:50));
            yn = y-offset+vOffset;
            
            t = t + aa.dts(counter);
            
            %Np = 100;
            %NN = floor(length(y)/Np);
            %t = mean(reshape(t(1:NN*Np),[Np,NN]));
            %yn =mean(reshape(yn(1:NN*Np),[Np,NN]));
            

            plot(t,yn,'color',colors(ind1,:),'LineWidth',0.5)
            
            xlim([0, g.t2-g.t1]*1e6)
            str = [a.ch_legend{i}(1:3) ' (' sprintf('%0.1f',D(sn)) ' km)'];
            ylabel('E-change (V/m)')
            xlabel('Time (\mus)')
           text(aa.textx-60,nanmean(yn(end-50:end)),str,'HorizontalAlignment','left','VerticalAlignment','bottom')
            
            if length(aa.spacing)== 1
                vOffset = vOffset + aa.spacing;
            else
                vOffset = vOffset + aa.spacing(counter);
            end
            
            voffsets = [voffsets vOffset];
            
                
        end
    end
    
end

set(gca,'ytick',voffsets)
box on

title(sprintf('Index = %3.3i     %0.6f s   (%0.1f, %0.1f, %0.1f) km',ind,t0,x0/1000,y0/1000,z0/1000))


delete(wbh)
set(fgh,'visible','on')


function [fittedX,fittedY,P,S] = best_fit(x,y)

coeffs = polyfit(x, y, 1);
% Get fitted values
fittedX = linspace(min(x), max(x), 200);
fittedY = polyval(coeffs, fittedX);

% Correlation
P = corr(x,y,'type','Pearson'); 
S = corr(x,y,'type','spearman');


function [radarFiles , radarFullFiles] = getRadarFiles(g,sen_set)

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



