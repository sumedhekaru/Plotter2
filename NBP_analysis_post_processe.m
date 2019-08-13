function NBP_analysis_post_processe
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

pbfat = ndata(:,10);
pbfax = ndata(:,11);
pbfay = ndata(:,12);
pbfaz = ndata(:,13);

L = sum(~isnan(pbfat))



%% Open the file with new types (5 types as of Nov 02,2014)

%   0 - go to next (Can't determine the type) (type 0)
%   1 - Catogory 1 - Clean bipolar pulses
%   2 - Catogory 2 - Bipolar w extra pulses after overshoot peak
%   3 - Catogory 3 - extra pulses befor and afer overshoot
%   4 - Catogory 4 - Bipolar with extra pulse between first peak and zero cross
%   5 - Catogory 5 - extra pulses before leading peak

fID = fopen('C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\Types-20141030\20110814-type_info.txt','r');
typeD = textscan(fID,'%f %f');
type = typeD{2};
fclose(fID);



%%



%return

% % Make a new base folder to save files
% if savePlots
%     dstr = datestr(now,'yyyy-mm-dd');
%     bf = [bf dstr '\'];
%
%     if ~exist(bf,'dir')
%         mkdir(bf)
%     end
% end
%
%
%% Time vs Altitude
% figure
% plot(pbfat,pbfaz/1000,'ro','markerfacecolor','r','markersize',2)
% ylim([0,ceil(max(pbfaz)/1000)])
% ylabel('Altitude (km)')
% xlabel('Seconds from mid night')
% legend([num2str(L) ' points'])
% if savePlots
%     saveas(gcf,[bf 'z-t.fig'])
%     disp('ssdf')
% end
%
% fprintf('Average altitude = %0.1f\n',nanmean(pbfaz/1000))

%% Altitude Histogram

% % Get data
% xlRange = 'A1:D243';
%
% [ndata, txtdata, data] = xlsread('C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\NBP-xyz.xlsx', 1, xlRange);
%
% %pbfat = ndata(:,10);
% %pbfax = ndata(:,11);
% %pbfay = ndata(:,12);
% pbfaz = ndata(:,4);
%
%
%
% % plot data
% figure
% hold all
% [counts,bins] = hist(pbfaz/1000,(500:1000:20000)/1000)
% %barh(bins,counts)
% ylabel('Altitude (km)')
% xlabel('Number of NBPs')
% % if savePlots
% %     saveas(gcf,[bf 'Z-hist.fig'])
% % end
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
%
% x = x*max(yy)/max(x);
%
% ylim([0 20])
% xlim([0 70])
%
%
% plot(x,y,'r-','Color',[0.7 0.74 0.71],'LineWidth',2)
% plot(yy,xx,'-r','LineWidth',2);
% [mm ind] = max(yy)
% xx(ind)
%
% box on
% legend('\it{Smith et al. 2004}','Current study')
% pH = plot(counts,bins,'ro');

%% XY plot
% figure
% hold all
% scatter(pbfax/1000,pbfay/1000,20,pbfat,'fill')
% florida_map
% cbh = colorbar;
% ylabel(cbh, 'Seconds from midnight')
% box on
% xlabel('East (km)')
% ylabel('North (km)')
% daspect([1 1 1])
% title([num2str(L) ' points'])
%
% if savePlots
%     saveas(gcf,[bf 'xy.fig'])
% end

%% XY plot according to flash types (Original 8 types)
% marker = {'ro','ko','bo','go','mo','yo','co','ks'};
% mfc = {'r','k','b','g','m','y','c','k'};
% 
% figure
% hold all
% 
% for i = 0:5
%     indx = find(i == type)
%     length(indx)
%     plot(pbfax(indx)/1000,pbfay(indx)/1000,marker{i+1},'markerfacecolor',mfc{i+1},'markersize',2)
% end
% 
% 
% 
% florida_map
% cbh = colorbar;
% ylabel(cbh, 'Seconds from midnight')
% box on
% xlabel('East (km)')
% ylabel('North (km)')
% daspect([1 1 1])
% legend({'0-not clear','A-Clean','B-pulses after neg peak',...
%     'C-pulses before and after neg peak','D-pulses before zero cross',...
%     'E-Pulses in leading peak'})
%% XY plot according to flash types (Later 4 types)
% figure
% hold all
%
% indx = find(ndata(:,16) == 1); plot(pbfax(indx)/1000,pbfay(indx)/1000,'ro','markerfacecolor','r','markersize',2);
% indx = find(ndata(:,16) == 6); plot(pbfax(indx)/1000,pbfay(indx)/1000,'ko','markerfacecolor','k','markersize',2);
% indx = find(ndata(:,16) == 2); plot(pbfax(indx)/1000,pbfay(indx)/1000,'bo','markerfacecolor','b','markersize',2);
% indx = find(ndata(:,16) == 3); plot(pbfax(indx)/1000,pbfay(indx)/1000,'go','markerfacecolor','g','markersize',2);
% plot(sns_x,sns_y,'mp','markerfacecolor','m','markersize',5)
% indx = find(ndata(:,16) == 4); plot(pbfax(indx)/1000,pbfay(indx)/1000,'go','markerfacecolor','g','markersize',2);
% indx = find(ndata(:,16) == 5); plot(pbfax(indx)/1000,pbfay(indx)/1000,'go','markerfacecolor','g','markersize',2);
% indx = find(ndata(:,16) == 7); plot(pbfax(indx)/1000,pbfay(indx)/1000,'go','markerfacecolor','g','markersize',2);
% indx = find(ndata(:,16) == 8); plot(pbfax(indx)/1000,pbfay(indx)/1000,'ko','markerfacecolor','k','markersize',2);
%
% florida_map
%
% legend('Type A','Type B','Type C','Noisy','Sensors')
%
% box on
% xlabel('East (km)')
% ylabel('North (km)')
% daspect([1 1 1])

%% Isolation
% figure
% iso = ndata(:,22);
% [y x] = hist(abs(iso),[0.5:1:9.5]);
% %[AX h1 h2] = plotyy(NaN,NaN,NaN,NaN);
% barh(x,y,1)
% set(gca,'XDir','reverse');
%
%
% legend([num2str(L) ' points'])
% xlabel('dT (s)')
% ylabel('Number of NBPs')
% title('NBP isolation 1')
%
% if savePlots
%     saveas(gcf,[bf 'pre-isol-0.1s.fig'])
% end
%
% figure
% iso2 = ndata(:,29);
% hist(iso2,[0.5:1:9.5])
% legend([num2str(L) ' points'])
% xlabel('dT (s)')
% ylabel('Number of NBPs')
% title('NBP isolation 2')
% set(gca,'YDir','reverse');
% if savePlots
%     saveas(gcf,[bf 'post-isol-0.1s.fig'])
% end
%
% figure
% plot(abs(iso),iso2,'ro','markerfacecolor','r','markersize',2)

% figure
% iso = ndata(:,22);
% hist(iso,[-0.010:0.001:0])
% legend([num2str(L) ' points'])
% xlabel('dT (s)')
% ylabel('Number of NBPs')
% title('NBP isolation 1')
% if savePlots
%     saveas(gcf,[bf 'pre-isol-0.001s.fig'])
% end
%
% figure
% iso = ndata(:,29);
% hist(iso,[0:0.001:0.010])
% legend([num2str(L) ' points'])
% xlabel('dT (s)')
% ylabel('Number of NBPs')
% title('NBP isolation 2')
% if savePlots
%     saveas(gcf,[bf 'post-isol-0.001s.fig'])
% end

%% xy plot with IDs

% % file name for data
% fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814-NBP_info2.xlsx';
%
% % base folder for saving plots
% bf = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\Results\';
%
% % Save plots?
% savePlots = 1;
%
% sheet = 1;
% xlRange = 'A4:AD307';
%
% ndata = xlsread(fn, sheet, xlRange);
%
% pbfat = ndata(:,10);
% pbfax = ndata(:,11)/1000;
% pbfay = ndata(:,12)/1000;
% pbfaz = ndata(:,13)/1000;
%
% cm = colormap;
% Lcm = length(cm);
% tmin = min(pbfat);
% tmax = max(pbfat);
% cind = round(1+ (pbfat-tmin)/(tmax-tmin)*(Lcm - 1));
%
% figure
% hold all
% for i = 1:length(pbfax)
%
%     if ~isnan(pbfat(i))
%         plot(pbfax(i),pbfay(i),'o','markeredgecolor',cm(cind(i),:),'markerfacecolor',cm(cind(i),:))
%         text(pbfax(i),pbfay(i),num2str(i))
%     end
% end
%
% florida_map
% daspect([1 1 1])
% box on
% xlabel('East (km)')
% ylabel('North (km)')
% colorbar

%% Isolation both sided
figure
hold all
ph1 = patch([.66 11 11 .66],[-1 -1 11 11],'g');
ph2 = patch([-1 11 11 -1],[.66 .66 11 11],'b');

 alpha(ph2,.3)
 alpha(ph1,0.3)


plot(ndata(:,29),abs(ndata(:,22)),'ro','markerfacecolor','r','markersize',4)
xlabel('Post pulse time (s)')
ylabel('Pre pulse time (s)')
axis equal
box on
legend(sprintf('%i Samples',L))
if savePlots
    saveas(gcf,[bf 'pre-post-isolation.fig'])
end

xlim([-.5 10.5])
ylim([-.5 10.5])

fprintf('Mean time that pulse happen before = (%0.3f+/-%0.3f)s\n',nanmean(ndata(:,22)),mad(ndata(:,22)))
fprintf('Mean time that pulse happen after  = (%0.3f+/-%0.3f)s\n',nanmean(ndata(:,29)),mad(ndata(:,29)))

% pecent of pulses that started flash
iso = cell2mat(data(:,29));
iso2 = cell2mat(data(:,22));
inds = isnan(iso);
iso(inds) = [];
iso2(inds) = [];

inds = isnan(iso2);
iso(inds) = [];
iso2(inds) = [];



ind = intersect(find(iso >= 0.66), find(abs(iso2) >= 0.66));
fprintf('Truly isolated = %0.1f%%\n',length(ind)/length(iso)*100);

ind = intersect(find(iso > 0.66), find(abs(iso2) < 0.66));
fprintf('After IC = %0.1f%%\n',length(ind)/length(iso)*100);

ind = intersect(find(iso < 0.66), find(abs(iso2) > 0.66));
fprintf('Before IC = %0.1f%%\n',length(ind)/length(iso)*100);

ind = intersect(find(iso < 0.66), find(abs(iso2) < 0.66));
fprintf('Within IC = %0.1f%%\n',length(ind)/length(iso)*100);


iso(ind) = nan;
iso2(ind) = nan;


ind = find(iso < 0.1);
length(ind)/L*100


ind = find(iso2 > -0.1);
length(ind)/L*100


%% Pulse Types
% figure
% hist(ndata(:,16),1:1:8)
% str = sprintf('1: Clean bipolar\n2: Like a positive RS\n3: Noisy front\n4: Overall noise \n5: Noisy back\n6: nice pulse with bouncing in the tail\n7: Noise between peaks\n8: Noisy pulse with bounce\n');
% legend(str,'Location','NorthWest')
% xlabel('Pulse type')
% ylabel('Number of NBPs')
% if savePlots
%     saveas(gcf,[bf 'Pulse_types.fig'])
% end
% L = sum(~isnan(ndata(:,16)));
% for i = 1:8
%     fprintf('%i\t%i\t%.1f\n',i,length(find(ndata(:,16) == i)),length(find(ndata(:,16) == i))/L*100)
% end
% 
% % isolated and non isolated percentages
% iso = cell2mat(data(:,29));
% iso2 = cell2mat(data(:,22));
% 
% for i = 1:8
%     inds = find(ndata(:,16) == i);
%     inds2 = intersect(find(iso(inds) >= 0.66), find(abs(iso2(inds)) >= 0.66));
%     fprintf('Type %i: iso = %3.1f\t non iso = %3.1f\n',i,length(inds2)/L*100,(length(inds) - length(inds2))/L*100)
% end

%% Pulse Types new (5types as of Nov 02,2014)
% figure
% hist(type,0:1:5)
% xlabel('Pulse type')
% ylabel('Number of NBPs')
% 
% 
% str = sprintf('0-not clear\nA-Clean\nB-pulses after neg peak\nC-pulses before and after neg peak\nD-pulses before zero cross\nE-Pulses in leading peak');
% legend(str)

%% Peak currents
% figure
% hist(abs(ndata(:,15)),0:10:120)
% title('Peak Current Distributions')
% xlabel('Peak currents (kA)')
% ylabel('Number of NBPs')
% legend(sprintf('%i NBPs',L))
% tools2fig
% 
% if savePlots
%     saveas(gcf,[bf 'Ip-histo.fig'])
% end
% 
% fprintf('***** Peak Currents *******')
% fprintf('Max = %.1f kA\n',nanmax(abs(ndata(:,15))))
% fprintf('Min = %.1f kA\n',nanmin(abs(ndata(:,15))))
% fprintf('Avg = %.1f kA\n',nanmean(abs(ndata(:,15))))
% fprintf('STD = %.1f kA\n',nanstd(abs(ndata(:,15))))
% fprintf('MOD = %.1f kA\n',mode(abs(ndata(:,15))))


%% Rise time and fall times
% fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\NBP_physical_properties.txt';
% 
% fID = fopen(fn);
% data = textscan(fID,'%f %f %f %f %f %f %f %f %f %f');
% fclose(fID);
% 
% indx = data{1};
% rt = data{4} - data{3};
% ft = data{6} - data{4};
% tt = data{6} - data{3};
% 
% 
% L = sum(~isnan(rt));
% [mmm mn] = max(rt)
% sprintf('%0.7f\n',data{4}(mn))
% nanmean(rt)
% 
% figure
% hist(rt*1e6,[0:20])
% xlabel('Rise Time (\mus)')
% ylabel('Number of pulses')
% legend(sprintf('N = %i',L))
% xlim([0 20])
% title('Rise Time')
% 
% 
% figure
% hist(ft*1e6,[20:5:120])
% xlabel('Fall Time (\mus)')
% ylabel('Number of pulses')
% legend(sprintf('N = %i',L))
% title('Fall Time')
% 
% 
% figure
% hist(tt*1e6,[20:5:120])
% xlabel('Total Time (\mus)')
% ylabel('Number of pulses')
% legend(sprintf('N = %i',L))
% title('Total Time')
% xlim([20 110])
% nanmean(tt)

%% Rise time and fall times of different types
% clc
% fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\NBP_physical_properties-enhansed2.txt';
% 
% fID = fopen(fn);
% data = textscan(fID,'%f %f %f %f %f %f %f %f %f %f %f %f');
% fclose(fID);
% 
% ind = data{1};
% 
% rt = data{11};
% tt = data{12};
% 
% 
% types = ndata(:,16);
% 
% indx = [find(types == 1)];
% %indx = [find(types == 1) ; find(types == 5); find(types == 7)];
% rt1 = rt(indx);
% tt1 = tt(indx);
% 
% figure
% subplot(4,2,1)
% hist(rt1,1:18)
% title('Clean NBPs')
% subplot(4,2,2)
% hist(tt1,10:10:120)
% title('Clean NBPs')
% 
% 
% indx = [find(types == 6) ; find(types == 8)];
% rt1 = rt(indx);
% tt1 = tt(indx);
% 
% subplot(4,2,3)
% hist(rt1,1:18)
% title('Bouncing NBPs')
% subplot(4,2,4)
% hist(tt1,10:10:120)
% title('Bouncing NBPs')
% 
% %indx = [find(types == 2); find(types == 3); find(types == 4)];
% indx = [find(types == 2)];
% rt1 = rt(indx);
% tt1 = tt(indx);
% 
% subplot(4,2,5)
% hist(rt1,1:18)
% title('RS like NBPs')
% subplot(4,2,6)
% hist(tt1,10:10:120)
% title('RS like NBPs')
% 
% 
% subplot(4,2,7)
% hist(rt,1:18)
% title('ALL NBPs')
% subplot(4,2,8)
% hist(tt,10:10:120)
% title('ALL NBPs')
% 
% L = sum(~isnan(tt1))
% tools2fig


%% Rise time and fall times of different catogories

% % File name for rice time fall time data
% fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\NBP_physical_properties.txt';
% 
% fID = fopen(fn);
% data = textscan(fID,'%f %f %f %f %f %f %f %f %f %f');
% fclose(fID);
% 
% rts = data{4} - data{3};
% fts = data{6} - data{4};
% tts = data{6} - data{3};
% 
% 
% % file name for NBP data
% fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814-NBP_info2.xlsx';
% sheet = 1;
% xlRange = 'A3:AD307';
% ndata = xlsread(fn, sheet, xlRange);
% 
% indx = ndata(:,1);
% types = ndata(:,16);
% 
% % To put in the title
% pulse_kinds = { ...
%     '1: Clean bipolar'
%     '2: Like a positive RS'
%     '3: Noisy front'
%     '4: Overall noise'
%     '5: Noisy back'
%     '6: nice pulse with bouncing in the tail'
%     '7: Noise between peaks'
%     '8: Noisy pulse with bounce'};
% 
% for type = 1:8
%     inds = find(types == type);
%     %indx(inds)
% 
%     fs = [];
%     amps = [];
%     coudnt_get = 0;
% 
%     rt = rts(indx(inds));
%     ft = fts(indx(inds));
%     tt = tts(indx(inds));
% 
%     L = length(rt);
% 
%     % Plot rise times
%     figure
%     hist(rt*1e6,0:20)
%     xlabel('Rise Time (\mus)')
%     ylabel('Number of pulses')
%     legend(sprintf('N = %i',L))
%     %xlim([0 16])
%     title([pulse_kinds{type} ': Rise Time'])
% 
% %     % Plot fall times
% %     figure
% %     hist(ft*1e6,[20:5:120])
% %     xlabel('Fall Time (\mus)')
% %     ylabel('Number of pulses')
% %     legend(sprintf('N = %i',L))
% %     title([pulse_kinds{type} ': Fall Time'])
% %
% %     % Plot total times
% %     figure
% %     hist(tt*1e6,[20:5:120])
% %     xlabel('Total Time (\mus)')
% %     ylabel('Number of pulses')
% %     legend(sprintf('N = %i',L))
% %     title([pulse_kinds{type} ': Total Time'])
% 
% 
% end


%% Aaverage Power calculations 1
% Please use average power calculations 2 witch will do the exactly the
% same thing but faster using .mat data instead of .fig data.
%
% % folder containing power curves
% dn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814_plots\MeanPowers\';
% 
% fns = dir([dn '*.fig']);
% 
% amps = [];
% fs  = [];
% 
% coudnt_get = 0;
% avgPs = nan(1,350);
% for i = 1:length(fns)
% 
%     fgh = open([dn fns(i).name]);
% 
%     % Get data
%     axesObjs = get(fgh, 'Children');  %axes handles
%     dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
%     %get(dataObjs{3},'Type')
%     f = get(dataObjs{3}(1), 'XData') ;
%     amp = get(dataObjs{3}(1), 'YData');
% 
%         % Get the power values from the title
%     tit = get(get(gca,'Title'),'String');
%     power = strrep(tit(2,:), '<P_{tot}> =', '');
%     avgPs(i) = str2double(strrep(power, 'W', ''));
% 
%     % save data for future use
%     %b.f = f;
%     %b.amp = amp;
%     %b.avgPower = avgPs(i);
% 
%     %[PATH,NAME,EXT] = fileparts([dn fns(i).name]);
%     %save([PATH NAME '.mat'],'-Struct','b')
% 
% 
%     % Normalize
%     amp = amp/max(amp);
% 
%     try
%         fs = [fs ;f];
%         amps = [amps; amp];
%     catch
%         coudnt_get = coudnt_get + 1
%     end
% 
%     delete(fgh)
% end
% 
% figure(500)
% hold all
% plot(fs',amps')
% 
% figure
% hold all
% box on
% plot(fs(1,:),mean(amps), 'LineWidth',2)
% mn = mean(amps);
% er = std(amps);
% ph = patch([fs(1,2:end) fliplr(fs(1,2:end))],[mn(2:end)+er(2:end) fliplr(mn(2:end))-fliplr(er(2:end))],'r','edgecolor','none');
% alpha(ph,0.2)
% 
% 
% figure
% avgPs(isnan(avgPs)) = [];
% %let's save power data
% a.avgPs = avgPs;
% save([dn 'MeanPower.mat'],'-Struct','a')
% hist(avgPs)

%% Aaverage Power calculations 2
% % folder containing power mat files
% dn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814_plots\MeanPowers\';
%
% fns = dir([dn '*.mat']);
%
% amps = [];
% fs  = [];
%
% coudnt_get = 0;
% for i = 1:length(fns)
%
%     b = open([dn fns(i).name]);
%
%     % save data for future use
%     f = b.f;
%     amp = b.amp;
%     % Normalize
%     amp = amp/max(amp);
%
%     try
%         fs = [fs ;f];
%         amps = [amps; amp];
%     catch
%         coudnt_get = coudnt_get + 1
%     end
% end
%
%
%
% figure
% hold all
% mn = mean(amps);
% er = std(amps);
% fill([fs(1,2:end) fliplr(fs(1,2:end))],[mn(2:end)+er(2:end) fliplr(mn(2:end))-fliplr(er(2:end))],[240/250 128/250 128/250],'edgecolor','none');
% %alpha(ph,0.2)
% plot(fs(1,:),mean(amps), 'LineWidth',2)
% xlabel('log_{10}(f in Hz)')
% ylabel('<Normalized power density>')
% box on
%
% figure
% hold all
% plot(fs',amps')
% box on
% xlabel('log_{10}(f in Hz)')
% ylabel('Normalized power density')

%% Mean Power distribution
% fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814_plots\MeanPower.mat';
%
% a = open(fn);
% min(a.avgPs)
% figure
% hold all
% x = log10(a.avgPs);
%
% % Histogram
% [ys xs] = hist(x);
% h = bar(xs,ys,1);
% set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
%
% % Gausian fit
% A = 37.26; %+-5.296
% B = 7.021; %+-0.064
% C = -0.7862; %+-0.15
% D = 6.067; % +-4.48   RMSE = 5.424
%
% xs = min(x):0.1:max(x);
% ys = A*exp(-(xs-B).^2/(C^2))+D;
% plot(xs,ys,'r','LineWidth',2)
%
%
% min(a.avgPs)
% max(a.avgPs)
%
%
%
% box on
% xlabel('log_{10}(<Power in watts>)')
% ylabel('Number of NBPs')
% legend(['\mu = ' sprintf('%0.2f',B) '    \sigma = ' sprintf('%0.2f',C/sqrt(2))])

%% Power distributions according to the types

% % file name for NBP data
% fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814-NBP_info2.xlsx';
% sheet = 1;
% xlRange = 'A3:AD307';
% ndata = xlsread(fn, sheet, xlRange);
% 
% % directory name for power data
% %dn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814_plots\MeanPowers\';
% dn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\PowerCurves5MHz-2\';
% 
% indx = ndata(:,1);
% types = ndata(:,16);
% 
% % To put in the title
% pulse_kinds = { ...
%     '1: Clean bipolar'
%     '2: Like a positive RS'
%     '3: Noisy front'
%     '4: Overall noise'
%     '5: Noisy back'
%     '6: nice pulse with bouncing in the tail'
%     '7: Noise between peaks'
%     '8: Noisy pulse with bounce'};
% 
% for type = 1:8
%     %inds = [find(types == 6); find(types == 8)] ;
%     inds = find(types == type);
%     indx(inds)
% 
%     fs = [];
%     amps = [];
%     coudnt_get = 0;
% 
%     for i = 1:length(inds)
% 
%         % Power file name
%         %pfn = sprintf('%sMeanPowersMeanPower-%3.3i.mat',dn,indx(i));
%         pfn = sprintf('%s%4.4i-power_dist.fig',dn,indx(inds(i)));
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
%             f = d.fx;
%             amp = d.dataMeanP;
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
%     figure
%     hold all
%     mn = mean(amps);
%     er = std(amps);
%     fh = fill([fs(1,2:end) fliplr(fs(1,2:end))],[mn(2:end)+er(2:end) fliplr(mn(2:end))-fliplr(er(2:end))],[240/250 128/250 128/250],'edgecolor','none');
%     set(get(get(fh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
%     %alpha(ph,0.2)
%     plot(fs(1,:),mean(amps), 'LineWidth',2)
%     xlabel('log_{10}(f in Hz)')
%     ylabel('<Normalized power density>')
%     legend(sprintf('Avg curve for %i NBPs',length(inds)))
%     title(pulse_kinds{type})
%     box on
%     
% end

%% NBPs basic statistics
% file name for data
% fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814-NBP_info2.xlsx';
%
% sheet = 1;
% xlRange = 'A3:AD307';
%
% ndata = xlsread(fn, sheet, xlRange);
%
% x = ndata(:,11);
% y = ndata(:,12);
% z = ndata(:,13);
%
% R = sqrt(x.^2+y.^2+z.^2)/1000;
%
% nanmin(R)
% nanmax(R)
% nanmean(R)
%
% % Horisontal locations
% figure
% plot(x/1000,y/1000,'ro','markerfacecolor','r','markersize',2)
% florida_map
% daspect([1 1 1])
% length(~isnan(z))
% xlabel('East (km)')
% ylabel('North (km)')


%% New time comparison
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
% 
% tfn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\Times\20110814-time_info-modified.txt';
% 
% fID = fopen(tfn);
% tdata = textscan(fID,'%f %f %f %f %f %f %f','headerlines',2);
% fclose(fID);
% 
% tdata = cell2mat(tdata);
% 
% types = ndata(:,16);
% 
% %     '1: Clean bipolar'
% %     '2: Like a positive RS'
% %     '3: Noisy front'
% %     '4: Overall noise'
% %     '5: Noisy back'
% %     '6: nice pulse with bouncing in the tail'
% %     '7: Noise between peaks'
% %     '8: Noisy pulse with bounce'};
% %
% 
% 
% figure
% tools2fig
% str1 = {};
% str2 = {};
% str3 = {};
% str4 = {};
% 
% typesChar = ['A','B','C','D','E','F','G','H'];
% 
% 
% 
% for i = 1:8
%     
%     %     if i == 1
%     %         inds = [find(types == 1); find(types == 5); find(types == 7)];
%     %         col = 'r';
%     %     elseif i ==2
%     %         inds = [find(types == 6); find(types == 8)];
%     %         col = 'b';
%     %     elseif i == 3
%     %         inds = [find(types == 2); find(types == 4); find(types == 3)];
%     %         col = 'g';
%     %     elseif i== 4
%     %         inds = 1:length(types);
%     %         col = 'k';
%     %     end
%     
%     
%     
%     
%     switch i
%         case 1
%             inds = find(types == 1);
%             col = 'r';
%         case 2
%             inds = find(types == 6);
%             col = 'g';
%         case 3
%             inds = find(types == 2);
%             col = 'b';
%         case 4
%             inds = find(types == 5);
%             col = 'c';
%         case 5
%             inds = find(types == 8);
%             col = 'm';
%         case 6
%             inds = find(types == 4);
%             col = 'y';
%         case 7
%             inds = find(types == 7);
%             col = 'k';
%         case 8
%             inds = find(types == 3);
%             col = 'k';
%     end
%     
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
%     % 10-90 rise time
%     subplot(2,2,1)
%     hold all
%     [y2, x2] = hist(tdata(inds,3),0.5:1.5:15);
%     xlabel('10-90% Rise time (\mus)')
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
%     str1 = [str1 sprintf('%s (%0.1f+-%0.1f)',typesChar(i),mn,stdd)];
%     legend(str1);
%     
%     % FWHM
%     subplot(2,2,2)
%     hold all
%     [y3, x3] = hist(tdata(inds,7),0.5:10);
%     xlabel('FWHM (\mus)')
%     
%     %spline fit
%     splinefitplot([0 x3],[0 y3],0:.1:10,col)
%     xlim([0 10])
%     %ylim([0 45])
%     box on
%     
%     mn = nanmean(tdata(inds,7));
%     stdd = nanstd(tdata(inds,7));
%     str2 = [str2 sprintf('%s (%0.1f+-%0.1f)',typesChar(i),mn,stdd)];
%     legend(str2);
%     
%     % 0 cross time
%     subplot(2,2,3)
%     hold all
%     [y4, x4] = hist(tdata(inds,4),0.5:4:35);
%     xlabel('Zero Cross Time (\mus)')
%     
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
%     str3 = [str3 sprintf('%s (%0.1f+-%0.1f)',typesChar(i),mn,stdd)];
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
%     str4 = [str4 sprintf('%s (%0.1f+-%0.1f)',typesChar(i),mn,stdd)];
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
% %subplot(2,2,1)
% %legend('A','B','C','D','E','F','G','H')
% str = {};
% % times in one plot
% % 10-90 rise time
% figure
% hold all
% inds = 1:length(types);
% [y2, x2] = hist(tdata(inds,3),0.5:1.5:15);
% xlabel('10-90% Rise time (\mus)')
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
% str = [str sprintf('%s (%0.1f+-%0.1f)',typesChar(i),mn,stdd)];
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
% str = [str sprintf('%s (%0.1f+-%0.1f)',typesChar(i),mn,stdd)];
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
% str = [str sprintf('%s (%0.1f+-%0.1f)',typesChar(i),mn,stdd)];
% legend(str);
% 
% % Negative peak time
% 
% [y5, x5] = hist(tdata(inds,5),0:4:32);
% xlabel('Negative peak Time (\mus)')
% 
% %spline fit
% splinefitplot(x5,y5,0:.1:32,'k')
% xlim([0 35])
% %ylim([0 32])
% box on
% 
% mn = nanmean(tdata(inds,5));
% stdd = nanstd(tdata(inds,5));
% str = [str sprintf('%s (%0.1f+-%0.1f)',typesChar(i),mn,stdd)];
% legend(str);



%% Find NBPs near sensors

% sn = 4;
% range = 10;
%
% r = sqrt((pbfax/1000-sns_x(sn)).^2+(pbfay/1000-sns_y(sn)).^2+pbfaz/1000.^2);
%
% for i = 1:length(pbfat)
%     if r(i) <= range
%         fprintf('%3.3i\t%12.6f\t%5.1f\t%5.1f\t%5.1f\t%5.1f\n',...
%             i,pbfat(i),pbfax(i)/1000,pbfay(i)/1000,pbfaz(i)/1000,r(i))
%     end
% end

%% Find NBPs within reversal distance of sensors

% for i = 1:length(pbfat)
%     % Distance to each sensor
%     d = sqrt((pbfax(i)/1000-sns_x(1:10)).^2+(pbfay(i)/1000-sns_y(1:10)).^2+pbfaz(i)/1000.^2);
%
%     % reversal distance (aproximately)
%     rd = sqrt(2)*pbfaz(i)/1000;
%     inds = find(d < rd);
%
%     if ~isempty(inds)
%         sensors = '';
%         for j = 1:length(inds)
%             sensors = sprintf('%s, %s(%0.1f)',sensors,sns_n{inds(j)},d(inds(j)));
%         end
%         sensors = sensors(3:end);
%
%         fprintf('%3.3i\t%12.6f\t%5.1f\t%5.1f\t%5.1f\t%5.1f\t%s\n',...
%              i,pbfat(i),pbfax(i)/1000,pbfay(i)/1000,pbfaz(i)/1000,rd,sensors)
%     end
%
%
% end


%% FUNCTIONS BELOW % DO NOT comment%

function splinefitplot(x,y,xs,col)

%spline fit
cs = csape(x,y);
yy = ppval(cs,xs);
plot(xs,yy,[col '-'],'linewidth',2)
H = plot(x(2:end),y(2:end),[col 'o']);
set(get(get(H,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
box on

