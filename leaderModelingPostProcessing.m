function leaderModelingPostProcessing

% Data file name
% F02 %
%fn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F02\1stRS_final_using_real_chi_squired\ModelDataOut_model_with_PBFA_10kHz_data.mat';
%fn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F02\20140509-K02\ModelDataOut_.mat';
%fn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F02\20140509-K24\\ModelDataOut_.mat';

% F18 1st RS
fn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F18\F18_final_without_K02\ModelDataOut_model_with_PBFA_10kHz_data.mat';

% F18 4th RS
%fn = 'C:\Users\Sumedhe\Desktop\LeaderModeling2013\F18-4th-RS\ModelDataOut_3rd attempt.mat';


% F01 1st RS
%fn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F01\1stRS\20140430\ModelDataOut_.mat';
%fn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F01\1stRS\Final\ModelDataOut_model_with_PBFA_10kHz_data.mat';


% F01 2nd RS
%fn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F01\2ndRS\20140502-step leader part only\ModelDataOut_.mat';
% fn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F01\2ndRS\20140512-LM8\ModelDataOut_.mat';


% For F01 3rd RS
%fn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F01\3rdRS\ModelDataOut_.mat';
fn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F01\3rdRS\20140514-LM8\ModelDataOut_.mat';


d = open(fn);
L = length(d.tuni);
b = d.b;
sen_set = open('sensor_setting.mat');

eps0 = 8.85418782e-12;

[PATH,NAME] = fileparts(fn);

%% Max and min dEs
for i = 1:4
    fprintf('%0.2e\t',min(abs(d.dEsm2(i,:))))
    fprintf('%0.2e\n',max(abs(d.dEsm2(i,:))))
end

figure
plot(diff(d.tuni))
mean(diff(d.tuni))
return


%% Plot cumulative charge on current plot
% hold all
% AX = plotyy(nan,nan,d.tuni(1:L-1),cumtrapz(d.Qs(1:L-1)));
% linkaxes(AX,'x')
% 
% FOrmatting axis
% set(AX(2),'xtick',[])
% set(AX, 'YColor', [0 0 0])
% assignin('base','AX',AX)
% set(AX(2),'YLim',[-10 5])
% %set(AX(2),'YTick',-10:3:5)
% %ylabel(AX(2),'Charge (C)')
% %set(AX(2),'FontSize',15)
% 
% return

%% Plot Cumilative charge and current in new figure
% figure
% plot(d.tuni(1:L-1),cumtrapz(d.Qs(1:L-1)))
% xlabel('Time (s)')
% ylabel('Cumulative charge (C)')
% tools2fig
% hold all
% % Plot current
% plot(d.tuni(1:L-1),d.Qs(1:L-1)/(20e-6)/1000)
% mean(d.Qs(1:L-1)/(20e-6)/1000)
% 
% % down sample to 10Khz
% tws = min(d.tuni):1/10000:max(d.tuni);
% 
% for i = 1:length(tws)-1
%     ts(i) = (tws(i)+tws(i+1))/2;
%     
%     lol = nnz(d.tuni < tws(i))+1;
%     ul = nnz(d.tuni < tws(i+1));
%     
%     I(i) = sum(d.Qs(lol:ul))*10000/1000;
% end
% 
% plot(ts,I)


%% Getting Charge Transfer just the points I want
% First plot all the HSVP points
% figure
% plot(b.N_x,240-b.N_y,'r.')
% daspect([1 1 1])
% 
% % After plot Obtain the after brushing from the figure
% h = findobj(gca,'Type','line');
% N_x=get(h,'Xdata');
% N_y=get(h,'Ydata');
% figure
% plot(N_x,N_y,'r.')
% daspect([1 1 1])
% 
% % Now find the charge
% newD = unique([N_x' N_y'],'rows');
% 
% q = 0;
% 
% for i=1:length(newD)
%     xinds = find(newD(i,1) == N_x);
%     yinds = find(newD(i,2) == N_y);
%     
%     inds = intersect(xinds,yinds);
%     
%     q = sum(b.qs(inds)) + q;
%     
% end
% 
% q



%% Plot q distribution plot
% figure
% qs = b.qs;
% q_dist_thr = -0.006;
% 
% Nsr = sum(qs<q_dist_thr);
% qs(find(qs < q_dist_thr)) = q_dist_thr;
% qs(find(qs > 0)) = 0;
% qs = abs(qs);
% 
% %figure
% hold all
% % data = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F18\background.mat');
% % imagesc(data.img,[0, 65535])
% scatter(b.N_x,240-b.N_y,4,qs,'filled','Marker','o')
% 
% cbh = colorbar;
% cbhy = get(cbh,'ylabel');
% set(cbhy,'String','dQ','fontsize',10)
% daspect([1 1 1])
% box on
% title(sprintf('Charge distribution along channels (shrinked %i points at %0.4fC )',Nsr,q_dist_thr))
% 
% 
% % Plotting PBFA data
% %b = ldar2scr(b.pbfa(:,6),b.pbfa(:,7),b);
% plot(b.pbfa(:,6),240-b.pbfa(:,7),'k.')
% 
% 
% indx = find(b.pbfa(:,5) == 0);
% b.pbfa(indx,6) = NaN;
% 
% scatter(b.pbfa(:,6),240-b.pbfa(:,7),abs(b.pbfa(:,5)/min(b.pbfa(:,5)))*1000,[1 0 0])
% 
% % Plotting horizontal line
% [scrX1 scrY1] = ldar2scr(b.ref_x,b.ref_y,1380,b)
 
%% changing axis to show meters q distribution plot

% fg = gcf;
% 
% dy = 1000; % grid dx in meters
% yL = ylim;
% 
% [x1 y1 z1] = scr2ldar(b.scr_x,240-yL(1),b);
% [x2 y2 z2] = scr2ldar(b.scr_x,240-yL(2),b);
% n = floor((z2-z1)/dy);
% zticks = 0:dy:n*dy;
% [scrX scrY] = ldar2scr(zeros(size(zticks))+b.ref_x,...
%     zeros(size(zticks))+b.ref_y,zticks,b);
% [x y z] = scr2ldar(scrX,240-scrY,b);
% set(gca,'YTick',240-scrY,'YTickLabel',zticks)
% ylabel('Altitude (m)')
% 
% dx = 500; % grid dx in meters
% xL = xlim;
% 
% [x1 y1] = scr2ldar(0,b.scr_y,b);
% [x2 y2] = scr2ldar(320,b.scr_y,b);
% hdist = sqrt((x1-x2)^2+(y1-y2)^2);
% n = floor(hdist/dx);
% xticklabs = 0:dx:n*dx;
% scrX  = 320/hdist*xticklabs;
% set(gca,'XTick',scrX,'XTickLabel',xticklabs);
% xlabel('Distance (m)')


%% Plotting q distribution on the video frame
qs = b.qs;
q_dist_thr = -0.003;
Nsr = sum(qs<q_dist_thr);
qs(find(qs < q_dist_thr)) = q_dist_thr;
qs(find(qs > 0)) = 0;
qs = abs(qs);

%figure
figure(gcf)
hold all
%data = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F18\background.mat');
%imagesc(data.img,[0, 65535])

% find color of each point
cm = colormap(jet);

colormap(gray); % Put image color map back to gray
%freezeColors



L = length(cm);
cinds =  round(1 + (L-1)*(qs - 0)/(max(qs) - 0));

length(qs)

for i = 1:L
    inds = find(cinds == i);
    plot(b.N_x(inds),b.N_y(inds),'o','markersize',2,'markerfacecolor',cm(i,:),'markeredgecolor',cm(i,:))
end

%scatter(b.N_x,240-b.N_y,4,qs,'filled','Marker','o')

colormap(jet)
cbh = colorbar;
cbhy = get(cbh,'ylabel');
set(cbhy,'String','dQ (mC)','fontsize',10)
v = caxis;
set(cbh,'YTick',v(1):range(v)/5:v(2));
set(cbh,'YTickLabel',0:max(qs)/5*1000:max(qs)*1000)

daspect([1 1 1])
box on
title(sprintf('Charge distribution along channels (shrinked %i points at %0.4fC )',Nsr,q_dist_thr))


% Plotting PBFA data
%b = ldar2scr(b.pbfa(:,6),b.pbfa(:,7),b);
plot(b.pbfa(:,6),b.pbfa(:,7),'k.')
indx = find(b.pbfa(:,5) == 0);
b.pbfa(indx,6) = NaN;
scatter(b.pbfa(:,6),b.pbfa(:,7),abs(b.pbfa(:,5)/min(b.pbfa(:,5)))*1000,[1 0 0])

xlim auto
ylim auto

% Manupulating Y axis to plot real values
dy = 1000; % grid dx in meters
yL = ylim;

[x1 y1 z1] = scr2ldar(b.scr_x,240-yL(1),b);
[x2 y2 z2] = scr2ldar(b.scr_x,240-yL(2),b);
n = floor((z2-z1)/dy);
zticks = 0:dy:n*dy;
[scrX scrY] = ldar2scr(zeros(size(zticks))+b.ref_x,...
    zeros(size(zticks))+b.ref_y,zticks,b);
[x y z] = scr2ldar(scrX,scrY,b);

set(gca,'YTick',fliplr(scrY),'YTickLabel',fliplr(zticks))
ylabel('Altitude (m)')

dx = 500; % grid dx in meters
xL = xlim;

[x1 y1] = scr2ldar(0,b.scr_y,b);
[x2 y2] = scr2ldar(320,b.scr_y,b);
hdist = sqrt((x1-x2)^2+(y1-y2)^2);
n = floor(hdist/dx);
xticklabs = 0:dx:n*dx;
scrX  = 320/hdist*xticklabs;
set(gca,'XTick',scrX,'XTickLabel',xticklabs);
xlabel('Distance (m)')

% Print the max altitude of video frame
fprintf('Max horizontal viewing distance = %0.1f\n',hdist)
[x1 y1 z1] = scr2ldar(b.scr_x,135,b);
fprintf('Top vertical viewing altitude = %0.1f\n',z1)


% Find out charges at PBFA locations
pbfa = b.pbfa;
% 1,2,3,4 t,x,y,z    5 - charges    6,7 - N_x,N_y

% Locations pairs you need to find charges
locPa = [202 -127.1
        227 -43.3 
        202.2   -21.9
         ];
     
[L1 ~] = size(locPa);

for i = 1:L1
    xind = find(round(pbfa(:,6)*10)/10 == locPa(i,1));
    %yind = find(round((240-pbfa(:,7))*10)/10 == locPa(i,2));
    yind = find(round((pbfa(:,7))*10)/10 == locPa(i,2));
    ind = intersect(xind,yind);
    % Print results
    fprintf('%0.1f\t%0.1f\t%0.6f\n',locPa(i,1),locPa(i,2),pbfa(ind,5))  
end


%% Find out charges at all PBFA locations and list them according to their altitudes

% pbfa = b.pbfa;
% % % 1,2,3,4 t,x,y,z    5 - charges    6,7 - N_x,N_y
% 
% fprintf('Nx\t\t Ny \t\t z (km)\t\t charge (C)\n')
% for i = 1:length(pbfa(:,1))
%     fprintf('%5.3i\t%5.3i\t%6.1f\t%6.4f\n',round(pbfa(i,6)),round(pbfa(i,7)),pbfa(i,4),pbfa(i,5))
% end

%% Some statistics
% fprintf('******************************************\ndT = %0.2f ms\n',(max(b.t)-min(b.t))*1000);
% fprintf('Number of frames = %i\n',max(b.frameN) - min(b.frameN)+1)
% fprintf('Number of HSVPs  = %i\n',length(b.frameN));
% [mm ind] = min(b.N_y);
% [x, y, z] = scr2ldar(b.N_x(ind),mm,b);
% fprintf('Max HSVP altitude = %0.3f km\n',z/1000);
% fprintf('Number of IB points = %i\n',sum(~isnan(b.pbfa(:,1))));
% fprintf('Average speed = %0.2e m/s\n',z/(max(b.t)-min(b.t)));
% fprintf('Total charge transfer = %0.3f C\n',-sum(b.qs))
% fprintf('Average Line Charge density = %0.3f mC/m\n',sum(b.qs)/z*1000)
% fprintf('Average current = %0.3f kA\n',sum(b.qs)/(max(b.t)-min(b.t))/1000)




%% FieldChange should produced from RS
% sns = [1 2 3 6];
% % F18 measured RS field changes
% MRS = [-67.9219 -73.76 -101.38 NaN -154.87 NaN NaN NaN NaN];
% 
% fprintf('\nField changes at each sensor\n')
% fprintf('Sensor\tMeasured\tCalculated\t%%Diff\n')
% 
% %1596
% factor = 1.00;
% lol = 1;
% ul  = 1596;
% 
% 
% x = b.x(lol:ul,1);
% y = b.y(lol:ul,1);
% z = b.z(lol:ul,1);
% q = b.qs(lol:ul,1);
% 
% data = [x  y z q];
% data = sortrows(data,3);
% x = data(:,1);
% y = data(:,2);
% z = data(:,3);
% q = data(:,4);
% 
% 
% % It is beleive that all return stroke does is nutralizing charges on the
% % leader. We will check that here. We will make a dipole by placing
% % negative charge on the ground (doesn't make any vertical E change)
% % and equal and positive charge on each PBHSV
% % locations. See what is the field change that produce.
% E1 = MRS - MRS;
% 
% for i = sns
%     E1(i) = sum(1/(2*pi*eps0)* factor* q.* z ./ ...
%         ((x-sen_set.x(i)).^2+(y-sen_set.y(i)).^2+z.^2).^1.5);
% 
% end
% 
% 



% Add linearly increasing charge
% ks = 0.000000000001:0.000000001:0.000002;
% er = size(ks);
% 
% 
% for j = 1:length(ks)
%     er_temp = 0;
%     
%     qn = q + ks(j)*z;
%     
%     for i = 1:length(sns)
%         E = sum(1/(2*pi*eps0)* factor* qn.* z ./ ...
%             ((x-sen_set.x(sns(i))).^2+(y-sen_set.y(sns(i))).^2+z.^2).^1.5);
%        
%         err = ((E - MRS(i))/(E + MRS(i)))^2;
%         
%         if ~isnan(err)
%             er_temp = er_temp + err;
%         end
%         %fprintf('%s\t\t%0.2f\t\t%0.2f\t\t%0.1f%%\n',sen_set.sen_IDs{sns(i)},MRS(i),E,abs(E-MRS(i))/abs(E+MRS(i))*200)
%     end
%     er(j) = er_temp;
% end
%  
% figure
% plot(ks,er)
% [mm ind] = min(er);
% 
% k = ks(ind)
% %k = 0;
% qn = q + k*z;
% for i = 1:length(sns)
%     E = sum(1/(2*pi*eps0)* factor* qn.* z ./ ...
%         ((x-sen_set.x(sns(i))).^2+(y-sen_set.y(sns(i))).^2+z.^2).^1.5);
%     
%     fprintf('%s\t\t%0.2f\t\t%0.2f\t\t%0.1f%%\n',sen_set.sen_IDs{sns(i)},MRS(i),E,abs(E-MRS(i))/abs(E+MRS(i))*200)
% end

 %qs(i) * k* b.zp./((b.xp-sen_set.x(sn(j))).^2+(b.yp-sen_set.y(sn(j))).^2+b.zp.^2).^1.5

%% Q - H
% figure
% %plot(round(b.z/10)*10,b.qs,'r.')
% 
% zs1 = unique(round(b.z/10)*10);
% zs2 = round(b.z/10)*10;
% qs1 = nan(size(zs1));
% 
% for i = 1:length(zs1)
%     ind = find(zs2 == zs1(i));
%     %qs1 = b.qs(ind)
%     qs1(i) = sum((b.qs(ind)));
% end
% hold all
% plot(qs1,zs1,'b.')
% box on
% xlabel('Charge (C)')
% ylabel('Altitude (m)')
% title(NAME,'interpreter','none')
% tools2fig
% 
% bins = 0:100:5000;
% q = zeros(size(bins));
% 
% 
% for i=1:length(b.qs);
%     n = sum(bins <= b.z(i));
%     q(n) = q(n)+b.qs(i);
% end
% 
% figure
% plot(q,bins)
% xlabel('Charge (C)')
% ylabel('Altitude (m)')
% title(NAME,'interpreter','none')
% %grid on
% tools2fig

%% Plot HSVP points on current plot
% t = b.align_time + (b.frameN - b.align_frame)./b.frame_rate;
% 
% [x y z] = scr2ldar(b.N_x,b.N_y,b);
% 
% lg = get(findobj(gcf,'Type','axes','Tag','legend'),'string')
% 
% 
% figure(gcf)
% h=findobj(gcf,'Type','axes');
% lg = get(findobj(gcf,'Type','axes','Tag','legend'),'string');
% %[AX,H1,H2]=plotyy(nan,nan,nan,nan);
% %plot (AX(2),DLS(:,10),DLS(:,8),'ko','MarkerFaceColor','k','MarkerSize',mz)
% plot(h(2),t,z,'go','markerfacecolor','g','markersize',1)
% legend([lg 'PBHSV'])


%% Delta-dE vs t plot
% figure
% hold all
% box on
% [L m] = size(d.E);
% 
% lg = {};
% offsets = [-0.05 0.05 0.1];
% for i = 1:L
%     plot(d.tuni(1:m-10),d.dEsc2(i,1:m-10)+offsets(i));
%     hLine = plot(d.tuni(1:m-10),d.dEsm2(i,1:m-10)+offsets(i),'--');  
%     
%     % Exclude line from legend
%     set(get(get(hLine,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); 
% 
%     lg = [lg sen_set.sen_IDs(d.sn(i))];
% 
% end
% 
% xlabel('Time (s)')
% ylabel('\DeltaE-change (V/m)')
% ylim([-0.12 0.12])
% tools2fig
% legend(lg)


%% Plot dE vs t plots
% figure
% hold all
% box on
% [L m] = size(d.E);
% 
% lg = {};
% offsets = [0 -2 -1];
% for i = 1:L
%     plot(d.tuni(1:m-1),d.E(i,1:m-1)+offsets(i));
%     
%     lol = nnz(d.rd.t(i,:) < d.tuni(1));
%     ul = length(d.rd.t(i,:)) - 5;
%     
%     hLine = plot(d.rd.t(i,lol:ul),d.rd.v(i,lol:ul)+offsets(i),'--');
%     
%     % Exclude line from legend
%     set(get(get(hLine,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); 
% 
%     lg = [lg sen_set.sen_IDs(d.sn(i))];
% 
% end
% 
% xlabel('Time (s)')
% ylabel('E-change (V/m)')
% %ylim([-0.12 0.12])
% tools2fig
% legend(lg)

%% ki-squired plot
% figure
% L = length(d.E);
% d.ki_squires(1:L-7)= d.ki_squires(1:L-7)*500;
% plot(d.tuni(1:L-7),d.ki_squires(1:L-7));
% xlabel('Time (s)')
% ylabel('\chi^2 (V^2/m^2)')
% title(sprintf('Max = %0.3f   Min = %0.3f  Avg = %0.8f ',...
%     max(d.ki_squires(1:L-7)),min(d.ki_squires(1:L-7)),nanmean(d.ki_squires(1:L-7))))
% tools2fig


%% xy points on current plot
% t = b.align_time + (b.frameN - b.align_frame)./b.frame_rate;
% [x, y, z] = scr2ldar(b.N_x,b.N_y,b);
% L1 = length(t);
% ts = [t; b.pbfa(:,1)];
% 
% lv = 77358.45; %min(ts);
% uv = 77358.685; %max(ts);
% cm = colormap;
% 
% L = length(cm);
% cv =round( 1 + (L - 1)*(ts - lv)/(uv -lv));
% cv1 = cv(1:L1);
% cv2 = cv(L1+1:end);
% 
% 
% 
% %figure
% hold all
% % %plot (DLS(:,10),DLS(:,8),'ko','MarkerFaceColor','k','MarkerSize',mz)
% % for i = 1:1:length(b.pbfa(:,2))
% %     i
% %     plot(b.pbfa(i,2)/1000,b.pbfa(i,3)/1000,'o','color',cm(cv2(floor(i/2)+1),:),'markerfacecolor',cm(cv2(floor(i/2)+1),:),'MarkerSize',2)
% % end
% % 
% for i = 1:L1
%     plot(x(i)/1000,y(i)/1000,'s','color',cm(cv1(i),:),'markerfacecolor',cm(cv1(i),:),'MarkerSize',2)
% end
% 
% %plot(b.pbfa(:,2)/1000,b.pbfa(:,3)/1000,'bo','markerfacecolor','b')
% %plot(x/1000,y/1000,'s','color',[1.0000    0.8750         0],'markerfacecolor',[1.0000    0.8750         0])
% 
% plot(b.ref_x/1000,b.ref_y/1000,'ks','markerfacecolor','k','markersize',4)



% box on
% %legend('PBFA', 'PBHSV','location','southwest')
% plot(sen_set.x(11)/1000,sen_set.y(11)/1000,'rp','markerfacecolor','r')
% text(sen_set.x(11)/1000+0.1,sen_set.y(11)/1000,'FFI')
% 
% plot(sen_set.x(6)/1000,sen_set.y(6)/1000,'rp','markerfacecolor','r')
% text(sen_set.x(6)/1000+0.1,sen_set.y(6)/1000,'K17')
% 
% plot(sen_set.x(2)/1000,sen_set.y(2)/1000,'rp','markerfacecolor','r')
% text(sen_set.x(2)/1000+0.1,sen_set.y(2)/1000,'K14')
% 
% plot(sen_set.x(3)/1000,sen_set.y(3)/1000,'rp','markerfacecolor','r')
% text(sen_set.x(3)/1000+0.1,sen_set.y(3)/1000,'K24')
% 
% daspect([1 1 1])
% %xlim([-22 -11])
% %ylim([-0 3])
% florida_map
% xlabel('East (km)')
% ylabel('North (km)')
% tools2fig

%% Vertical rho distribution
% d
% d.b
% 
% figure
% plot(b.N_x,240-b.N_y,'ro');
% q = nan(1,240);
% for i = 1:240
%     inds = find(240-b.N_y == i);
%     q(i) = sum(b.qs(inds))
%     
% end
% 
% figure
% plot(q,1:240);

%% Estimate path length and thrue line charge density
% %lineD = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F02\leader-path-data.mat');
% %lineD = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F18\F18-leader-path-data.mat');
% %lineD =open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F01\1stRS\\F01-1st-leader-path-data.mat');
% %lineD = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F01\2ndRS\F01-2nd-leader-path-data.mat');
% lineD = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F01\3rdRS\\F01-3rd-leader-path-data.mat');
% 

% %RS path
 %lineD = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F02\F02-RS_path.mat');
 %lineD = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F18\F18-RS-path-data.mat');
 %lineD = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F01\1stRS\\F01-1st-RS-path-data.mat');
 %lineD = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F01\2ndRS\F01-2nd-RS-path-data.mat');
 lineD = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F01\3rdRS\\F01-3rd-RS-path-data.mat');

% figure
% plot(lineD.x,lineD.y)
% totL = 0;
% [x y z] = scr2ldar(lineD.x,lineD.y,b);
% 
% for i = 1:length(lineD.x)-1
%     if ~isnan(x(i+1)) && ~isnan(x(i))
%         totL = totL + sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2+(z(i+1)-z(i))^2);
%     end    
% end
% totL
% sum(b.qs)
% lamda = sum(b.qs)/totL
% daspect([1 1 1])
% t = range(b.t)*1000
% v = totL/range(b.t)
%% Generating the movie
% stepOut = d.stepOut;
% 
% L = length(d.tuni);
% 
% fg = figure;
% set(fg,'Position',[100 50 1200 1000])
% 
% 
% subplot(2,2,3)
% hold all
% %xlim([77427.03524  77427.0513])
% %ylim([-350 50])
% box on
% %set(gca,'XTick',77427.03524:0.002:77427.0513)
% %set(gca,'XTickLabel',0:2:16)
% xlabel('Time (ms)')
% ylabel('E-change (V/m)')
% 
% subplot(2,2,1)
% hold all
% %xlim([77427.03524  77427.0513])
% %ylim([-5 1])
% box on
% %set(gca,'XTick',77427.03524:0.002:77427.0513)
% %set(gca,'XTickLabel',0:2:16)
% xlabel('Time (ms)')
% ylabel('\DeltaE-change (V/m)')
% 
% subplot(2,2,[2 4])
% hold all
% xlim([50 200])
% ylim([0 450])
% box on
% daspect([1 1 1])
% %set(gca,'YTick',[31 78.69 126.38 174.07 221.77 269.46 317.16 364.85 412.54])
% %set(gca,'YTickLabel',[0 1000 2000 3000 4000 5000 6000 7000])
% %set(gca,'XTick',[166.9 214.61])
% %set(gca,'XTickLabel',[3500 4500])
% xlabel('Distance (m)')
% ylabel('Altitude (m)')
% 
% cbh = colorbar;
% cbhy = get(cbh,'ylabel');
% set(cbhy,'String','-dQ','fontsize',10)
% caxis([0,6e-3])
% 
% 
% Ecs1 = 0; Ecs2 = 0; Ecs3 = 0; Ecs4 = 0;
% Ems1 = 0; Ems2 = 0; Ems3 = 0; Ems4 = 0;
% dEcs1 = 0; dEcs2 = 0; dEcs3 = 0; dEcs4 = 0;
% dEms1 = 0; dEms2 = 0; dEms3 = 0; dEms4 = 0;
% 
% qs = abs(b.qs);
% %q_dist_thr = -0.002;
% 
% % Nsr = sum(qs<q_dist_thr);
% % qs(find(qs < q_dist_thr)) = q_dist_thr;
% % qs(find(qs > 0)) = 0;
% % qs = abs(qs);
% 
% pbfaQs = zeros(size(b.pbfa(:,6)));
% 
% %find(round(b.pbfa(:,6)*10)/10== 195.6)
% %b.pbfa(69,6) = NaN;
% 
% t1 = d.tuni(1);
% t2 = d.tuni(end);
% i = 2;
% 
% F =  moviein(length((t2-t1)/20e-6));
% 
% 
% 
% while  t1 < t2
% 
%     if d.tuni(i) < t1
%     
%         Ecs1(i) = stepOut(i).Ecs(1);    Ecs2(i) = stepOut(i).Ecs(2);
%         Ecs3(i) = stepOut(i).Ecs(3);    %Ecs4(i) = stepOut(i).Ecs(4);
%         
%         Ems1(i) = stepOut(i).Ems(1);    Ems2(i) = stepOut(i).Ems(2);
%         Ems3(i) = stepOut(i).Ems(3);    %Ems4(i) = stepOut(i).Ems(4);
%         
%         dEcs1(i) = stepOut(i).dEcs(1);    dEcs2(i) = stepOut(i).dEcs(2);
%         dEcs3(i) = stepOut(i).dEcs(3);    %dEcs4(i) = stepOut(i).dEcs(4);
%         
%         dEms1(i) = stepOut(i).dEms(1);    dEms2(i) = stepOut(i).dEms(2);
%         dEms3(i) = stepOut(i).dEms(3);    %dEms4(i) = stepOut(i).dEms(4);
%         
%         
%         
%         subplot(2,2,3)       
%         plot(d.tuni(1:i),Ems1,'r')
%         plot(d.tuni(1:i),Ems2-60,'m')
%         plot(d.tuni(1:i),Ems3-120,'b')
%         %plot(d.tuni(1:i),Ems4-210,'k')
%         
%         plot(d.tuni(1:i),Ecs1,'k--')
%         plot(d.tuni(1:i),Ecs2-60,'b--')
%         plot(d.tuni(1:i),Ecs3-120,'r--')
%         %plot(d.tuni(1:i),Ecs4-210,'m--')
%         
%         subplot(2,2,1)
%         plot(d.tuni(1:i),dEms1,'r')
%         plot(d.tuni(1:i),dEms2-1,'m')
%         plot(d.tuni(1:i),dEms3-2,'b')
%         %plot(d.tuni(1:i),dEms4-3.6,'k')
%         
%         plot(d.tuni(1:i),dEcs1,'k--')
%         plot(d.tuni(1:i),dEcs2-1,'b--')
%         plot(d.tuni(1:i),dEcs3-2,'r--')
%         %plot(d.tuni(1:i),dEcs4-3.6,'m--')
%         
%         
%         subplot(2,2,[2 4])
%         plot(b.pbfa(:,6),240-b.pbfa(:,7),'k.')
%         scatter(b.N_x(stepOut(i).ind),240-b.N_y(stepOut(i).ind),4,qs(stepOut(i).ind),'filled','Marker','o')
%         
%         
%         
%         pbfaQs(stepOut(i).pbfaInd) = pbfaQs(stepOut(i).pbfaInd)+stepOut(i).dQ;
%         pbfaQs2 = pbfaQs;
%         indx = find(pbfaQs == 0);
%         pbfaQs2(indx) = NaN;
%         
%         
%         try; delete(hsc); catch; end;
%         hsc = scatter(b.pbfa(:,6),240-b.pbfa(:,7),abs(pbfaQs2/min(b.pbfa(:,5)))*1000,[1 0 0]);
%                
%         i = i+1;
%     end
%     
%     F(i) = getframe(gcf);
%     
%     pause(0.01)
%     t1 = t1 + 20e-6;
% end
% 
% %M = movie(F);
% movie2avi(F(2:i-1),'C:\Users\Sumedhe\Desktop\leader_animation-10fps.avi',...
%     'compression','none',...
%     'fps',10, ...
%     'quality',100 ...
%    )

function [N_x,N_y]=ldar2scr(x,y,z,b)

f = 8.0e-3;  %Focal length
p = 20.0e-6; %Fixel Size

N_x=b.scr_x-(f*tan(atan((y-b.cam_y)./(x-b.cam_x ))- atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x )) ))/p;

N_y=b.scr_y - (z - b.cam_z).*f./(p.*sqrt((x-b.cam_x ).^2+(y-b.cam_y ).^2 ))...
    ./sqrt(1+(tan(atan((y-b.cam_y)./(x-b.cam_x ))- atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x )))).^2);

function [x y z] = scr2ldar(N_x,N_y,b)


f = 8.0e-3;  %Focal length
p = 20.0e-6; %Fixel Size

r = sqrt((((b.ref_x-b.cam_x)).^2+(b.ref_y-b.cam_y ).^2 ).*(1+(p^2.*(b.scr_x-N_x ).^2)./f^2 ));

theta = atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x ))+atan(p*(b.scr_x-N_x)/f);

x = b.cam_x + r.*cos(theta);
y = b.cam_y + r.*sin(theta);
z = b.cam_z + r.*((b.scr_y-N_y).*p)./f;
