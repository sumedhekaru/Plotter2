function pbfa_post_processing
sen_set = open('sensor_setting.mat');
%% Error plots
d1 = open('C:\Users\sumedhe\Desktop\Combined_method_200m_res_small_test.mat')
%d2 = open('C:\Users\sumedhe\Desktop\PBFA_errors\reviewer_suggested_bigMap500m.mat');

% figure
% 
% return



figure
subplot(2,2,1)
[ch ch]=contourf(d1.x/1000,d1.y/1000,d1.dz3);
set(ch,'edgecolor','none')
hold all
florida_map
daspect([1 1 1])
c = colorbar;
%caxis([.99 1.02])
caxis([0 1100])
xlabel('East (km)')
ylabel('North (km)')
ylabel(c','Z error (m)')
%ylabel(c','\chi^2')
title('Method 1')
plot(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
        'ro','markerfacecolor','r','MarkerSize',2)
    text(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
         {'K02','K14','K24','BCC','K17','EDW','STC','FLT','OVD','FFI'})
tools2fig
%xlim([-30 10])
%ylim([-20 20])



%figure
subplot(2,2,2)
[ch ch]=contourf(d1.x/1000,d1.y/1000,d1.dz1);
set(ch,'edgecolor','none')
hold all
florida_map
daspect([1 1 1])
c = colorbar;
xlabel('East (km)')
ylabel('North (km)')
ylabel(c','z error (m)')
%ylabel(c','\chi^2')
title('Method 2:Z-using sensors between 6--\infty km')
plot(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
    'ro','markerfacecolor','r','MarkerSize',2)
text(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
    {'K02','K14','K24','BCC','K17','EDW','STC','FLT','OVD','FFI'})
%caxis([0.7  0.75])
caxis([0 1100])
tools2fig
%xlim([-30 10])
%ylim([-20 20])

% %figure
% subplot(2,3,63)
% [ch ch]=contourf(d1.x/1000,d1.y/1000,d1.dz);
% set(ch,'edgecolor','none')
% hold all
% florida_map
% daspect([1 1 1])
% c = colorbar;
% xlabel('East (km)')
% ylabel('North (km)')
% ylabel(c','z error (m)')
% %ylabel(c','\chi^2')
% title('PBFA Combined method')
% plot(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
%     'ro','markerfacecolor','r','MarkerSize',2)
% text(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
%     {'K02','K14','K24','BCC','K17','EDW','STC','FLT','OVD','FFI'})
% %caxis([0.5 1.5])
% caxis([0 5000])
% tools2fig
% xlim([-30 10])
% ylim([-20 20])

%d1 = open('C:\Users\sumedhe\Desktop\PBFA_errors\Combined_method_500m_res_big.mat');

d1 = open('C:\Users\sumedhe\Desktop\Combined_method_500m_res_large_test.mat')

subplot(2,2,3)
[ch ch]=contourf(d1.x/1000,d1.y/1000,d1.dz3);
set(ch,'edgecolor','none')
hold all
florida_map
daspect([1 1 1])
c = colorbar;
%caxis([.99 1.02])
caxis([0 5000])
xlabel('East (km)')
ylabel('North (km)')
ylabel(c','Z error (m)')
%ylabel(c','\chi^2')
title('New method 1')
plot(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
        'ro','markerfacecolor','r','MarkerSize',2)
    text(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
         {'K02','K14','K24','BCC','K17','EDW','STC','FLT','OVD','FFI'})
tools2fig

%d2 = open('C:\Users\sumedhe\Desktop\PBFA_errors\reviewer_suggested_bigMap500m.mat');


%figure
subplot(2,2,4)
[ch ch]=contourf(d1.x/1000,d1.y/1000,d1.dz1);
set(ch,'edgecolor','none')
hold all
florida_map
daspect([1 1 1])
c = colorbar;
xlabel('East (km)')
ylabel('North (km)')
ylabel(c','z error (m)')
%ylabel(c','\chi^2')
title('Method 2: Z-using sensors between 6--\infty km')
plot(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
    'ro','markerfacecolor','r','MarkerSize',2)
text(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
    {'K02','K14','K24','BCC','K17','EDW','STC','FLT','OVD','FFI'})
%caxis([0.7  0.75])
caxis([0 5000])
tools2fig

if d1.dz2 == d1.dz1
    disp('oh yreah')
end

% %figure
% subplot(2,3,6)
% [ch ch]=contourf(d1.x/1000,d1.y/1000,d1.dz);
% set(ch,'edgecolor','none')
% hold all
% florida_map
% daspect([1 1 1])
% c = colorbar;
% xlabel('East (km)')
% ylabel('North (km)')
% ylabel(c','z error (m)')
% %ylabel(c','\chi^2')
% title('PBFA Combined method')
% plot(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
%     'ro','markerfacecolor','r','MarkerSize',2)
% text(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
%     {'K02','K14','K24','BCC','K17','EDW','STC','FLT','OVD','FFI'})
% %caxis([0.5 1.5])
% caxis([0 5000])
% tools2fig

% figure
% [ch ch]=contourf(d1.x/1000,d1.y/1000,d2.dz - d1.dz);
% set(ch,'edgecolor','none')
% hold all
% florida_map
% daspect([1 1 1])
% c = colorbar;
% xlabel('East (km)')
% ylabel('North (km)')
% ylabel(c','\DeltaZ error')
% title('\DeltaZ_{error} = New - Old')
% plot(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
%     'ro','markerfacecolor','r')
% text(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
%     {'K02','K14','K24','BCC','K17','EDW','STC','FLT','OVD','FFI'})
% caxis([-1200 1200])
% tools2fig
