function read_PBFA_CG_filr

fn = 'C:\Users\sumedhe\Desktop\PBFA vs CGLSS\PBFA_CG_data_20110814_1700-2400_SN8-9.txt';


fID = fopen(fn,'r');
data = textscan(fID,'%f %f %f %f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fID);

% Colums
% 1 - Index 
% 2 - PBFA time     3 - PBFA x      4 - PBFA y      5 - PBFA Ip
% 6 - Number of sensors used for PBFA calculations
% 7 - CGLSS time    8 - CGLSS x     9 - CGLSS y     10 - CGLSS Ip
% 11 - number of sensors used by CGLSS

opt.title = 'X comparison';
opt.xlabel = 'CGLSS x (km)';
opt.ylabel = 'PBFA x (km)';
best_fit_line(data{8}/1000,data{3}/1000,opt)


opt.title = 'Y comparison';
opt.xlabel = 'CGLSS y (km)';
opt.ylabel = 'PBFA y (km)';
best_fit_line(data{9}/1000,data{4}/1000,opt)


opt.title = 'Time from midnight comparison';
opt.xlabel = 'CGLSS t (s)';
opt.ylabel = 'PBFA t (s)';
best_fit_line(data{7},data{2},opt)
   
opt.title = 'Peak current comparision';
opt.xlabel = 'CGLSS Ip(kA)';
opt.ylabel = 'PBFA Ip (kA)';
best_fit_line(abs(data{10}),abs(data{5}),opt)


opt.title = 'Range comparison';
opt.xlabel = 'CGLSS R (km)';
opt.ylabel = 'PBFA R (km)';
best_fit_line(((data{8}.^2+data{9}.^2).^0.5)/1000,...
              ((data{3}.^2+data{4}.^2).^0.5)/1000,opt)




figure
y = ((data{8}.^2+data{9}.^2).^0.5)/1000;
n = ceil(max(y)/5);
hist(y,n)
title('PBFA CG Range histogram 2011-08-14 17:00 - 24:00UT')
xlabel('Number of PBFA CG points')
ylabel('Distance (km)')


figure
r  = ((data{8}.^2+data{9}.^2).^0.5)/1000;
dr = abs(((data{8}.^2+data{9}.^2).^0.5)/1000 - ((data{3}.^2+data{4}.^2).^0.5)/1000);
plot(r,dr,'ro','MarkerSize',4,'MarkerFaceColor','r')
xlabel('R (km)')
ylabel('dR (km)')
title('Range difference between PBFA and CGLSS 2011-08-14 17:00 - 24:00UT')

figure
r  = ((data{8}.^2+data{9}.^2).^0.5)/1000;
dr = abs(((data{8}.^2+data{9}.^2).^0.5)/1000 - ((data{3}.^2+data{4}.^2).^0.5)/1000)./r*100;
plot(r,dr,'ro','MarkerSize',4,'MarkerFaceColor','r')
xlabel('R (km)')
ylabel('dR (%)')
title('Percent Range difference between PBFA and CGLSS 2011-08-14 17:00 - 24:00UT')



figure
n = ceil(max(dr)/1);
hist(dr,n)
title('Range difference between CGLSS and PBFA 2011-08-14   17:00 - 24:00UT')
xlabel('dR (km)')
ylabel('Number of points')



