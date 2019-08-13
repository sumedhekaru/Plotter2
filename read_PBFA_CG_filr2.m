function read_PBFA_CG_filr2

fn = 'C:\Users\sumedhe\Desktop\2011-08-14-PBFA_CG_data_test_5MHz.txt';

fID = fopen(fn,'r');
data = textscan(fID,'%f %f %f %f %f %f %f %s %f %f %f %f','Headerlines',1);
fclose(fID);

% Colums
% 1 - Index 
% 2 - PBFA time     3 - PBFA x      4 - PBFA y      5 - PBFA z
% 6 - PBFA Ip
% 7 - Number of sensors used for PBFA calculations
% 8 - which sensors used
% 9 - CGLSS time    10 - CGLSS x     11 - CGLSS y     12 - CGLSS Ip
% 13 - number of sensors used by CGLSS

opt.title = 'X comparison';
opt.xlabel = 'CGLSS x (km)';
opt.ylabel = 'PBFA x (km)';

% CGLSS range
data{8} = sqrt(data{10}.^2 + data{11}.^2);

data = cell2mat(data);


% choose data between 0 30km
data = sortrows(data,8);
ul = sum(data(:,8) < 50000); 
newdata = data(1:ul,:);

% Choose atleast 8 sensor solutions
newdata = sortrows(newdata,7);
lol = sum(newdata(:,7) < 8);
newdata = newdata(lol:end,:);
best_fit_line(newdata(:,10)/1000,newdata(:,3)/1000,opt)




opt.title = 'Y comparison';
opt.xlabel = 'CGLSS y (km)';
opt.ylabel = 'PBFA y (km)';
best_fit_line(newdata(:,11)/1000,newdata(:,4)/1000,opt)

figure
hold all
plot(newdata(:,3)/1000,newdata(:,4)/1000,'go','markerfacecolor','g','markersize',3)
plot(newdata(:,10)/1000,newdata(:,11)/1000,'ro','markerfacecolor','r','markersize',3)
legend('PBFA','CGLSS')
florida_map
daspect([1 1 1]);
box on

% plotting the circle
ang=0:0.01:2*pi; 
xp=50*cos(ang);
yp=50*sin(ang);
plot(xp,yp);

%return

% opt.title = 'Time from midnight comparison';
% opt.xlabel = 'CGLSS t (s)';
% opt.ylabel = 'PBFA t (s)';
% best_fit_line(data{7},data{2},opt)
%    
% opt.title = 'Peak current comparision';
% opt.xlabel = 'CGLSS Ip(kA)';
% opt.ylabel = 'PBFA Ip (kA)';
% best_fit_line(abs(data{10}),abs(data{5}),opt)
% 
% 
% opt.title = 'Range comparison';
% opt.xlabel = 'CGLSS R (km)';
% opt.ylabel = 'PBFA R (km)';
% best_fit_line(((data{8}.^2+data{9}.^2).^0.5)/1000,...
%               ((data{3}.^2+data{4}.^2).^0.5)/1000,opt)
% 
% 
% 
% 
% figure
% y = ((data{8}.^2+data{9}.^2).^0.5)/1000;
% n = ceil(max(y)/5);
% hist(y,n)
% title('PBFA CG Range histogram 2011-08-14 17:00 - 24:00UT')
% xlabel('Number of PBFA CG points')
% ylabel('Distance (km)')
% 
% 
% figure
% r  = ((data{8}.^2+data{9}.^2).^0.5)/1000;
% dr = abs(((data{8}.^2+data{9}.^2).^0.5)/1000 - ((data{3}.^2+data{4}.^2).^0.5)/1000);
% plot(r,dr,'ro','MarkerSize',4,'MarkerFaceColor','r')
% xlabel('R (km)')
% ylabel('dR (km)')
% title('Range difference between PBFA and CGLSS 2011-08-14 17:00 - 24:00UT')

% figure
% r  = ((data(:,8).^2+data(:,9).^2).^0.5)/1000;
% dr = abs(((data(:,8).^2+data(:,9).^2).^0.5)/1000 - ((data(:,3).^2+data(:,4).^2).^0.5)/1000)./r*100;
% 
% 
% plot(r,dr,'ro','MarkerSize',4,'MarkerFaceColor','r')
% xlabel('R (km)')
% ylabel('dR (%)')
% title('Percent Range difference between PBFA and CGLSS 2011-08-14 17:00 - 24:00UT')
% 
% 
% 
% figure
% n = ceil(max(dr)/1);
% hist(dr,n)
% title('Range difference between CGLSS and PBFA 2011-08-14   17:00 - 24:00UT')
% xlabel('dR (km)')
% ylabel('Number of points')

% Distance between two points
figure

%dr = sqrt((newdata(:,10)-newdata(:,3)).^2+(newdata(:,11)-newdata(:,4)).^2); %abs(((data{8}.^2+data{9}.^2).^0.5)/1000 - ((data{3}.^2+data{4}.^2).^0.5)/1000);
dr = abs(newdata(:,2) - newdata(:,9))*1e6;
dr(find(dr>33)) = [];
length(dr)
max(dr)
min(dr)
std(dr)
mean(dr)
median(dr)
mode(round(dr))
hist(dr,100)





