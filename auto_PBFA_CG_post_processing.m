function auto_PBFA_CG_post_processing

fn = 'C:\Users\sumedhe\Desktop\PBFA_CG_auto_comp\PBFA_CG_data_TOAM6-new.txt';

fID = fopen(fn);
data = textscan(fID,'%f %f %f %f %f %f %f %s %f %f %f %f %f %f','Headerlines',2);
fclose(fID);

% Columns]
% 1 - Index
% 2,3,4,5 - (t,x,y,z) PBFA
% 6 - [USER, METHOD, NPBFA] 3-digits
% 7 - Ip PBFA
% 8 - sensors used (comma separated
% 9 - PBFA ki squired
% 10,11,12 - (t,x,y) CGLSS
% 13 - Ip CGLSS
% 14 - nCGLSS

L = length(data{7});

figure
plot(abs(data{7}),abs(data{13}),'ro','markerfacecolor','r')
hold all
plot([0,90],[0,90])
xlabel('PBFA Ip')
ylabel('CGLSS Ip')
legend(sprintf('%i Points',L),'y = x Line')
title (' y= mx fit')



%% Calculate Peak Curents accoding to new calibration factors
sen_set = open('sensor_setting.mat');

% Read Hilbert Transform Peak file
fID = fopen('C:\Users\sumedhe\Desktop\PBFA_CG_auto_comp\PBFA_CG_data_TOAM6-new_hilb_peaks.txt');
data2 = textscan(fID,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fID)

data2 = cell2mat(data2);

% y = mx values
% m =  [472.4 387.1 418.7 NaN NaN 411.9 676.4  204.2 165.6 209.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];

% y = mx + c
m = [478.9 389.4  422.2 NaN NaN 415.5 684.8 207.2 166.6 211.3 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];

c = [-0.5  -0.2   -0.3  NaN NaN -0.3  -0.4  -0.5   -0.2 -0.4 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];

Ips = nan(L,1);

for i = 1:L
    r = sqrt((sen_set.x - data{3}(i)).^2 + ...
             (sen_set.y - data{4}(i)).^2 + ...
             (sen_set.z - data{5}(i)).^2)/1000;
         size(data2(i,2:21))
         size(r)
     Ip = m.*data2(i,2:21).*((r/100).^1.13)+c;
     indx = find(r > 30);
     Ip = -nanmean(Ip(indx));
     Ips(i) = Ip;
         
end

figure
hold all
plot(abs(Ips),abs(data{13}),'ro','markerfacecolor','r')
plot([0,90],[0,90])
xlabel('PBFA Ip')
ylabel('CGLSS Ip')
legend(sprintf('%i Points',L),'y = x Line')
title('y = mx + c fit')


