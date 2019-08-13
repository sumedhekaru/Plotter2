function Ip_vs_Ep_calibration

% To drow the best fit line, plrese refer to Ip_vs_Ep_best_fit_line

output_fn = 'C:\Users\sumedhe\Desktop\Ip_Ep_data_20110814_ch1-1MHz-wlith-PBFAxy-20131121.txt';
fID = fopen(output_fn);
data = textscan(fID,'%f %f %s %f %f %f %f');
fclose(fID);


sID = 'OVD';
L = length(data{1});

Ip = nan(1,L);
Eprn = nan(1,L);
k = 0;
for i=1:L
    str = data{3}(i);
    
    if strcmp(str{1}(1:3),sID) && data{6}(i) > 30
        k = k+1;
        Ip(k) = data{7}(i);
        Eprn(k) = data{5}(i)*(data{6}(i)/100).^1.13;
    end
end

figure
plot(abs(Eprn),abs(Ip),'ro','markerfacecolor','r')
box on
ylabel('CGLSS I_p (kA)')
xlabel('E_{p(norm)} (V/m)')
title(sID)
tools2fig
