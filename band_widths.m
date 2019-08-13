
%% 1s Slow antenna
SA1s =[ ...
    1000 8.2
    2000 8.2
    5000 8.2
    6000 8.2
    10000 8.2
    15000 8.2
    20000 8.2
    30000 8.2
    50000 8.2
    100000 8.0
    200000 8.2
    500000 8.0
    1e6    8.0
    1.1e6  7.8
    1.2e6  7.6
    1.3e6  7.5
    1.4e6  7.4
    1.5e6  7.3
    1.6e6  7.2
    1.7e6  7.2
    2.0e6  6.9
    2.1e6  6.8
    2.2e6  6.6
    2.3e6  6.4
    2.4e6  6.3
    2.5e6  6.0
    2.6e6  5.9
    2.7e6  5.7
    2.8e6  5.5
    2.9e6  5.4
    3.0e6  5.1
    3.2e6  4.7
    3.4e6  4.3
    3.6e6 3.9
    3.8e6 3.7
%     4.0e6 3.66  % out liers
%     4.5e6 3.66
%     5.0e6 3.34
%     6.0e6 2.89
];
%1s slow antenenna was not stable when f>6.0MHz

figure
semilogx(SA1s(:,1),SA1s(:,2)/max(SA1s(:,2)),'ro','markerfacecolor','r')
xlabel('Frequency (Hz)')
ylabel('Normalized Amplitude')
title ('Frequency Response of 1s Slow Antenna')


A = 1.509E-14;
B = 0.9953;
C = -4.761;
x = 3:.01:-log(B/A+1)/C;
y = A*(1-exp(-C*x))+B;

hold all
semilogx(10.^x,y,'LineWidth',2)
freq = 10^(-1/C*log(1-(1/sqrt(2)-B)/A));
plot([freq freq],ylim,'--g','LineWidth',2)
legend('Real Data','Best fit',['3dB cut off ' num2str(freq/1e6) 'MHz'],'Location','southwest')

%% 10s slow ant
SA10s =[ ...
    1000    1.05
    5000    1.06
    10000   1.05
    20e3    1.05
    50e3    1.04
    100e3   1.05
    500e3   1.07
    1e6     1.14
    1.5e6   1.15
    2.0e6   1.1
    3.0e6   1.05
    3.5e6   0.897
%     4.5e6   1.30
%     5.0e6   1.403
%     6.0e6   1.6
%     7.0e6   2.01
    ];


figure
semilogx(SA10s(:,1),SA10s(:,2)/max(SA10s(:,2)),'ro','markerfacecolor','r')
xlabel('Frequency (Hz)')
ylabel('Normalized Amplitude')
title ('Frequency Response of 1s Slow Antenna')


%% 100us fast antenna

FA100 = [ ...
    200 .7
    500 1.45
    700 1.81
    800 1.99
    900 2.17
    1000    2.29
    1.1e3   2.43
    1.2e3   2.57
    1.3e3   2.67
    1.4e3   2.77
    1.5e3   2.87
    1.6e3   2.95
    %1.7e3   3.60
    1.8e3   3.12
    1.9e3   3.20
    2.0e3   3.26
    2.2e3   3.38
    2.4e3   3.48
    2.6e3   3.52
    2.8e3   3.60
    3e3     3.64
    3.5e3   3.76
    4.0e3   3.82
    4.5e3   3.86
    5.0e3   3.88
    10e3    4.02
    50e3    4.14
    100e3   4.14
    200e3   3.98
    300e3   3.78
    400e3   3.50
    500e3   3.26
    600e3   2.97
    700e3   2.81
    800e3   2.55
    900e3   2.35
    1e6 2.19
    1.2e6   1.91
    1.4e6   1.61
    1.6e6   1.57
    1.8e6   1.45
    2.0e6   1.31
    2.5e6   .945
    3.0e6   .719
    3.5e6   .487
%     4.0e6   .487
%     5.0e6   0.518
%     6.0e6   .497
    ];
    
 figure
 x = FA100(:,1);
 y = FA100(:,2)/max(FA100(:,2));
 
semilogx(x,y,'bo','markerfacecolor','r')
xlabel('Frequency (Hz)')
ylabel('Normalized Amplitude')
title ('Frequency Response of 100\mus Fast Antenna')

xx = 2.1:.01:6.71;
yy = spline(log10(x),y,xx);

hold all
semilogx(10.^xx,yy,'LineWidth',2)

dif = abs(yy - 1/sqrt(2));

[m ind1] = min(dif);

f1 = 10^xx(ind1);
dif(ind1) = 100;

[m ind2] = min(dif);
f2 = 10^xx(ind2);

plot([f1 f1],ylim,'--g','LineWidth',2)
plot([f2 f2],ylim,'--m','LineWidth',2)

legend('Real Data','smooth fit',['3dB High ' num2str(f1/1e6) 'MHz'],...
    ['3dB Low ' num2str(f2/1000) 'kHz'],'Location','south')

