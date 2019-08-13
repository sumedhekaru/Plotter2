function bandwidth2

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

y1 = [    0.9953
    1.005
    1
    1.005
    1
    0.9953
    0.9953
    1.
    0.995
    0.995
    0.993
    0.98
    0.95
    0.94
    0.93
    0.92
    0.93
    0.92
    0.89
    0.88
    0.83
    0.82
    0.79
    0.78
    0.77
    0.73
    0.71
    0.69
    0.67
    0.64
    0.63
    0.54
    0.51
    0.49
    0.38];

x1 = log10(SA1s(:,1));
A = 1.509E-14;
B = 0.9953;
C = -4.761;
x = 3:.01:-log(B/A+1)/C;
y = A*(1-exp(-C*x))+B;

figure
semilogx(SA1s(:,1),y1,'ro','markerfacecolor','r')
xlabel('Frequency (Hz)')
ylabel('Normalized Amplitude')
title ('Frequency Response of 10s Slow Antenna')



hold all
semilogx(10.^x,y,'LineWidth',2)
freq = 2.6121e6;
plot([freq freq],ylim,'--g','LineWidth',2)
legend('Real Data','Best fit',['3dB cut off ' num2str(freq/1e6) 'MHz'],'Location','southwest')
ylim([0,1.05])