function MTLE_pulse_length

% x = 6612:1:8000;
% lamda = 100;
% pr = 90; % percent decay distance
% 
% y = exp(-(x-x(1))/lamda);
% ind = sum(y >= (1-pr/100));
% 
% figure
% plot(x,y)
% L = x(ind)-x(1)


%=============================
%% peak current of MTLEI
A = 750
H1 = 6400;
H2 = 7030;
alpha = 100;

Ip = A*exp((H2-H1)/alpha)/1000

A = 420
H1 = 7600;
H2 = 8400;
alpha = 113;

Ip = A*exp((H2-H1)/alpha)/1000
