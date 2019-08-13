t = 0:0.5e-6:0.1;
y = 10*sin(2*pi*5000*t) + 5*sin(2*pi*500000*t)+ 15*sin(2*pi*250000*t);
figure
plot(t,y)



L = length(t);                     % Length of signal
NFFT = 2^nextpow2(L);
[Pxx, fx] = pwelch(y,[], NFFT/2, NFFT, 2.5e6);
figure
plot(fx,Pxx*1.5)
Ptot = trapz(fx,Pxx)