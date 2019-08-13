function LDAR_frequency(fn,t1,t2,rc,x0,y0,z0,sen_set)
% Plotting frequency vs. time


%fn='F:/data/ldar2/2010/07/11/ldar2_20101920000.txt';             % File name for idar 

%t1= 0;%7230;                                % Lower time bound
%t2= 1500;%7234;
%rc=50000;                               % Bound radus



[CG,CAL,DLS]=ldarExtract(fn,t1,t2,rc,x0,y0,z0,0);

t=[CG(:,10);DLS(:,10)];
t=sort(t);

f=(t(2:end)-t(1:end-1)).^(-1)/1000;    % frequency in KHz


figure;
plot(t(1:end-1),f')

tit=sprintf('LDAR Frequncy     \n%s-%s-%s     UT: %s - %s',...
    fn(end-31:end-28),fn(end-26:end-25),fn(end-23:end-22),sec2hhmmss(t1),sec2hhmmss(t2));

title(tit)
xlabel('Time (s)')
ylabel('Frequency (kHz)')
box on
grid on

hold on
[y,index]=max(f);
plot(t(index),y,'ro','markerfacecolor','r')
text(t(index),1.05*y,num2str(y))
hold off