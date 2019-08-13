function Ip_vs_Ep_best_fit_line

% Draw the line from Ip_vs_Ep calibration and brush obious outlier data points and
% then run this progam.

h = findobj(gcf,'Type','line');
x = get(h,'Xdata');
y = get(h,'Ydata');

x = x(~isnan(y));
y = y(~isnan(y));

length(x)
length(y)

h=get(gca,'Title');
tit=get(h,'String');



figure
plot(x,y,'ro','markerfacecolor','r')
box on
ylabel('CGLSS I_p (kA)')
xlabel('E_{p(norm)} (V/m)')
cc = corr2(x,y);
title(sprintf('%s CorrCoeff - %0.3f',tit,cc))

% best fit line trough (0,0)
L = length(x);
x = x(1:L);
y = y(1:L);

 a = x(:)\y(:);
 hold all
 mx = max(x);
 plot([0 mx],[0 mx*a],'linewidth',2)
 
 
 % Best fit line
   P = polyfit(x,y,1);
   xx = [0 mx];
   yy = P(1)*xx+P(2);
   plot(xx,yy,'linewidth',2)
 
 legend(sprintf('%i points',L),...
        sprintf('y = %0.1f x',a),...
        sprintf('y = %0.1f x + %0.1f',P(1),P(2)))
 
 
tools2fig





