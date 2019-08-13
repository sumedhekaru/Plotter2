function peakCurrBestFit(Ip,Ep,R,N,str)

% % Additionally plot N = [these values] will plot
% N = [1 1.3];
%     
y = Ip;

max_corr = 0;

for n = 0.1:0.01:2.0
    
    x = Ep.*R.^n;
    

    CC=abs(corrcoef(x,y));
    %fprintf('%1.2f\t%.4f\n',n,CC(1,2))
    
    if max_corr < CC(1,2)
        max_corr = CC(1,2);
        best_n = n;
               
    end
    
end

% Plot the best fit
plot_data(y,Ep,R,best_n,str)

% Plot additional
if ~isempty(N)    
    for i=1:length(N)          
        plot_data(y,Ep,R,N(i),str)
    end
end
    

   
function plot_data(y,Ep,R,n,str)
    
x = Ep.*R.^n;

%linear fit
[P S] = polyfit(x,y,1); 

% Corrilation coeff
CC=abs(corrcoef(x,y));

x = reshape(x,1,length(x));
y = reshape(y,1,length(y));


% Extent of the data.
mx = min(x);
Mx = max(x);
my = min(y);
My = max(y);

% Scale factors for plotting.
sx = 0.05*(Mx-mx);
sy = 0.05*(My-my);

% Plot the data, the fit, and the roots.
figure
plot(x,y,'ro','MarkerSize',4,'MarkerFaceColor','r');

hold all
xfit = mx-sx:0.01:Mx+sx;
yfit = polyval(P,xfit);
plot(xfit,yfit,'b-','LineWidth',2);
 grid on

axis([mx-sx Mx+sx my-sy My+sy])

alpha = 0.05;

% Add prediction intervals to the plot.
[Y,DELTA] = polyconf(P,xfit,S,'alpha',alpha );%,'predopt','curve','simopt','on');
plot(xfit,Y+DELTA,'b--');
plot(xfit,Y-DELTA,'b--')

title(sprintf('Best fit for %s n = %0.2f     CorrCoeff = %.4f',...
    str,n,CC(1,2)))
xlabel('E_p * r^n')
ylabel('I_p')

legend(sprintf('%i points',length(x)) , ...
       sprintf('y = %0.4fx + %0.4f',P(1),P(2)),...
       sprintf('%0.0f%% Confidant',(1-alpha)*100),...
       'Location','SouthEast');

% figure
% box on; hold all;
% x = Ep.*R.^best_n;
% plot(x,y,'ro','markerfacecolor','r','MarkerSize',3)
% x1 = min(x);
% x2 = max(x);
% y1 = best_m(1)*x1 + best_m(2);
% y2 = best_m(1)*x2 + best_m(2);
% plot([x1 x2],[y1 y2],'LineWidth',2)
% title(sprintf('Best fit for %s n = %0.2f     CorrCoeff = %.4f   y = %0.4fx + %0.4f ',...
%     str,best_n,max_corr,best_m(1),best_m(2)))
% xlabel('E_p * r^n')
% ylabel('I_p')
% % 
% 
% 
% 
% if ~isempty(N)
%     
%     for i=1:length(N)
%         n = N(i);
%         x = Ep.*R.^n;
%         
%         CC=corrcoef(x,y);
%         m = polyfit(x,y,1);
%         
%         
%         figure
%         box on; hold all;
%         plot(x,y,'go','markerfacecolor','g','MarkerSize',3)
%         x1 = min(x);
%         x2 = max(x);
%         y1 = m(1)*x1 + m(2);
%         y2 = m(1)*x2 + m(2);
%         plot([x1 x2],[y1 y2],'LineWidth',2)
%         
%         title(sprintf('Fitted for %s n = %0.2f     CorrCoeff = %.4f   y = %0.4fx + %0.4f ',...
%             str,n,CC(1,2),m(1),m(2)))
%         xlabel('E_p * r^n')
%         ylabel('I_p')
%         
%         
%     end
% end
% 
% 
% 
