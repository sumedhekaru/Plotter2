function out = best_fit_line(x,y,opt)

% this function was written by sumedhe to find bext fit lines trugh points.
% Minimum required input x,y - coordinates of points
% Opttinal input structure
%       opt.title - title
%       opt.xlabel - xlabel
%       opt.ylabel - ylabel
%       opt.figure - handle to a figure to plot graph. If not given new
%                    figure will be created.
% Output structure
%       out.fig
%       out.P - poly fit parameters
%       out.corrCoeff - corrilation coefficient


    
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
try 
    figure(opt.figure)
catch
    fg = figure;
end     

plot(x,y,'ro','MarkerSize',4,'MarkerFaceColor','r');

hold all
xfit = mx-sx:0.1:Mx+sx;
yfit = polyval(P,xfit);
plot(xfit,yfit,'b-','LineWidth',2);
 grid on

axis([mx-sx Mx+sx my-sy My+sy])

try
    alpha = opt.alpha;
catch
    alpha = 0.05;
end

% Add prediction intervals to the plot.
[Y,DELTA] = polyconf(P,xfit,S,'alpha',alpha );%,'predopt','curve','simopt','on');
plot(xfit,Y+DELTA,'b--');
plot(xfit,Y-DELTA,'b--')

try title(sprintf('%s     CorrCoeff = %0.4f',opt.title,CC(1,2))); end
try; xlabel(opt.xlabel); end
try; ylabel(opt.ylabel); end
 
legend(sprintf('%i points',length(x)) , ...
       sprintf('y = %0.4fx + %0.4f',P(1),P(2)),...
       sprintf('%0.0f%% Confidant',(1-alpha)*100),...
       'Location','SouthEast');
   
out.figure = fg;
out.P = P;
out.corrCoeff = CC(1,2);
   