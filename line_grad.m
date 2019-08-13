function line_grad
h= findobj(gcf,'type','line');
% for i=1:length(h)
%     x= get(h(i),'xdata');
%     y= get(h(i),'ydata');
%     
%     if length(x)==2
%         m= diff(y)/diff(x);
%         fprintf('%0.6f\n',m)
%     end
%     
% end

maxyN = 9875;

for i=1:length(h)
    x= get(h(i),'xdata');
    y= get(h(i),'ydata');
    
    
    if length(x)>100
        
        y = y(1:maxyN);
        x = x(1:maxyN);
        
        figure
        plot(x,y)
        
        
        
        % reduce sampling rate
        nn = 10;
        n = floor(length(y)/nn);
        
        ys = mean(reshape(y(1:nn*n),n,nn));
        xs = mean(reshape(x(1:nn*n),n,nn));
        
        ms = (ys(2:end)-ys(1:end-1))./(xs(2:end)-xs(1:end-1));
        xxs = (xs(1:end-1)+xs(2:end))/2;
        
        hold all
        figure
        plot(xxs,ms)
        
        
    end
    
end