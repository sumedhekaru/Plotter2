
% x = 1:1:5;
% y = [ 136 19 7 172
%     0 10  0 0];       
figure   
bar(y','stacked')
   
 x = 1:1:4;
    y = [ 136 19 7 54
        0 10  0 0];
    tot = sum(y(:));
    figure
    h = bar(x,y','stacked');
    sum(y(:))
    ybuff=4;
    for i=1:length(h)
        XDATA=get(get(h(i),'Children'),'XData');
        YDATA=get(get(h(i),'Children'),'YData');
        for j=1:size(XDATA,2)
            x=XDATA(1,j)+(XDATA(3,j)-XDATA(1,j))/2;
            y=YDATA(2,j)+ybuff;
            YDATA(2,j)/tot;
            t=sprintf('%0.1f%%',YDATA(2,j)/tot*100);
            text(x,y,t,'Color','k','HorizontalAlignment','center')
        end
    end
    ylim([0 185])