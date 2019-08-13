a = open('C:\Users\Sumedhe\Desktop\video_intensity_data2.mat')
intensity = a.intensity;
frames = a.frames;

t = (70754.06032045-323e-6)+(frames+62990-1)*20e-6;

xlimv = xlim;

[AX,H1,H2]=plotyy(nan,nan,nan,nan);
hold(AX(2), 'on')
%plot(AX(2),t,intensity,'r.-')

stairs(AX(2),t,intensity,'linewidth',1)
linkaxes(AX,'x')
set(AX(2),'xtick',[])
set(AX, 'YColor', [0 0 0])
set(AX,'xlim',xlimv)

ylimits = get(AX(1),'YLim');
yinc = (ylimits(2)-ylimits(1))/5;

set(AX(1),'YTick',ylimits(1):yinc:ylimits(2))
set(AX(2),'YTick',0:1/5:1)
ylabel(AX(2),'Cumulative Light Intensity')


% Get current legend