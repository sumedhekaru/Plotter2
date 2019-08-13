function change_color_bar_y_label

h1 = findobj(get(gcf,'Children'),'Tag','Colorbar');

% Change ylabel
h = get(h1,'ylabel');
str = {'Time (ms after 21:29:18.45 UT)'};
set(h,'String',str,'fontsize',10)

% % change ytick positions and names
% set(h1,'YTick',1:63/5:64)
% %yticklabel = {'0' '0.2' '0.4' '0.6' '0.8' '1.0' '1.2' '1.4' '1.6'};
% yticklabel = {'0' '3' '6' '9' '12' '15'};
% %yticklabel = {'0' '0.1' '0.2' '0.3' '0.4' '0.5'};
% set(h1,'YTickLabel',yticklabel)

realYlims = [77358.45 77358.655];
tickVals = [77358.45:0.05:77358.655];
tickLabs = [0:50:200];


% find ticks
set(h1,'YTick',1+64*(tickVals-realYlims(1))/(realYlims(2)-realYlims(1)))
set(h1,'YtickLabel',tickLabs)