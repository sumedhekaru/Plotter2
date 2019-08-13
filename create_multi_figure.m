function create_multi_figure

h1 = openfig('E:\Sumedhe\Documents\PBFA_paper\Locating IBP2\Original matlab figures\fig_very_begining1.fig','reuse'); % open figure
ax1 = gca; % get handle to axes of figure
ax1=findall(gcf,'Type','axes');
h2 = openfig('E:\Sumedhe\Documents\PBFA_paper\Locating IBP2\Original matlab figures\fig_very_begining1.fig','reuse');
ax2 = gca;
ax2 = findall(gcf,'Type','axes');
% test1.fig and test2.fig are the names of the figure files which you would % like to copy into multiple subplots

h3 = figure; %create new figure
s1 = subplot(2,1,1); %create and get handle to the subplot axes
hold all
s2 = subplot(2,1,2);
hold all

for i = 1:length(ax1)
    fig1 = get(ax1(i),'children'); %get handle to all the children in the figur
    copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
end
return
    
fig1 = get(ax1,'children'); %get handle to all the children in the figure
fig2 = get(ax2,'children');

copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
copyobj(fig2,s2);