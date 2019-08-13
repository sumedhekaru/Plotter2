function plot_xy
% This function is writtent to plot xy data (in km) on a current plot. 
% Note: 
%      Date, time, and witch plots to plot will be loaded from Plotter2
%      Click on the plot (the one you wanted to overlay)to make it active
%      plot
%      Legend is not controlled, either you have to do it manually, or
%      start with a graph with apropriate legend.
%   
% Modification history
%      2015-03-23 Created by Sumedhe Karunarathne

%% User inputs

cl = 'w'; % color value: letter like "r" or 3 value RGB vector like [1 0 1]
%all the data irrespective of type will be plotted in one color
% Marker size
mz = 4;

% Plot on current figure (if not a new figure will be opened)
pcf = 1;

%% Start the program

% Get current plotter2 data and settings
try
    h=guidata(findall(0,'Tag','plotter2'));
catch
    disp('Run plotter2 first! Date info coming from plotter2')
    return
end

g = h.g;
sen_set = h.sen_set;


% Turn off all LP/CH plots
g.lpgraphs = g.lpgraphs - g.lpgraphs;
g.chgraphs = g.chgraphs - g.chgraphs;



a = load_data(g,sen_set);

if pcf
    figure(gcf)
else
    figure
    xlabel('East (km)')
    ylabel('North (km)')
    title(sprintf('xy data\n%s-%s-%s  UT:%s - %s',g.YYYY{:},g.MM{:},g.DD{:}, sec2hhmmss(g.t1),sec2hhmmss(g.t2)))
    box on
    daspect([1 1 1])
end


hold all

% Plot ldar
if ~isnan(a.DLS(1,1))
    plot(a.DLS(:,6)/1000,a.DLS(:,7)/1000,'o','markersize',mz,'color',cl,'markerfacecolor',cl)
end

% Plot CG
if ~isnan(a.CG(1,1))
    plot(a.CG(:,6)/1000,a.CG(:,7)/1000,'s','markersize',mz,'color',cl,'markerfacecolor',cl)
end


% Plot PBFA
if ~isnan(a.PBFA(1,1))
    plot(a.PBFA(:,3)/1000,a.PBFA(:,4)/1000,'d','markersize',mz,'color',cl,'markerfacecolor',cl)
end

% Plot PBFA-auto
if ~isnan(a.PBFA2(1,1))
    plot(a.PBFA2(:,3)/1000,a.PBFA2(:,4)/1000,'v','markersize',mz,'color',cl,'markerfacecolor',cl)
end

% Plot PBFA-old
if ~isnan(a.PBFA_old(1,1))
    plot(a.PBFA_old(:,3)/1000,a.PBFA_old(:,4)/1000,'^','markersize',mz,'color',cl,'markerfacecolor',cl)
end


if ~pcf
    florida_map
end
