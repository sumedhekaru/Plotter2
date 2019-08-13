function get_real_proparties
% Get IBP proparties from real data. You should click the data range in the
% current figure (means plot data first).
% Procedure:
%       (1) Open the figure you want to get pusle parameters
%       (2) Run this program
%       (3) Click on the pule start and end (the data range)
%       (4) It will get max and min points automatically
%       (5) Two more clicks to get start time and end time
%       (6) The mat data will be seved in the base folder below.


% Save data in a base folder
% flash type (CG or IC)
%bf = 'C:\Users\Sumedhe\Desktop\IBP_modeling_2013\2748\';
%bf = 'C:\Users\Sumedhe\Desktop\NBP-timing-JGR-2015\Modeling\20150819\';
bf = '/home/daqop/Desktop/WinDesktop/NBP-timing-JGR-2015/Modeling/20150819/';

dateV = '2015-08-19';


% Get plotter2 data




%% Obtaining data from current figure
if ~exist(bf,'dir')
    mkdir(bf)
end

fg = gcf;

h = findobj(fg,'Type','line');
xdata=get(h,'Xdata');
ydata=get(h,'Ydata');

try
    xdata = xdata{end};
    ydata = ydata{end};
catch
    % do nothing
end

% get legend info

lg = get(legend(gca),'String');
lg = lg{1};
lg = lg(1:3);

figure
tools2fig
plot(xdata,ydata);
legend(lg)


%% Obtain IBP data only

% [x1,~] = ginput(1);
% [x2,~] = ginput(1);
% 
% lol = nnz(xdata < x1) + 1;
% ul  = nnz(xdata < x2);

lol = 1;
ul = length(xdata);

xdata = xdata(lol:ul);
ydata = ydata(lol:ul);

cla
plot(xdata,ydata); hold all;
legend(lg)

% Get positive and negative peaks
[ymax, ind1] = max(ydata);
[ymin, ind2] = min(ydata);

tmax = xdata(ind1);
tmin = xdata(ind2);

plot(tmax,ymax,'ro')
plot(tmin,ymin,'ro')


% Get start time and end time
[st_t,~] = ginput(1);
[temp, st_ind] = min(abs(xdata -st_t));
plot(xdata(st_ind),ydata(st_ind),'ro')
[en_t,~] = ginput(1);
[temp, end_ind] = min(abs(xdata -en_t));
plot(xdata(end_ind),ydata(end_ind),'ro')


% offset (static data)
offset = ydata(end_ind) - ydata(st_ind);


data.t = xdata;
data.y = ydata - ydata(st_ind);
data.offset = offset;
data.max = ymax;
data.min = ymin; 
data.max_t = tmax - st_t;
data.min_t = tmin - st_t;
data.t_dur = en_t - st_t;

fprintf('Data for sensor %s\n',lg)
fprintf('\toffset = %0.2f\n',data.offset)
fprintf('\tmax    = %0.2f\n',data.max)
fprintf('\tmin    = %0.2f\n',data.min)
fprintf('\tmax_t  = %0.2f\n',data.max_t*1e6)
fprintf('\tmin_t  = %0.2f\n',data.min_t*1e6)
fprintf('\tt_dur  = %0.2f\n',data.t_dur*1e6)

% save data in the base folder
xx = [bf lg '-' dateV '.mat']
save([bf lg '-' dateV '.mat'],'-Struct','data')







