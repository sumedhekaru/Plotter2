function get_ldar_xyz_from_figure

%% User input

% 1 Sellect the date and the time range you want.
% 2 Create the LDAR2 plot type 12. (Do not plot CGLSS or background map)
% Choose a file name (below) to save data

% File name to save new data
fn = 'C:\Users\sumedhe\Desktop\test.txt';



%% get plotter2 data

% Ask user reall want to replace

if exist(fn,'file')
    answer = questdlg(['The file' fn ' already exist. Do you really want to replace? If the answer is "no", change the file name and rerun the program.'], ...
                      'Warning by Sumedhe',...
                      'yes','no','yes');
else
    answer = 'yes';
end

if ~strcmp(answer,'yes')
    return
end
  

fID = fopen(fn,'w');


handles=guidata(findall(0,'Tag','plotter2'));
g =  handles.g;
sen_set = handles.sen_set;

x0 = 0;
y0 = 0;
z0 = 0;


% Load LDAR2 data in this t1 t2

% Plot LDAR2
if g.mm < 30
    ext=0;
else
    ext=30;
end


dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
    sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);

[CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(sen_set.ldar_r),...
    x0,y0,z0,0);
    
    
% Get x and y from the plot
h = gcf;
axesObjs = get(h, 'Children');
dataObjs = get(axesObjs, 'Children');


h = findobj(gca,'Type','line');
xdata=get(h,'Xdata');
ydata=get(h,'Ydata');

try
    xdata = xdata{end};
    ydata = ydata{end};
catch
    % do nothing
end

figure
plot(xdata,ydata,'ko','markerfacecolor','k','markersize',2)
daspect([1 1 1])

for i = 1:length(xdata) 
    
    ind = find(DLS(:,6)/1000== xdata(i));
    
    if ~isempty(ind)
        fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\n',...
            DLS(ind,10),DLS(ind,6),DLS(ind,7),DLS(ind,8))
        fprintf(fID,'%0.7f\t%0.1f\t%0.1f\t%0.1f\n',...
            DLS(ind,10),DLS(ind,6),DLS(ind,7),DLS(ind,8));
    end
        
end

fclose(fID);


