function get_line_data_from_plots

% File name to save new data
%fn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F02\K02_10KHz_leader_data.mat';
% fn = 'C:\Users\Sumedhe\Desktop\Ip_cal_data.mat'
% 
% if exist(fn,'file')
%     answer = questdlg(['The file' fn ' already exist. Do you really want to replace? If the answer is "no", change the file name and rerun the program.'], ...
%                       'Warning by Sumedhe',...
%                       'yes','no','yes');
% else
%     answer = 'yes';
% end
% 
% 
% 
% 
% if ~strcmp(answer,'yes')
%     return
% end

[fn,pn] = uigetfile('C:\Users\daqop\Desktop\CGLSS-HSV comparisons\*.fig','Open figure file');

if isequal(fn,0)
   disp('User selected Cancel')
   return
end

fg = open([pn fn]);

h = findobj(fg,'Type','line');
xdata=get(h,'Xdata');
ydata=get(h,'Ydata');

try
    xdata = xdata{end};
    ydata = ydata{end};
catch
    % do nothing
end

figure
plot(xdata,ydata)

%ynew = moving(ydata,1000);
%hold all
%plot(xdata,ynew);

data.t = xdata;
data.v = ydata

save([pn fn(1:end-3) 'mat'],'-struct','data')