function get_line_data_from_plot_for_IBP_modeling

% 1. First plot the data you want
% 2. Make sure you have sat offset to zero
% 2. Change user inputs below aproprately and then run

%% User inputs
bf = 'C:\Users\Sumedhe\Desktop\NBP_modeling_2016\201110814_70844\';
sns = 'OVD';

%% Start the program

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

figure
ff1 = plot(xdata,ydata);

%ynew = moving(ydata,1000);
%hold all
%plot(xdata,ynew);

data.t = xdata;
data.y = ydata;


fname = [bf sns '.mat'];

if exist(fname,'file')
    % Construct a questdlg with three options
    choice = questdlg('File exist. Would you like to replace?', ...
        'File exist!','Yes','No','Cancel','No');
    % Handle response
    if strcmp(choice,'Yes')
        save(fname,'-struct','data')  
    end
else
    save(fname,'-struct','data')

end

disp(fname)
