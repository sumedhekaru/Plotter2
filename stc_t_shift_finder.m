function stc_t_shift_finder
% This is an attempt and I have no idea whether this is going to
% Work. To do this find a return strock and make all time arrangments
% in plotter 2. Then run this.

% Load plotter2 data
try
    h=guidata(findall(0,'Tag','plotter2'));
    handles.g = h.g;
catch
    fprintf('\nE-field plotting parameters are \ncoming from Plotter2. \n\nPlease run plotter2 first.\n')
end



for i = -285:-275
    settings=open('sensor_setting.mat');
    settings.t_shift(22:24) = [i i i];
    save('sensor_setting.mat','-Struct','settings')
    plot_all3(handles.g)
    set(gcf,'position',[30 366 570 450])
end