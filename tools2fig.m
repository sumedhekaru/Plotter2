function tools2fig()

f = uimenu('Label','Plotter');
uimenu(f,'Label','Update Y Grids','Callback','plotter_tools(1)','Accelerator','Y');
uimenu(f,'Label','Update X Grids','Callback','plotter_tools(2)','Accelerator','X');
uimenu(f,'Label','Update XY Grids','Callback','plotter_tools(29)');
uimenu(f,'Label','Update Time Range','Callback','plotter_tools(3)','Accelerator','T');
uimenu(f,'Label','Save Figure As...','callback','plotter_tools(16)');
uimenu(f,'Label','Recall figure size','callback','set_figure_size');
uimenu(f,'Label','Find Delta t','Callback','plotter_tools(4)','Accelerator','D');
uimenu(f,'Label','Find Delta y','Callback','plotter_tools(10)','Accelerator','F');
uimenu(f,'Label','hh:mm:ss Data Curser','Callback','plotter_tools(5)','Accelerator','K');
uimenu(f,'Label','ss.ssssss Data Curser','Callback','plotter_tools(6)','Accelerator','M');
uimenu(f,'Label','Location Data Curser','Callback','plotter_tools(26)');
uimenu(f,'Label','Delete PBFA point','Callback','plotter_tools(30)');
%uimenu(f,'Label','CGLSS Data Curser','Callback','plotter_tools(17)');
%uimenu(f,'Label','PBFA Data Curser','Callback','plotter_tools(18)');
%uimenu(f,'Label','PBFA-A Data Curser','Callback','plotter_tools(20)');
%uimenu(f,'Label','PBFA-O Data Curser','Callback','plotter_tools(21)');
%uimenu(f,'Label','LDAR2 Data Curser','Callback','plotter_tools(19)');
%uimenu(f,'Label','NLDN2 Data Curser','Callback','plotter_tools(22)');
%uimenu(f,'Label','LINET Data Curser','Callback','plotter_tools(23)');
uimenu(f,'Label','X Zoom!','Callback','plotter_tools(7)','Accelerator','G');
uimenu(f,'Label','Y Zoom!','Callback','plotter_tools(8)','Accelerator','H');
uimenu(f,'Label','XY Zoom!','Callback','plotter_tools(9)','Accelerator','J');
uimenu(f,'Label','Calculate PBFA points','Callback','peak_modifier','Accelerator','L');
uimenu(f,'Label','Pulse Modeling','Callback','pulse4');
uimenu(f,'Label','Vedio Framing','callback','plotter_tools(11)');
uimenu(f,'Label','Fix for Publishing','callback','fix_fig');
uimenu(f,'Label','Fix the title','callback','plotter_tools(15)');
uimenu(f,'Label','Find Rise & Fall Times','callback','rise_fall_times','Accelerator','R');
uimenu(f,'Label','PBFA Auto Calculations','callback','plotter_tools(24)','Accelerator','B');
uimenu(f,'Label','Nearest Pulses','callback','plotter_tools(25)','Accelerator',',');
uimenu(f,'Label','Get Pulse Intervals','callback','get_pulse_intervals','Accelerator','I');
uimenu(f,'Label','Get Pulse Intervals2','callback','get_pulse_intervals2','Accelerator','9');
uimenu(f,'Label','Show Pulse Intervel data','callback','show_pulse_intervals','Accelerator','U');
uimenu(f,'Label','Name figure a, b, c','callback','plotter_tools(28)');
uimenu(f,'Label','Add annotations','callback','anotation(gcf)');

