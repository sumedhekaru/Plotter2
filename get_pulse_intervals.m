function get_pulse_intervals
% This function was written to obtain pulse group intervals of IC flashes.
%
%   Requirements: The IC flash should be plotted and it should be
%                 the current figure just before this run.
%
%   Usage : get_pulse_intervals
%           Once run, click on the plot using left mouse clicks (order or
%           number of clicks doesn't matter. Once finish right click
%           anywhere on the plot to exit.
%
%           The pulse time and amplitude of the pulses will be written on a
%           file in Plotter2 folder with name
%           'PulseIntervalData-yyyymmdd.csv' format. If the file already
%           exsist, it will update the that file.
%
%   History:
%       2014-02-21 Created by Sumedhe Karunarathne
%       2014-02-25 Modified to get peak to peak info



% get plotter2 data
h=guidata(findall(0,'Tag','plotter2'));
g = h.g;

% Open file to write data
fn = sprintf('PulseIntervalData-%s%s%s.txt',g.YYYY{:},g.MM{:},g.DD{:});
fID = fopen(fn,'a+');

% Variables
t = [];
h = [];

% Hold the plot
hold all

% Get lengend
lg = get(legend(gca),'String');
sn = lg{1}(1:3);

% Get user input
[t1,E1,button] = ginput(1);

% Do untill right click
while button == 1
    
    % Give visual feed back of the point chosed
    h1 = plot(t1,E1,'sr');
    
    % Get peak to peak info
    Epp = get_peak2peak(t1);
    
    % save data
    h = [h h1];    
    t = [t, t1];
    
    % Print info
    fprintf('%s\t%13.6f\t%5.1f\n',sn,t1,Epp)
    fprintf(fID,'%s\t%13.6f\t%5.1f\n',sn,t1,Epp);
    
    % Get user input
    [t1,E1,button] = ginput(1);    
end

% Delete visual feedbacks
delete(h)
fclose(fID);



function Epp = get_peak2peak(t1)

dobj = get(get(gcf, 'Children'), 'Children');

for i = 1:length(dobj)
    
    type = get(dobj{i},'Type');
    
    if strcmp(type,'line')
        t = get(dobj{end}, 'XData');
        t = t{end};
        E = get(dobj{end}, 'YData');
        E = E{end};
        
        if  length(t) > 50
            break
        end
    end
end



lol = sum(t < (t1-50e-6))+1;
ul  = sum(t < (t1+50e-6));

% figure
% plot(t(lol:ul),E(lol:ul))
Epp = range(E(lol:ul));
    