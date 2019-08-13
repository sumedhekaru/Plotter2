function output_txt = curser(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');

sec=pos(1);
   hour   = fix(sec/3600);      % get number of hours
   sec    = sec - 3600*hour;    % remove the hours
   minute = fix(sec/60);        % get number of minutes
   sec    = sec - 60*minute;    % remove the minutes
   second = sec;

t=[num2str(hour),':',num2str(minute),':',num2str(second,'%.6f')];

output_txt = {['X: ',t],...
    ['Y: ',num2str(pos(2))]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),'%.6f')];
end
