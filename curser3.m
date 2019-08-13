function output_txt = curser3(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');

[x y z] = scr2ldar(pos(1),pos(2));

t=sprintf('%.8f',pos(1));

output_txt = {sprintf('X : %i Y : %i',pos(1),pos(2)),...
    sprintf('x : %.1fm',x),...
    sprintf('y : %.1fm',y),...
    sprintf('z : %.1fm',z)};


clipboard('copy', [x y z])

function [x y z] = scr2ldar(N_x,N_y)

h=guidata(findall(0,'Tag','hsv_analyzer'));

b=h.b;

f = 8.0e-3;  %Focal length
p = 20.0e-6; %Fixel Size


r = sqrt((((b.ref_x-b.cam_x)).^2+(b.ref_y-b.cam_y ).^2 ).*(1+(p^2.*(b.scr_x-N_x ).^2)./f^2 ));

theta = atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x ))+atan(p*(b.scr_x-N_x)/f);

x = b.cam_x + r.*cos(theta);
y = b.cam_y + r.*sin(theta);
z = b.cam_z + r.*((b.scr_y-N_y).*p)./f;

