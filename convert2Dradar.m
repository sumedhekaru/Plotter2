function r = convert2Dradar(x1,x2,y1,y2,x,y)


% slope and intercept of the line (need later)
m = (y2 - y1)/(x2 - x1);
b = (x2*y1 - x1*y2)/(x2 - x1);

% Snap the location coordinates to line
x = (m*y+x-m*b)/(m*m + 1);
y = (m*m*y+m*x+b)/(m*m + 1);

% Make them 2D
r = sqrt((x - x1).^2 + (y - y1).^2);
