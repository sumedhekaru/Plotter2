function test_line_angled

[x1,y1] = ginput(1);
[x2,y2] = ginput(1);
hold all;
plot([x1 x2],[y1 y2],'k')