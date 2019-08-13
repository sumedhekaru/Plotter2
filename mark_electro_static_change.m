function mark_electro_static_change

hold all
xl = xlim;
[x1, y1] = ginput(1);
plot(xl,[y1 y1],'--k')

[x2, y2] = ginput(1);
plot(xl,[y2 y2],'--k')

text(xl(1)+range(xl)*0.55,(y1+y2)/2,sprintf('%0.1f V/m',y2-y1))
clipboard('copy',y2-y1)



