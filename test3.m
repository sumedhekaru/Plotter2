function test3

figure
plot(rand(1,100))
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'units','centimeters')
A = get(gcf,'outerposition');
set(gcf,'paperunits','centimeters','paperposition',A)
print('-dpng', 'C:\Users\sumedhe\Desktop\text.png', '-r150');
delete(gcf)