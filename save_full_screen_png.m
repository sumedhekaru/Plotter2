function save_full_screen_png(fn,dpi)

if nargin < 2
    dpi = 150;
end

set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'units','centimeters')
A = get(gcf,'outerposition');
set(gcf,'paperunits','centimeters','paperposition',A)
print('-dpng', fn, sprintf('-r%i',dpi));
