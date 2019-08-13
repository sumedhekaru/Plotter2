function extract_data_from_figure
% This function was writtent to extract HSV frame intensity data and write
% to a txt file.

[fn,pn] = uigetfile('C:\Users\daqop\Desktop\CGLSS-HSV comparisons\*.fig','Open figure file');

if isequal(fn,0)
   disp('User selected Cancel')
   return
end

fg = open([pn fn]);

h = findobj(gca,'Type','line');
x=get(h,'Xdata');
y=get(h,'Ydata');
tit = get(get(gca,'title'),'string');
delete(fg);

x = x{2};
y = y{2};

savefn = tit(1:end-4);
%savefn = [savefn 'txt'];
%direc  = 'C:\data\Video Data\HSV data 2011_KSC\intensity_data\';
fID  = fopen([pn '\' fn '.txt'],'a+');

for i=1:length(x)
    fprintf(fID,'%i\t%0.7f\n',x(i),y(i));
end

fclose(fID);


