function xml2xy
% 
% % Creating south dakota map
% % (1) Downloaded the shp file from http://www2.census.gov/cgi-bin/shapefiles/county-files?county=46103
% % (2) Opened is using Quantum GIS (added as a vector layer)
% % (3) Right clicked on the file and the save as 'kml'.
% % (4) Ran this program to create mat files with lat lon data
% 
% % fn = 'C:\Users\Sumedhe\Desktop\sdshape\test2.txt';
% % 
% % fid = fopen(fn);
% % data = textscan(fid,'%f, %f','HeaderLines',0);
% % fclose(fid)
% % 
% % 
% % figure
% % plot(data{1},data{2},'ro-')
% clc
% 
% %fid = fopen('C:\Users\Sumedhe\Desktop\sdshape\test.kml');
% %fid = fopen('C:\Users\Sumedhe\Desktop\sdshape\votingDistrict.kml');
% % fid = fopen('C:\Users\Sumedhe\Desktop\sdshape\cousub.kml');
% fid = fopen('C:\Users\Sumedhe\Downloads\test_ms_map.kml');
% 
% figure
% hold all
% 
% tline = fgets(fid);
% 
% d.x = [];
% d.y = [];
% while ischar(tline)
% 
%     tline = fgets(fid);
%     
%     f1 = strfind(tline,'<coordinates>');
%     f2 = strfind(tline,'</coordinates>');
%     
%     if ~isempty(f2)
%         
%         for i = 1:length(f1)
%             
%             data = tline(f1(i)+13:f2(i)-1);
%             [xdata ydata] = get_data(data);
%             
%             if ~isempty(xdata)
%                 d.x = [d.x xdata NaN];
%                 d.y = [d.y ydata NaN];
%             end
%                 
%         end
%     end
%     
% end
% 
% plot(d.y,d.x)
% d.x
% d.y
% 
% 
% save('C:\Users\Sumedhe\Desktop\Plotter2\ms_map.mat','-Struct','d')
% 
% fclose(fid);


% %% Make it x,y
% d = open('C:\Users\Sumedhe\Desktop\Plotter2\ms_map.mat');
% figure
% %[x,y] = latlon2xy(d.y,-d.x,0,34.363882, 89.535407);
% 
% daspect([1 1 1])
% plot(d.x,d.y)
% % save('C:\Users\Sumedhe\Desktop\Plotter2\ms_map.mat','-Struct','d')






function [xdata,ydata] = get_data(data)


    x = strfind(data,' ');
    L = length(x);

    xdata = nan(1,L);
    ydata = nan(1,L);


    y = str2num(data(1:x(1)-1));
    xdata(1) = y(1);
    ydata(1) = y(2);

    for i=1:L-1
        temp = data(x(i)+1:x(i+1)-1);
        y = str2num(temp);
        xdata(i+1) = y(1);
        ydata(i+1) = y(2);
    end

