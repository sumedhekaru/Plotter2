function camera_calibration

% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\photos_2015_06_04\IMG_3849.jpg';
% L1 = 1390;
% L2 = 1173;
% 
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\photos_2015_06_04\IMG_3850.jpg';
% L1 = 2009;
% L2 = 1158;
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\photos_2015_06_04\IMG_3851.jpg';
% L1 = 2309;
% L2 = 782;
% 
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3852.jpg';
% L1 = 1754;
% L2 = 1113;
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3853.jpg';
% L1 = 1615;
% L2 = 1002;
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3854.jpg';
% L1 = 1858;
% L2 = 758;
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3855.jpg';
% L1 = 1553;
% L2 = 979;
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3856.jpg';
% L1 = 1363;
% L2 = 439;
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3857.jpg';
% L1 = 1249;
% L2 = 853;
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3858.jpg';
% L1 = 1276;
% L2 = 880;
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3859.jpg';
% L1 = 936;
% L2 = 398;
% 
% % For Sumedhe's height
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3852.jpg';
% L1 = 2076;
% L2 = 994;
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3853.jpg';
% L1 = 1924;
% L2 = 890;
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3854.jpg';
% L1 = 2403;
% L2 = 582;
% 
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3855.jpg';
% L1 = 2254;
% L2 = 893;
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3857.jpg';
% L1 = 2411;
% L2 = 797;
% 
% a.pfn = 'C:\Users\Sumedhe\Desktop\Cloud_photo\2015_06_12\IMG_3858.jpg';
% L1 = 2438;
% L2 = 829;
% 
% 
% exif = exifread(a.pfn);
% 
% 
% % Pole info
% d = 64*12+167;  % Distance in inches
% h = 57;         % Height in inches
% 
% % Sumedhe info
% d = 54*12;
% d = 38*12;
% d = 22*12;
% h = 66.25;
% 
% 
% I = imread(a.pfn);
% % Camera info
% a.p = 1.8e-6; %Fixel Size
% 
% 
% % Estimaate focul length
% f = d * (L1 - L2) * a.p / h;
% 
% % Plot the image
% figure; 
% imagesc(I);
% hold all;
% plot([1 3264],[L1 L1],'r')
% plot([1 3264],[L2 L2],'r')
% str = sprintf('Calculated focal length = %0.3f mm',f*1000);
% text(100, L2-100,str,'color','w')
% str = sprintf('EXIF focal length = %0.3f mm',exif.FocalLength);
% text(100, L2-200,str,'color','w')
% 
% 
% 
% 
% daspect([1 1 1])
% hold all
% 
% tools2fig
% 
data = [0 0 % Origin?
           6 6.407 % pole 
           23.8 25.127 %pole
           42.8 45.087 %pole
           18.2 18.926 % Pole
           18.2 19.050 % sumedhe
           17.3 18.100 % pole
           17.3 18.205 % sumedhe
           31.2 32.479 % pole
           31.2 32.061 % suemedhe
           16.4 16.648 % pole
           16.4 16.862 % sumedhe
           26.0 27.282 % pole
           11.3 11.692 % pole
           11.3 11.577 % Sumedhe
           11.3 11.541 % Sumedhe
           15.2 15.885 % pole           
           ];

f_photo = data(:,1);
f_cal = data(:,2);

figure
plot(f_cal,f_photo,'ro')
title('Canon PowerShot S5 IS focal length calibration')
xlabel('EXIF - focal length (mm)');
ylabel('Calculated focal length (mm)')