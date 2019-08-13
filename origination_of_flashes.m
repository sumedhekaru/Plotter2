function origination_of_flashes
% origination of flashes pjoject : plot flash origination data in x,y and z
%fn = '/home/daqop/Desktop/WinDesktop/Nadee/Nadee-official-sumpc/Origination-of-flashes-in-samespot/Flash-origination-07222011_V04.xlsx';
% fn = 'C:/Users/Sumedhe/Desktop/Nadee/Nadee-official-sumpc/Origination-of-flashes-in-samespot/Flash-origination-07222011_V04.xlsx';
%sheet = 1; % for over the sea
%sheet =2 ; % over the land all together
% sheet = 3; % over the land t-storm 1
% change this as per your data
%xlRange = 'A1:J9'; %for over the sea
%xlRange = 'A1:J47';  %over the land both storms together
%xlRange = 'A1:J23'; %over the land storm 1 seperately

 

% Echange sensor sites 
sns_x = [ -0.7524, 0.3372, -0.3149, -1.1658, -0.2701, -2.7640, -6.0091, 0.1825, -5.7394, -2.0637]*10;
sns_y = [  1.6555, 0.4446, -0.6838, -1.7020,  0.2631,  4.9254, -3.3983, -5.3008, 1.1923, 0.1569]*10;
sns_n = { 'K02', 'K14',    'K24',    'BCC',    'K17',    'EDW',    'STC',    'FLT',    'OVD',    'FFI'};

 %% Plot over the sea storm
% fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011_V04.xlsx';
% 
%   sheet = 1; % for over the sea
%   xlRange = 'A1:J7'; %for over the sea
%   [ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
% %  % Graph option 1 : x (east vs y (north)
%   figure
% %  %subplot(no of rows (eg 3), no of columns(eg 2), desired position of figure (1:2)),% 1:2 means use space 1 and 2 together
% %  
%  subplot(3,2,1:2)
%  hold all;
%  % use color map for data
%  %cm contains color values for RGB in 3 columns 
%  % cm has 64 colors , so L is 64
%  cm = jet;
%  L = length(cm); % 64
%  % data will be colored according to time pbfat
%  %pbfat is a matrix 
%  pbfat = ndata(:,6);
%  % cl matrix will hold all the color data for pbfat (following equation for cl is derived using gradient
%  % method to calculate unknown color value for a known time pbfat)
%  cl = 1+ round((L-1) *(pbfat - pbfat(1))/(pbfat(end) - pbfat(1)));
% 
% 
% for i= 1:7
%  
%  pbfax = ndata(i,7)/1000;
%  pbfay = ndata(i,8)/1000;
%  pbfaz = ndata(i,9)/1000;
% 
%   plot(pbfax,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:))
%  text(pbfax+0.25,pbfay,num2str(i))
% end
% plot([20 35 35 20 20], [30 30 20 20 30],'k--');
% 
% 
% % plot florida map and sensor sites
%  plot(sns_x,sns_y,'kp','markerfacecolor','k')
%  text(sns_x,sns_y,sns_n)
%  % Florida map
%  florida_map
%  xlim([-40, 40])
%  ylim([-15, 40])
%  xlabel('East(km)')
%  ylabel('North(km)')
%  box on
%  daspect([1,1,1])
%  set(gca,'position',[0.1122    0.7183    0.7417    0.2624])
%   subplot(3,2,2)
% cg = find(strcmp(txtdata(:,3),'CG'));
% ic = find(strcmp(txtdata(:,3),'IC'));
% 
%  for i= 1:9
%  pbfat = ndata(i,6);
%  pbfaz = ndata(i,9)/1000;
% 
%  
%  plot(ndata(cg(i),6),ndata(cg(i),9),'rx','markerfacecolor','r')
%  hold all
%  plot(ndata(ic(i),6),ndata(ic(i),9),'bx','markerfacecolor','b')
%  
% end 
%  % Graph option 2: x (east) vs Z (altitude)
% 
%  subplot(3,2,3)
%  hold all;
%  for i= 1:7
%  pbfat = ndata(i,6);
%  pbfax = ndata(i,7)/1000;
%  pbfay = ndata(i,8)/1000;
%  pbfaz = ndata(i,9)/1000;
%  %indx = ndata(i,4)
%  % color ceded data will be plotted cm(raw, column), seperate columns for
%  % RGB , all we want to give is row number cl(i)
%  plot(pbfax,pbfaz,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:))
%  % numbering the data point
%  text(pbfax+0.35,pbfaz,num2str(i))
%   
%  end
%  xlim([23, 33])
%  ylim([4, 10])
%  xlabel('East(km)')
%  ylabel('Altitude(km)')
%  box on
%  daspect([1,1,1])
% set(gca, 'position',[0.1045    0.3007    0.4299    0.4332]);
%   % Graph option 3: Histogram
%  
%  subplot(3,2,4)
% %  hist(ndata(:,9)/1000,0.5:1:14.5)
% %  xlim([0,15])
% %  ylim([0,5])
% %  set(gca,'xgrid','on','ygrid','off')
% %  set(gca,'xtick',[0:1:15])
% %  xlabel('altitude(km)')
% %  ylabel('Number of flashes')
% 
% %disp(txtdata(:,3))
% cg = find(strcmp(txtdata(:,3),'CG'));
% ic = find(strcmp(txtdata(:,3),'IC'));
% 
% [frequency, bin] = hist(ndata(cg,9)/1000,4:1:12);
% [frequency2, bin] = hist(ndata(ic,9)/1000,4:1:12);
% 
%  %bar(bin',[frequency;frequency2]',1,'grouped')
%  bar(bin',[frequency;frequency2]')
%  xlim([4,10])
%  ylim([0,6])
%  set(gca,'xgrid','on','ygrid','off')
%  set(gca,'xtick',[4:1:10])
%  set(gca,'xticklabel',{'3',' ','5',' ','7',' ','9',' ','11',' '})
%   xlabel('Altitude(km)')
%   ylabel('Number of flashes')
%   legend ('CG','IC')
%  set(gca, 'position',[0.6069    0.4109    0.2352    0.2472]);
% %   
%  
% % Graph option 4 : x (east vs y (north)
% % setup
%  
%  subplot(3,2,5)
%  hold all;
%  cm = jet;
%  L = length(cm);
%  pbfat = ndata(:,6);
%  cl = 1+ round(63 *(pbfat - pbfat(1))/(pbfat(end) - pbfat(1)));
% 
% for i= 1:7
%  
%  pbfax = ndata(i,7)/1000;
%  pbfay = ndata(i,8)/1000;
%  pbfaz = ndata(i,9)/1000;
% 
%   plot(pbfax,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:))
%  text(pbfax+0.35,pbfay,num2str(i))
% end
% 
% % plot florida map and sensor sites
%  plot(sns_x,sns_y,'kp','markerfacecolor','k')
%  text(sns_x,sns_y,sns_n)
%  % Florida map
%  florida_map
%   xlim([23, 33])
%  ylim([20, 28])
%  xlabel('East(km)')
%  ylabel('North(km)')
%  box on
%  daspect([1,1,1])
%  set(gca, 'position',[0.0944    0.0401    0.4413    0.3224]);
% 
% % Graph option 5: z (altitude) vs y (north)
% % 
%  subplot(3,2,6)
%  hold all;
%  for i= 1:7
%  pbfat = ndata(i,6);
%  pbfax = ndata(i,7)/1000;
%  pbfay = ndata(i,8)/1000;
%  pbfaz = ndata(i,9)/1000;
%  indx = ndata(i,4)
%  
%  plot(pbfaz,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:))
%  text(pbfaz+0.25,pbfay,num2str(i))
%  
%  
%  end
%  xlim([4, 10])
%  ylim([20, 28])
%  xlabel('Altitude(km)')
%  ylabel('North(km)')
%  box on
%  daspect([1,1,1])
%  set(gca, 'position',[0.4596    0.0457    0.5297    0.3151]);

 %% plot over the land flashes (two thumnderstorms together)
% fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011_V04.xlsx';
%   sheet =2 ; % over the land all together
%  xlRange = 'A1:J47';  %over the land both storms together
%  [ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
% %  
% %  %  % Graph option 1 : x (east vs y (north)
% % 
%    figure
% % %  %subplot(no of rows (eg 3), no of columns(eg 2), desired position of figure (1:2)),% 1:2 means use space 1 and 2 together
% % %  
% %   subplot(3,2,1:2)
%    hold all;
% % %  % use color map for data
% % %  %cm contains color values for RGB in 3 columns 
% % %  % cm has 64 colors , so L is 64
%    cm = jet;
%    L = length(cm); % 64
% % %  % data will be colored according to time pbfat
% % %  %pbfat is a matrix 
%    pbfat = ndata(:,6);
% % %  % cl matrix will hold all the color number data for pbfat (following equation for cl is derived using gradient
% % %  % method to calculate unknown color value for a known time pbfat)
%    cl = 1+ round((L-1) *(pbfat - pbfat(1))/(pbfat(end) - pbfat(1)));
% % % 
%   for i= 1:47
%  
%    pbfax = ndata(i,7)/1000;
%    pbfay = ndata(i,8)/1000;
%    pbfaz = ndata(i,9)/1000;
%    plot(pbfax,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
%    %text(pbfax+0.25,pbfay,num2str(i))
%   end
%   plot([-14 -20 -20 -14 -14], [-12 -12 -2 -2 -12],'k--');% tstorm- 2
%    %plot([-4 -10 -10 -4 -4], [2 2 12 12 2],'k--'); %tstorm - 1
% % % 
% % % 
% % % % plot florida map and sensor sites
%    plot(sns_x,sns_y,'kp','markerfacecolor','k')
%    text(sns_x,sns_y,sns_n)
% % %  % Florida map
%    florida_map
%    xlim([-30, 10])
%    ylim([-20, 20])
%    xlabel('East(km)')
%    ylabel('North(km)')
%    box on
%    daspect([1,1,1])
%   set(gca,'position',[0.0360    0.7706    0.7750    0.2157])
% % 
% %  % Graph option 2: x (east) vs Z (altitude)
% % 
%   subplot(3,2,3)
%   hold all;
%   for i= 1:47
%   pbfat = ndata(i,6);
%   pbfax = ndata(i,7)/1000;
%   pbfay = ndata(i,8)/1000;
%   pbfaz = ndata(i,9)/1000;
% %  indx = ndata(i,4)
% %  % color ceded data will be plotted cm(raw, column), seperate columns for
% %  % RGB , all we want to give is row number cl(i)
%   plot(pbfax,pbfaz,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
% %  % numbering the data point
%   %text(pbfax+0.35,pbfaz,num2str(i))
% %   
%  end
%   xlim([-20, -3])
%   ylim([0, 15])
%   xlabel('East(km)')
%   ylabel('Altitude(km)')
%   box on
%   daspect([1,1,1])
%   set(gca,'position',[0.0830    0.4976    0.3347    0.2157])
% % 
% %   % Graph option 3: Histogram
% %  
%   subplot(3,2,4)
%   hist(ndata(:,9)/1000,0.5:1:14.5)
%   xlim([0,15])
%   ylim([0,15])
%   set(gca,'xgrid','on','ygrid','off')
%   set(gca,'xtick',[0:1:15])
%   set(gca,'xticklabel',{'0',' ','2',' ','4',' ','6',' ','8',' ','10',' ','12',' ','14',' '})
%   xlabel('Altitude(km)')
%   ylabel('Number of flashes')
%   set(gca,'position',[0.5374    0.5127    0.2902    0.2153])
% %  
% %  
% % % Graph option 4 : x (east vs y (north)
% % % setup
% %  
%   subplot(3,2,5)
%  hold all;
%   cm = jet;
%   L = length(cm);
%   pbfat = ndata(:,6);
%   cl = 1+ round(63 *(pbfat - pbfat(1))/(pbfat(end) - pbfat(1)));
% % 
%  for i= 1:47
% %  
%   pbfax = ndata(i,7)/1000;
%   pbfay = ndata(i,8)/1000;
%   pbfaz = ndata(i,9)/1000;
% % 
%    plot(pbfax,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
% %  text(pbfax+0.35,pbfay,num2str(i))
%  end
% % 
% % % plot florida map and sensor sites
%   plot(sns_x,sns_y,'kp','markerfacecolor','k')
%   text(sns_x,sns_y,sns_n)
% %  % Florida map
%   florida_map
%    xlim([-20, -3])
%   ylim([-15, 15])
%   xlabel('East(km)')
%   ylabel('North(km)')
%   box on
%   daspect([1,1,1])
%   set(gca,'position',[-0.0376    0.0359    0.5737    0.4271])
% %  
% % 
% % % Graph option 5: z (altitude) vs y (north)
% % % 
%   subplot(3,2,6)
%   hold all;
%   for i= 1:47
%   pbfat = ndata(i,6);
%   pbfax = ndata(i,7)/1000;
%   pbfay = ndata(i,8)/1000;
%   pbfaz = ndata(i,9)/1000;
% %  indx = ndata(i,4)
% %  
%   plot(pbfaz,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
% %  text(pbfaz+0.25,pbfay,num2str(i))
% %  
% %  
%   end
%   xlim([0, 15])
%   ylim([-15, 15])
%   xlabel('Altitude(km)')
%   ylabel('North(km)')
%   box on
%   daspect([1,1,1])
%  set(gca,'position',[0.4120    0.0417    0.5394    0.4276])
 
  %% plot over the land storms seperately
%  % Thunderstorm 1
% %  
% fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011_V04.xlsx';
% 
%   sheet = 3; % over the land t-storm 1
%   xlRange = 'A1:J24'; %over the land storm 1 seperately
%   [ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
% % 
%   figure
% 
%  % Graph option 1 : x (east vs y (north)  %  %subplot(no of rows (eg 3), no of columns(eg 2), desired position of figure (1:2)),% 1:2 means use space 1 and 2 together
%   
%   subplot(3,2,1:2)
%    hold all;
%   % use color map for data
%   %cm contains color values for RGB in 3 columns 
%   % cm has 64 colors , so L is 64
%    cm = jet;
%    L = length(cm); % 64
%   % data will be colored according to time pbfat
%   %pbfat is a matrix 
%    pbfat = ndata(:,6);
%   % cl matrix will hold all the color number data for pbfat (following equation for cl is derived using gradient % method to calculate unknown color value for a known time pbfat)
%    cl = 1+ round((L-1) *(pbfat - pbfat(1))/(pbfat(end) - pbfat(1)));
%  
%  for i= 1:24
% % %  
%   pbfax = ndata(i,7)/1000;
%   pbfay = ndata(i,8)/1000;
%   pbfaz = ndata(i,9)/1000;
% % % 
%  plot(pbfax,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
%  %text(pbfax+0.25,pbfay,num2str(i))
% end
% % plot([20 35 35 20 20], [35 35 20 20 35],'k--');
% 
%  % plot florida map and sensor sites
%  plot(sns_x,sns_y,'kp','markerfacecolor','k')
%  text(sns_x,sns_y,sns_n)
%  % Florida map
%  florida_map
%  xlim([-20, 5])
%  ylim([-5, 20])
%  xlabel('East(km)')
%  ylabel('North(km)')
%  box on
%  daspect([1,1,1])
%  %set(gca,'position',[0.0360    0.7706    0.7750    0.2157])
% % % 
%  % Graph option 2: x (east) vs Z (altitude)
% % % 
%    subplot(3,2,3)
%    hold all;
%    for i= 1:24
%    pbfat = ndata(i,6);
%    pbfax = ndata(i,7)/1000;
%    pbfay = ndata(i,8)/1000;
%    pbfaz = ndata(i,9)/1000;
% % %  indx = ndata(i,4)
% % %  % color ceded data will be plotted cm(raw, column), seperate columns for
% % %  % RGB , all we want to give is row number cl(i)
%   plot(pbfax,pbfaz,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
% % %  % numbering the data point
%   text(pbfax,pbfaz,num2str(i))
%   end
%   xlim([-10, -4])
%   ylim([4, 12])
%   xlabel('East(km)')
%   ylabel('Altitude(km)')
%   box on
%   daspect([1,1,1])
%  
%    %set(gca,'position',[0.0830    0.4976    0.3347    0.2157])
% % % 
% % %   % Graph option 3: Histogram
% % %  
% %  subplot(3,2,4)
% %   hist(ndata(:,9)/1000,0.5:1:14)
% %   xlim([0,14])
% %   ylim([0,15])
% %   set(gca,'xgrid','on','ygrid','off')
% %   set(gca,'xtick',[0:1:14])
% %   set(gca,'xticklabel',{'0',' ','2',' ','4',' ','6',' ','8',' ','10',' ','12',' ','14'})
% %   xlabel('Altitude(km)')
% %   ylabel('Number of flashes')
% %   %set(gca,'position',[0.5374    0.5127    0.2902    0.2153])
% 
%  % disp(txtdata(:,3))
%    subplot(3,2,4)
%    cg = find(strcmp(txtdata(:,3),'CG'));
%    ic = find(strcmp(txtdata(:,3),'IC'));
% %   
%     [frequency, bin] = hist(ndata(cg,9)/1000, 0.5:1:14)
%    [frequency2, bin2] = hist(ndata(ic,9)/1000, 0.5:1:14)
%     size(frequency)
%     size(frequency2)
%     bar(bin',[frequency; frequency2]',1,'grouped')
%     xlim([4,12])
%     ylim([0,8])
%     set(gca,'xgrid','on','ygrid','off')
%     set(gca,'xtick',[3:1:12])
%     %set(gca,'xticklabel',{'3',' ','5',' ','7',' ','9',' ','11',' '})
%     xlabel('Altitude(km)')
%     ylabel('Number of flashes')
%     legend ('CG','IC')
% 
% 
%  % Graph option 4 : x (east vs y (north)
%  % setup
%  
%   subplot(3,2,5)
%   hold all;
%   cm = jet;
%   L = length(cm);
%   pbfat = ndata(:,6);
%    cl = 1+ round(63 *(pbfat - pbfat(1))/(pbfat(end) - pbfat(1)));
% % % 
%  for i= 1:24
% % %  
%    pbfax = ndata(i,7)/1000;
%    pbfay = ndata(i,8)/1000;
%    pbfaz = ndata(i,9)/1000;
% % % 
%    plot(pbfax,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
%    %text(pbfax+0.35,pbfay,num2str(i))
%    text(pbfax,pbfay,num2str(i))
%     end
% % % 
% % % % plot florida map and sensor sites
%    plot(sns_x,sns_y,'kp','markerfacecolor','k')
%    text(sns_x,sns_y,sns_n)
% % %  % Florida map
%    florida_map
%    xlim([-10, -4])
%    ylim([2, 12])
%    xlabel('East(km)')
%    ylabel('North(km)')
%    box on
%    daspect([1,1,1])
%    set(gca,'position',[-0.0376    0.0359    0.5737    0.4271])
% 
% % % % Graph option 5: z (altitude) vs y (north)
% 
%    subplot(3,2,6)
%    hold all;
%    for i= 1:24
%    pbfat = ndata(i,6);
%    pbfax = ndata(i,7)/1000;
%    pbfay = ndata(i,8)/1000;
%    pbfaz = ndata(i,9)/1000;
% % %  indx = ndata(i,4)
% % %  
%    plot(pbfaz,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
%   %text(pbfaz+0.25,pbfay,num2str(i))
%   text(pbfaz,pbfay,num2str(i))
%    end
%    xlim([4, 12])
%    ylim([2,12])
%    xlabel('Altitude(km)')
%    ylabel('North(km)')
%    box on
%    daspect([1,1,1])
%    set(gca,'position',[0.4120    0.0417    0.5394    0.4276])
   %%
   %plot over the land storm 1 in 2 colors only CG- blue, IC red
%    
%    fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011_V04.xlsx';
% 
%   sheet = 3; % over the land t-storm 1
%   xlRange = 'A1:J24'; %over the land storm 1 seperately
%   [ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
% 
  %figure

 % Graph option 1 : x (east vs y (north)  %  %subplot(no of rows (eg 3), no of columns(eg 2), desired position of figure (1:2)),% 1:2 means use space 1 and 2 together
  
  %subplot(3,2,1:2)
  % hold all;
  % use color map for data
  %cm contains color values for RGB in 3 columns 
  % cm has 64 colors , so L is 64
%    cm = jet;
%    L = length(cm); % 64
%   % data will be colored according to time pbfat
%   %pbfat is a matrix 
%    pbfat = ndata(:,6);
%   % cl matrix will hold all the color number data for pbfat (following equation for cl is derived using gradient % method to calculate unknown color value for a known time pbfat)
%    cl = 1+ round((L-1) *(pbfat - pbfat(1))/(pbfat(end) - pbfat(1)));
%    
%    indx = find(strcmp(txtdata(:,3),'IC'));
%    indx2 = find(strcmp(txtdata(:,3),'CG'));
%    
%     % ndata(indx,3)
%     figure
%     subplot(4,3,1)
%     hold all
%     %xxx = ndata(indx2,7)
%     plot(ndata(indx,7)/1000,ndata(indx,8)/1000,'ko','markerfacecolor','r')
%     %text((ndata(indx,7))/1000+.25,ndata(indx,9)/1000,num2str(ndata(indx,3)))
%     plot(ndata(indx2,7)/1000,ndata(indx2,8)/1000,'ko','markerfacecolor','b')
% %     text((ndata(indx2,7))/1000+.25,ndata(indx2,8)/1000,num2str(ndata(indx2,3)))
%     xlim([-25 40])
%     ylim([-5 40])
%     title('Over the land- storm 1');
%  
%  for i= 1:24
% % %  
%   pbfax = ndata(i,7)/1000;
%   pbfay = ndata(i,8)/1000;
%   pbfaz = ndata(i,9)/1000;
% % % 
%  plot(pbfax,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
%  %text(pbfax+0.25,pbfay,num2str(i))
% end
% plot([20 35 35 20 20], [35 35 20 20 35],'k--');

 % plot florida map and sensor sites
%  plot(sns_x,sns_y,'kp','markerfacecolor','k')
%  text(sns_x,sns_y,sns_n)
%  % Florida map
%  florida_map
%  %xlabel('East(km)')
%  ylabel('North(km)')
%  box on
%  daspect([1,1,1])
%  
%  
%  subplot(4,3,4)
%  hold all
%  %xxx = ndata(indx2,7)
%  plot(ndata(indx,7)/1000,ndata(indx,8)/1000,'ko','markerfacecolor','r')
%  text((ndata(indx,7))/1000+.25,ndata(indx,8)/1000,num2str(ndata(indx,3)))
%  plot(ndata(indx2,7)/1000,ndata(indx2,8)/1000,'ko','markerfacecolor','b')
%  text((ndata(indx2,7))/1000+.25,ndata(indx2,8)/1000,num2str(ndata(indx2,3)))
%  xlim([-10.5 -3])
%  ylim([2 12])
%  %plot(sns_x,sns_y,'kp','markerfacecolor','k')
%  %text(sns_x,sns_y,sns_n)
%  % Florida map
%  florida_map
% %  xlabel('East(km)')
%  ylabel('North(km)')
%  box on
%  daspect([1,1,1])
%  
%   subplot(4,3,7)
%  hold all
%  %xxx = ndata(indx2,7)
%  plot(ndata(indx,7)/1000,ndata(indx,9)/1000,'ko','markerfacecolor','r')
%  text((ndata(indx,7))/1000+.25,ndata(indx,9)/1000,num2str(ndata(indx,3)))
%  plot(ndata(indx2,7)/1000,ndata(indx2,9)/1000,'ko','markerfacecolor','b')
%  text((ndata(indx2,7))/1000+.25,ndata(indx2,9)/1000,num2str(ndata(indx2,3)))
%  xlim([-10.5 -3])
%  ylim([2 12])
%  %plot(sns_x,sns_y,'kp','markerfacecolor','k')
%  %text(sns_x,sns_y,sns_n)
%  % Florida map
%  florida_map
%  xlabel('East(km)')
%  ylabel('Altitude(km)')
%  box on
%  daspect([1,1,1])
%  
%  
%  subplot(4,3,10)
%  cg = find(strcmp(txtdata(:,3),'CG'));
%  ic = find(strcmp(txtdata(:,3),'IC'));
%  %
%  [frequency, bin] = hist(ndata(cg,9)/1000, 0.5:1:14);
%  [frequency2, bin2] = hist(ndata(ic,9)/1000, 0.5:1:14);
%  size(frequency)
%  size(frequency2)
%  b = bar(bin',[frequency; frequency2]',1,'grouped');
%  set(b(1),'facecolor','b')
%  set(b(2),'facecolor','r')
%  xlim([4,12])
%  ylim([0,8])
%  set(gca,'xgrid','on','ygrid','off')
%  set(gca,'xtick',[3:1:12])
%  %set(gca,'xticklabel',{'3',' ','5',' ','7',' ','9',' ','11',' '})
%  xlabel('Altitude(km)')
%  ylabel('Number of flashes')
%  legend ('CG','IC')
 
 
 
 
 %%  over the land Storm 2 in 2 colors
%  
%  fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011_V04.xlsx';
%  sheet = 4; % over the land t-storm 2
%  xlRange = 'A1:J23'; %over the land storm 1 seperately
%  [ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
% %  ndata
% % size (ndata)
% indx = find(strcmp(txtdata(:,3),'IC'));
% indx2 = find(strcmp(txtdata(:,3),'CG'));
% 
% figure
%     subplot(4,3,2)
%     hold all
%     %xxx = ndata(indx2,7)
%     plot(ndata(indx,7)/1000,ndata(indx,8)/1000,'ko','markerfacecolor','r')
%     %text((ndata(indx,7))/1000+.25,ndata(indx,9)/1000,num2str(ndata(indx,3)))
%     plot(ndata(indx2,7)/1000,ndata(indx2,8)/1000,'ko','markerfacecolor','b')
% %     text((ndata(indx2,7))/1000+.25,ndata(indx2,8)/1000,num2str(ndata(indx2,3)))
%     xlim([-25 40])
%     ylim([-25 40])
%     title('Over the land- storm 2');
% 
%  % plot florida map and sensor sites
%  plot(sns_x,sns_y,'kp','markerfacecolor','k')
%  text(sns_x,sns_y,sns_n)
%  % Florida map
%  florida_map
%  %xlabel('East(km)')
%  ylabel('North(km)')
%  box on
%  daspect([1,1,1])
%  
%  
%  subplot(4,3,5)
%  hold all
%  %xxx = ndata(indx2,7)
%  plot(ndata(indx,7)/1000,ndata(indx,8)/1000,'ko','markerfacecolor','r')
%  text((ndata(indx,7))/1000+.25,ndata(indx,8)/1000,num2str(ndata(indx,3)))
%  plot(ndata(indx2,7)/1000,ndata(indx2,8)/1000,'ko','markerfacecolor','b')
%  text((ndata(indx2,7))/1000+.25,ndata(indx2,8)/1000,num2str(ndata(indx2,3)))
%  xlim([-25 -5])
%  ylim([-15 5])
%  %plot(sns_x,sns_y,'kp','markerfacecolor','k')
%  %text(sns_x,sns_y,sns_n)
%  % Florida map
%  florida_map
% %  xlabel('East(km)')
%  ylabel('North(km)')
%  box on
%  daspect([1,1,1])
%  
%   subplot(4,3,8)
%  hold all
%  %xxx = ndata(indx2,7)
%  plot(ndata(indx,7)/1000,ndata(indx,9)/1000,'ko','markerfacecolor','r')
%  text((ndata(indx,7))/1000+.25,ndata(indx,9)/1000,num2str(ndata(indx,3)))
%  plot(ndata(indx2,7)/1000,ndata(indx2,9)/1000,'ko','markerfacecolor','b')
%  text((ndata(indx2,7))/1000+.25,ndata(indx2,9)/1000,num2str(ndata(indx2,3)))
%  xlim([-25 -5])
%  ylim([2 12])
%  %plot(sns_x,sns_y,'kp','markerfacecolor','k')
%  %text(sns_x,sns_y,sns_n)
%  % Florida map
%  florida_map
%  xlabel('East(km)')
%  ylabel('Altitude(km)')
%  box on
%  daspect([1,1,1])
%  
%  
%  subplot(4,3,11)
%  cg = find(strcmp(txtdata(:,3),'CG'));
%  ic = find(strcmp(txtdata(:,3),'IC'));
%  %
%  [frequency, bin] = hist(ndata(cg,9)/1000, 0.5:1:14);
%  [frequency2, bin2] = hist(ndata(ic,9)/1000, 0.5:1:14);
%  size(frequency)
%  size(frequency2)
%  b = bar(bin',[frequency; frequency2]',1,'grouped');
%  set(b(1),'facecolor','b')
%  set(b(2),'facecolor','r')
%  xlim([3,12])
%  ylim([0,8])
%  set(gca,'xgrid','on','ygrid','off')
%  set(gca,'xtick',[3:1:12])
%  %set(gca,'xticklabel',{'3',' ','5',' ','7',' ','9',' ','11',' '})
%  xlabel('Altitude(km)')
%  ylabel('Number of flashes')
%  legend ('CG','IC')
%  
%  
%  return
 %%  over the sea Storm in 2 colors
%  
 fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011_V04.xlsx';
% 
   sheet = 1; % over the See storm
   xlRange = 'A1:J7'; %over the sea seperately
   [ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
%    %x= txtdata(:,3)
%    %txtdata(:,3)
%  ndata
% size (ndata)
indx = find(strcmp(txtdata(:,3),'IC'));
indx2 = find(strcmp(txtdata(:,3),'CG'));

figure
    subplot(4,3,3)
    hold all
    %xxx = ndata(indx2,7)
    plot(ndata(indx,7)/1000,ndata(indx,8)/1000,'ko','markerfacecolor','r')
    %text((ndata(indx,7))/1000+.25,ndata(indx,9)/1000,num2str(ndata(indx,3)))
    plot(ndata(indx2,7)/1000,ndata(indx2,8)/1000,'ko','markerfacecolor','b')
%     text((ndata(indx2,7))/1000+.25,ndata(indx2,8)/1000,num2str(ndata(indx2,3)))
    xlim([-25 40])
    ylim([-5 40])
    title('Over the sea storm');

 % plot florida map and sensor sites
 plot(sns_x,sns_y,'kp','markerfacecolor','k')
 text(sns_x,sns_y,sns_n)
 % Florida map
 florida_map
 %xlabel('East(km)')
 ylabel('North(km)')
 box on
 daspect([1,1,1])
 
 
 subplot(4,3,6)
 hold all
 %xxx = ndata(indx2,7)
 plot(ndata(indx,7)/1000,ndata(indx,8)/1000,'ko','markerfacecolor','r')
 text((ndata(indx,7))/1000+.25,ndata(indx,8)/1000,num2str(ndata(indx,3)))
 plot(ndata(indx2,7)/1000,ndata(indx2,8)/1000,'ko','markerfacecolor','b')
 text((ndata(indx2,7))/1000+.25,ndata(indx2,8)/1000,num2str(ndata(indx2,3)))
 xlim([25 35])
 ylim([18 28])
 %plot(sns_x,sns_y,'kp','markerfacecolor','k')
 %text(sns_x,sns_y,sns_n)
 % Florida map
 florida_map
%  xlabel('East(km)')
 ylabel('North(km)')
 box on
 daspect([1,1,1])
 
  subplot(4,3,9)
 hold all
 %xxx = ndata(indx2,7)
 plot(ndata(indx,7)/1000,ndata(indx,9)/1000,'ko','markerfacecolor','r')
 text((ndata(indx,7))/1000+.25,ndata(indx,9)/1000,num2str(ndata(indx,3)))
 plot(ndata(indx2,7)/1000,ndata(indx2,9)/1000,'ko','markerfacecolor','b')
 text((ndata(indx2,7))/1000+.25,ndata(indx2,9)/1000,num2str(ndata(indx2,3)))
 xlim([25 35])
 ylim([2 12])
 %plot(sns_x,sns_y,'kp','markerfacecolor','k')
 %text(sns_x,sns_y,sns_n)
 % Florida map
 florida_map
 xlabel('East(km)')
 ylabel('Altitude(km)')
 box on
 daspect([1,1,1])
 
 
 subplot(4,3,12)
 cg = find(strcmp(txtdata(:,3),'CG'));
 ic = find(strcmp(txtdata(:,3),'IC'));
 %
 [frequency, bin] = hist(ndata(cg,9)/1000, 0.5:1:14);
 [frequency2, bin2] = hist(ndata(ic,9)/1000, 0.5:1:14);
 size(frequency)
 size(frequency2)
 b = bar(bin',[frequency; frequency2]',1,'grouped');
 set(b(1),'facecolor','b')
 set(b(2),'facecolor','r')
 xlim([3,10])
 ylim([0,8])
 set(gca,'xgrid','on','ygrid','off')
 set(gca,'xtick',[3:1:10])
 %set(gca,'xticklabel',{'3',' ','5',' ','7',' ','9',' ','11',' '})
 xlabel('Altitude(km)')
 ylabel('Number of flashes')
 legend ('CG','IC')
 
 
 return


   %
 
 %% plot over the land storms seperately
 % Thunderstorm 2
 
fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011_V04.xlsx';
 sheet = 4; % over the land t-storm 2
 xlRange = 'A1:J23'; %over the land storm 1 seperately
 [ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
 ndata
size (ndata)


 figure
 
 %  % Graph option 1 : x (east vs y (north)
  %  %subplot(no of rows (eg 3), no of columns(eg 2), desired position of figure (1:2)),% 1:2 means use space 1 and 2 together
%  
  subplot(3,2,1)
  hold all;
%  % use color map for data
%  %cm contains color values for RGB in 3 columns 
%  % cm has 64 colors , so L is 64
  cm = jet;
  L = length(cm); % 64
%  % data will be colored according to time pbfat
%  %pbfat is a matrix 
  pbfat = ndata(:,6);
%  % cl matrix will hold all the color number data for pbfat (following equation for cl is derived using gradient
%  % method to calculate unknown color value for a known time pbfat)
  cl = 1+ round((L-1) *(pbfat - pbfat(1))/(pbfat(end) - pbfat(1)));
% 
 for i= 1:23
%  
  pbfax = ndata(i,7)/1000;
  pbfay = ndata(i,8)/1000;
  pbfaz = ndata(i,9)/1000;
% 
  plot(pbfax,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
  %text(pbfax+0.25,pbfay,num2str(i))
 end
% plot([20 35 35 20 20], [35 35 20 20 35],'k--');
% 
% 
% % plot florida map and sensor sites
  plot(sns_x,sns_y,'kp','markerfacecolor','k')
  text(sns_x,sns_y,sns_n)
%  % Florida map
  florida_map
  xlim([-30, -5])
  ylim([-15, 10])
  xlabel('East(km)')
  ylabel('North(km)')
  box on
  daspect([1,1,1])
  %set(gca,'position',[0.0360    0.7706    0.7750    0.2157])
  
  subplot(3,2,2)
  hold all;
  for i= 1:23
%  
  pbfax = ndata(i,7)/1000;
  pbfay = ndata(i,8)/1000;
  pbfaz = ndata(i,9)/1000;
  pbfat = ndata(i,6);
% 
  plot(pbfax,pbfay,'ro','markersize',5)
  %text(pbfax+0.25,pbfay,num2str(i))
 end
  
% 
%  % Graph option 2: x (east) vs Z (altitude)
% 
  subplot(3,2,3)
  hold all;
  for i= 1:23
  pbfat = ndata(i,6);
  pbfax = ndata(i,7)/1000;
  pbfay = ndata(i,8)/1000;
  pbfaz = ndata(i,9)/1000;
%  indx = ndata(i,4)
%  % color ceded data will be plotted cm(raw, column), seperate columns for
%  % RGB , all we want to give is row number cl(i)
  plot(pbfax,pbfaz,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
%  % numbering the data point
  text(pbfax,pbfaz,num2str(i))
%   
 end
  xlim([-20, -14])
  ylim([3, 12])
  xlabel('East(km)')
  ylabel('Altitude(km)')
  box on
  daspect([1,1,1])
  %set(gca,'position',[0.0830    0.4976    0.3347    0.2157])
% 
%   % Graph option 3: Histogram
%  
   subplot(3,2,4)
%   hist(ndata(:,9)/1000,0.5:1:14)
%   xlim([3,12])
%   ylim([0,12])
%   set(gca,'xgrid','on','ygrid','off')
%   set(gca,'xtick',[3:1:12])
%   set(gca,'xticklabel',{'3',' ','5',' ','7',' ','9',' ','11',' '})
%   xlabel('Altitude(km)')
%   ylabel('Number of flashes')
%   %set(gca,'position',[0.5374    0.5127    0.2902    0.2153])
  
  %[frequency, bin]
  
   %disp(txtdata(:,1))
  
  cg = find(strcmp(txtdata(:,1),'CG'));
  ic = find(strcmp(txtdata(:,1),'IC'));
  
   [frequency, bin] = hist(ndata(cg,9)/1000, 0.5:1:14);
   [frequency2, bin2] = hist(ndata(ic,9)/1000, 0.5:1:14);
%    size(frequency)
%    size(frequency2)
   bar(bin',[frequency; frequency2]',1,'grouped')
   xlim([3,12]);
   ylim([0,10]);
   set(gca,'xgrid','on','ygrid','off');
   set(gca,'xtick',[3:1:12]);
   %set(gca,'xticklabel',{'3',' ','5',' ','7',' ','9',' ','11',' '})
   xlabel('Altitude(km)');
   ylabel('Number of flashes');
   legend ('CG','IC');
  
  
  
%  
%  
% % Graph option 4 : x (east) vs y (north)
% % setup
%  
  subplot(3,2,5)
  hold all;
  cm = jet;
  L = length(cm);
  pbfat = ndata(:,6);
  cl = 1+ round(63 *(pbfat - pbfat(1))/(pbfat(end) - pbfat(1)));
% 
 for i= 1:23
%  
  pbfax = ndata(i,7)/1000;
  pbfay = ndata(i,8)/1000;
  pbfaz = ndata(i,9)/1000;
% 
   plot(pbfax,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
   %text(pbfax+0.35,pbfay,num2str(i))
   text(pbfax,pbfay,num2str(i))
 end
% 
% % plot florida map and sensor sites
  plot(sns_x,sns_y,'kp','markerfacecolor','k')
  text(sns_x,sns_y,sns_n)
%  % Florida map
  florida_map
  xlim([-20, -14])
  ylim([-12, -2])
  xlabel('East(km)')
  ylabel('North(km)')
  box on
  daspect([1,1,1])
  %set(gca,'position',[-0.0376    0.0359    0.5737    0.4271])
%  
% 
% % Graph option 5: z (altitude) vs y (north)
% % 
  subplot(3,2,6)
  hold all;
  for i= 1:23
  pbfat = ndata(i,6);
  pbfax = ndata(i,7)/1000;
  pbfay = ndata(i,8)/1000;
  pbfaz = ndata(i,9)/1000;
%  indx = ndata(i,4)
%  
  plot(pbfaz,pbfay,'o','color',cm(cl(i),:),'markerfacecolor',cm(cl(i),:),'markersize',5)
  text(pbfaz,pbfay,num2str(i))
%  
%  
  end
  xlim([3, 12])
  ylim([-12,-2])
  xlabel('Altitude(km)')
  ylabel('North(km)')
  box on
  daspect([1,1,1])
  set(gca,'position',[0.4120    0.0417    0.5394    0.4276])
 


%% alternative method to plot florida map and E change sensor sites.
 %sen_set = open('sensor_setting.mat');
% hold all;
% plot(sen_set.x(1:10)/1000,sen_set.y(1:10)/1000,'k*');
% text(sen_set.x(1:10)/1000,sen_set.y(1:10)/1000,sen_set.sen_IDs{1:10});
% daspect([1 1 1]);



%% Thunderstorm 1 IC and CG color coded by type and Z-T (altitude, time)plot
%  fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011_V04.xlsx';
% % 
%    sheet = 3; % over the land t-storm 1
%    xlRange = 'A1:J24'; %over the land storm 1 seperately
%    [ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
%    %x= txtdata(:,3)
%    
%   
%     indx = find(strcmp(txtdata(:,3),'IC'));
%     indx2 = find(strcmp(txtdata(:,3),'CG'));
%     % ndata(indx,3)
%     figure
%     plot(ndata(indx,6),ndata(indx,9)/1000,'r*')
%     text((ndata(indx,6))+25,ndata(indx,9)/1000,num2str(ndata(indx,3)))
%     hold all
%     plot(ndata(indx2,6),ndata(indx2,9)/1000,'b*')
%     text((ndata(indx2,6))+25,ndata(indx2,9)/1000,num2str(ndata(indx2,3)))
%     xlabel('Time in seconds');
%     ylabel('Altitude (km)')
%     legend('IC','CG')
%     x= ndata(:,6)
   
   %% Thunderstorm 2 IC and CG color coded by type and Z-T (altitude, time)plot
% fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011_V04.xlsx';
% 
%    sheet = 4; % over the land t-storm 2
%    xlRange = 'A1:J23'; %over the land storm 2 seperately
%    [ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
%    %x= txtdata(:,3)
%    %txtdata(:,3)
%   
%     indx = find(strcmp(txtdata(:,3),'IC'));
%     indx2 = find(strcmp(txtdata(:,3),'CG'));
%     % ndata(indx,3)
%     figure
%     plot(ndata(indx,6),ndata(indx,9)/1000,'r*')
%     text((ndata(indx,6))+25,ndata(indx,9)/1000,num2str(ndata(indx,3)))
%     hold all
%     plot(ndata(indx2,6),ndata(indx2,9)/1000,'b*')
%     text((ndata(indx2,6))+25,ndata(indx2,9)/1000,num2str(ndata(indx2,3)))
%     xlabel('time in seconds');
%     ylabel('Altitude (km)')
%     legend('IC','CG')
%    x= ndata(:,6)

 %% Thunderstorm 3 IC and CG color coded by type and Z-T (altitude, time)plot
% fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011_V04.xlsx';
% 
%    sheet = 1; % over the See storm
%    xlRange = 'A1:J7'; %over the sea seperately
%    [ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
%    %x= txtdata(:,3)
%    %txtdata(:,3)
%   
%     indx = find(strcmp(txtdata(:,3),'IC'));
%     indx2 = find(strcmp(txtdata(:,3),'CG'));
%     % ndata(indx,3)
%     figure
%     plot(ndata(indx,6),ndata(indx,9)/1000,'r*')
%     text((ndata(indx,6))+25,ndata(indx,9)/1000,num2str(ndata(indx,3)))
%     hold all
%     plot(ndata(indx2,6),ndata(indx2,9)/1000,'b*')
%     text((ndata(indx2,6))+25,ndata(indx2,9)/1000,num2str(ndata(indx2,3)))
%     xlabel('time in seconds');
%     ylabel('Altitude (km)')
%     legend('IC','CG')

%% Long time LDAR plot( to create a plot for more than 5 mins with ,LDAR,CGLSS and PBFA)

% User inputs: date(uopdate plotter date), time,sheetno, xlrange, LDAR box
%check LDAR box in plotter and update it as per your storm
% The date coming from Plotter2. Don't worry about hh, mm, ss.

% fn = 'C:\Users\Sumedhe\Desktop\Nadee\Nadee-official-sumpc\Origination-of-flashes-in-samespot\Flash-origination-07222011_V04.xlsx';
  

%%Over the land storm 1
% addust LDAR box in plotter to ( x = -10, y =0, xlength ->15, y length -20)
% sheet = 3; % over the land 1
%xlRange = 'A1:J24'; %for over the land1
%[ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
% t1 = 59297;
% t2 = 61942;

%%Over the land storm 2
%check LDAR box in plotter and update it to ( x = -22, y =-15, xlength -> 15, y length -> 15)
% The date coming from Plotter2. Don't worry about hh, mm, ss.

% sheet = 4; % over the land 2
% xlRange = 'A1:J23'; %for over the land2
% [ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
% t1 = 61382;
% t2 = 63495;

%%Over the sea storm 1
%check LDAR box in plotter and update it to ( x = 20, y =20, xlength -> 20, y length -> 20)
% The date coming from Plotter2. Don't worry about hh, mm, ss.

% sheet = 1; % over the sea
% xlRange = 'A1:J7'; %for over sea
% [ndata, txtdata, data] = xlsread(fn, sheet, xlRange);
% t1 = 78303;
% t2 = 79309.5;

% Start the program
% try h=guidata(findall(0,'Tag','plotter2'));
% catch; disp('Run plotter2 first!'); return;
% end
% 
% g = h.g;
% 
% tSt = floor(t1/(30*60))*(30*60);
% tEn = floor(t2/(30*60))*(30*60);
% ts = tSt:30*60:tEn;
% L = length(ts);
% 
% g.t1 = t1;
% g.t2 = t2;
% sen_set = h.sen_set;
% loc_dir = sen_set.loc_dir;
% 
% dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
%     -datenum(str2double(g.YYYY),0,0));
% 
% x0 = 0; y0 = 0; z0 = 0;
% DLSs = [];
% CGs = [];
% Pt = [];
% Pz = [];
% 
% for t = ts
%     [~, g.hh ,g.mm] = sec2hhmmss(t);
%    
% 
%     if g.mm < 30
%         ext=0;
%     else
%         ext=30;
%     end
%     
%     ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
%         loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);
%     
%     
%     % load data
%     [CG,~,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(sen_set.ldar_r),...
%          x0,y0,z0,0);
%      
%      DLSs = [DLSs; DLS];
%      CGs = [CGs; CG];
%      
%      
%      pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
%         loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
%         g.MM{:},g.DD{:},g.hh,ext);
%     
%     try
%         PBFA = pbfaExtract(pbfa_fn,g.t1,g.t2,str2double(sen_set.ldar_r),...
%             x0,y0,z0);
%         if ~isnan(PBFA(1))
%             Pt = [Pt; PBFA(:,6)];
%             Pz = [Pz; PBFA(:,5)/1000];
%         end
%     end
%     
%     pbfa_fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
%         loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
%         g.MM{:},g.DD{:},g.hh,ext);
%     
%     try
%         PBFA = pbfaExtract(pbfa_fn,g.t1,g.t2,str2double(sen_set.ldar_r),...
%             x0,y0,z0);
%         if ~isnan(PBFA(1))
%             Pt = [Pt; PBFA(:,6)];
%             Pz = [Pz; PBFA(:,5)/1000];
%         end
%     end
% 
% 
% end
% 
% figure
% hold all
% plot(DLSs(:,10),DLSs(:,8)/1000,'ko','MarkerFaceColor','k','MarkerSize',3)
% plot(CGs(:,10),CGs(:,8)/1000,'ks','MarkerFaceColor','g','MarkerSize',5)
% plot(Pt,Pz,'ro','MarkerFaceColor','r','MarkerSize',3)
% % over the land storm 1 , begining points
% cg = find(strcmp(txtdata(:,3),'CG'));
% ic = find(strcmp(txtdata(:,3),'IC'));
% 
% 
% 
% plot(ndata(ic,6), ndata(ic,9)/1000,'rs')
% plot(ndata(cg,6), ndata(cg,9)/1000,'bs')

% for i=1:length(ndata(:,6))
% plot([0,0]+ndata(i,6),ylim,'k--')
% end
% box on
% xlabel('Seconds from mid-night')
% ylabel('Altitude (km)')
% title(sprintf('%s -- %s UT',sec2hhmmss(t1),sec2hhmmss(t2)))
% legend('LDAR','CGLSS','PBFA','IC initiation','CG initiation')
% 



%% Plot vertical radar data- Over the land storm 1
% This program is for generating the same line(AA') in different radar PPI , so
% tht vertical radar plane generated will be the same for each time period
% copy  and paste x'y'z point of the flash origination from excell sheeet
% per each 5 min period

%first using the plotter generate the PPI for the desired 'time' and 'angle' -we use angle 5(3.2 deg)here
 % User input
%    a = guidata(gcf);
%   fn = a.fileName;
%   hold all
%   xlim([-20 0]);
%   ylim([-5 20]);
% 
%  %***********IC/CG flash- time 16.25-16.30***************
% % xyz_CG = [29662.4	22075.9	4755.4]; no CG during this time
%  xyz_IC = [1 -7976	3568	7922
%            4 -7467	4277	8345
%            5 -9110	5029	9237
%            7 -8860	3655	9685
%            8 -7029	4458	10338
%             ];
% %         
% %    
% %         
% % % %cordinates for over the sea flash in '2011-07-22-Originationof flash project'
% % % %line AA'-option 1
% x1= -10.8395;
% x2= -3.7778;
% y1= 9.4444;
% y2= -2.9012;
% % x1= -11.5299;
% % x2= -1.8284;
% % y1= 2.4254;
% % y2= 6.6791;
% 
% %sumedhe
% x1= -30.3189
% x2= -0.5239
% y1= -7.0387
% y2= 11.8223
% 
% 
% % %line BB'-option1
% % x1b= -9.7910;
% % x2b= -3.8209;
% % y1b= 11.6269;
% % y2b= -1.6269;
% 
% x1b= -12.6493;
% x2b= -4.6642;
% y1b= 19.6642;
% y2b= -4.2164;
% % 
% %sumedhe
% x1b= -27.8588
% x2b= -3.5308
% y1b= 23.5763
% y2b= -13.7358
% 
% hold all
%   %for CG
% %   plot(xyz_CG(:,1)/1000,xyz_CG(:,2)/1000,'ko','markerfacecolor','b','markersize',6)
%   %for IC
%   plot(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)
%   text(xyz_IC(:,2)/1000+0.2,xyz_IC(:,3)/1000+0.2,num2str(xyz_IC(:,1)))
% % addfm
%   plot([x1 x2],[y1 y2],'ks-');
%   plot([x1b x2b],[y1b y2b],'ks--');
%  
% %  % Plot data
% % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%  plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(2)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
%  r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6);
%  text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
%  
%  plot_locations=0;
% %   %to generate RHI on BB
%  plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
%  radarTools(1)
% delete(fID)
% r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_IC = convert2D([x1b y1b x2b y2b],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6);
%  text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))

%************ %***********IC/CG flash- time 16.30-16.35***************
% % % no IC for this range
% xyz_CG =[ 9   -7808	  2942	7095
%           11  -8020	  3028	8074 ];
%       
%   % %line AA'-option1
% % x1= -12.5373
% % x2= -4.7761
% % y1= 10.6119
% % y2= -1.9254
% 
% x1= -11.4552;
% x2= -1.6045;
% y1= 1.3806;
% y2= 5.7090;
% 
% 
% % %line BB'-option1
% % x1b= -9.7910
% % x2b= -3.8209
% % y1b= 11.6269
% % y2b= -1.6269
% 
% x1b= -12.6493;
% x2b= -4.6642;
% y1b= 19.6642;
% y2b= -4.2164;
% 
% hold all
%   
%   %for CG
%   plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
%   text(xyz_CG(:,2)/1000+0.2,xyz_CG(:,3)/1000+0.2,num2str(xyz_CG(:,1)))
%   plot([x1 x2],[y1 y2],'ks-')
%   plot([x1b x2b],[y1b y2b],'ks--')
%addfm
%  
% %  % Plot data
% % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%  plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_CG = convert2D([x1b y1b x2b y2b],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))    



%************ %***********IC/CG flash- time 16.35-16.40***************
%no CG or IC in this period
% this matrix is just used as a dummy to generate graph witout errors
%  xyz_CG =[ 9   -7808	  2942	7095
%            11  -8020	  3028	8074 ];
% 
% x1= -11.0821;
% x2= -1.1567;
% y1= -0.3358;
% y2= 4.1418;
% 
% x1b= -12.6493;
% x2b= -4.6642;
% y1b= 19.6642;
% y2b= -4.2164;
% 
% hold all
% 
% plot([x1 x2],[y1 y2],'ks-')
% plot([x1b x2b],[y1b y2b],'ks--')
% 
% 
% plot_locations = 0;
% % %to generate RHIon AA
% plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_CG = convert2D([x1b y1b x2b y2b],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))  




%************ %***********IC/CG flash- time 16.40-16.45***************

% xyz_IC = [17 -6120	3300	8581
%            18 -6324	3494	9413
%            20 -6425	3671	8483
%            23 -6409	3298	8711
%            24 -6589	3360	8738
%            26 -7858	3731	8764];
% 
% 
% xyz_CG =[  19 -6420	4006	6071
%            22 -6696	3804	4559
%            25 -6063.5	3994	5451.4
%            27 -6292	4442	5157 ];
%       
%   % %line AA'-option 1
% % x1= -12.5373
% % x2= -4.7761
% % y1= 10.6119
% % y2= -1.9254
% 
% x1= -11.3433;
% x2= -1.7279;
% y1= 1.6035;
% y2= 5.9587;
% 
% % %line BB'-option 1
% % x1b= -9.7910
% % x2b= -3.8209
% % y1b= 11.6269
% % y2b= -1.6269
% 
% 
% x1b= -12.6493;
% x2b= -4.6642;
% y1b= 19.6642;
% y2b= -4.2164;
% 
% hold all
%   
%   %for CG
%   plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
%   text(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,num2str(xyz_CG(:,1)))
%   plot(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)
%   text(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,num2str(xyz_IC(:,1)))
%   plot([x1 x2],[y1 y2],'ks-')
%   plot([x1b x2b],[y1b y2b],'ks--')
%addfm
%  
% %  % Plot data
% % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%  plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6)
% text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
%  
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% r_IC = convert2D([x1b y1b x2b y2b],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% r_CG = convert2D([x1b y1b x2b y2b],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6)
% text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))    

 
 %************ %***********IC/CG flash- time 16.45-16.50***************

% % no IC for this range
% xyz_CG =[  28  -7029.4	4211.6	5347.8
%            29  -7483	5102	5728
%            30  -6299.8	5234.3	6051.4
%             ];
%       
%   % %line AA'-option1
% % x1= -12.5373
% % x2= -4.7761
% % y1= 10.6119
% % y2= -1.9254
% 
% x1= -11.7537;
% x2= -2.4254;
% y1= 2.7239;
% y2= 6.9030;
% 
% % %line BB'-option1
% % x1b= -9.7910
% % x2b= -3.8209
% % y1b= 11.6269
% % y2b= -1.6269
% 
% 
% x1b= -12.6493;
% x2b= -4.6642;
% y1b= 19.6642;
% y2b= -4.2164;
% 
% hold all
%   
%   %for CG
%   plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
%   text(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,num2str(xyz_CG(:,1)))
% %   plot(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)
% %   text(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,num2str(xyz_IC(:,1)))
%   plot([x1 x2],[y1 y2],'ks-')
%   plot([x1b x2b],[y1b y2b],'ks--')
%addfm
%  
% %  % Plot data
% % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%  plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6)
% % text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
%  
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_IC = convert2D([x1b y1b x2b y2b],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% r_CG = convert2D([x1b y1b x2b y2b],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6)
% % text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))    
% 

%    

%************ %***********IC/CG flash- time 16.50-16.55***************
% No IC


% xyz_CG =[ 31 -6346	5657	5233
%           32 -7434	6234	5584 ];
      
  % %line AA'-option 1
% x1= -12.5373
% x2= -4.7761
% y1= 10.6119
% y2= -1.9254

% x1= -12.1269;
% x2= -2.7985;
% y1= 3.6194;
% y2= 7.8731;



% %line BB'-option 1
% x1b= -9.7910
% x2b= -3.8209
% y1b= 11.6269
% y2b= -1.6269

% x1b= -12.6493;
% x2b= -4.6642;
% y1b= 19.6642;
% y2b= -4.2164;
% 
% 
% hold all
%   
%   %for CG
%   plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
%   text(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,num2str(xyz_CG(:,1)))
% %   plot(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)
% %   text(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,num2str(xyz_IC(:,1)))
%   plot([x1 x2],[y1 y2],'ks-')
%   plot([x1b x2b],[y1b y2b],'ks--')
%addfm
%  
% %  % Plot data
% % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%  plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6)
% % text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
%  
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_IC = convert2D([x1b y1b x2b y2b],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% r_CG = convert2D([x1b y1b x2b y2b],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6)
% % text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))    

%           
%           
%           

%************ %***********IC/CG flash- time 16.55-17.00***************
% No flashes
%************ %***********IC/CG flash- time 17.00-17.05***************
% No flashes
% %For time 17.00-17.05
% x1= -13.1716
% x2= -2.2761
% y1= 5.6343
% y2= 10.5597

% %************ %***********IC/CG flash- time 17.05-17.10***************
% % no CG
% xyz_IC =[ 38 -7382	7957	6355.7
%           ];
% %       
%   % %line AA'-option 1
% % x1= -12.5373
% % x2= -4.7761
% % y1= 10.6119
% % y2= -1.9254
% 
% x1= -13.0224;
% x2= -2.2015;
% y1= 5.7836;
% y2= 10.1866;
% 
% % %line BB'-option 1
% % x1b= -9.7910
% % x2b= -3.8209
% % y1b= 11.6269
% % y2b= -1.6269
% 
% x1b= -12.6493;
% x2b= -4.6642;
% y1b= 19.6642;
% y2b= -4.2164;
% 
% hold all
%   
%   %for CG
% %   plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% %   text(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,num2str(xyz_CG(:,1)))
%   plot(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)
%   text(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,num2str(xyz_IC(:,1)))
%   plot([x1 x2],[y1 y2],'ks-')
%   plot([x1b x2b],[y1b y2b],'ks--')
%addfm
%  
% %  % Plot data
% % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%  plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
%  r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% %r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6)
% text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
%  
% % plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
% %  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
%  r_IC = convert2D([x1b y1b x2b y2b],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% % r_CG = convert2D([x1b y1b x2b y2b],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6)
% text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
% % plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
% %  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))    
% % 
% %    

% %************ %***********IC/CG flash- time 17.10-17.12***************
% % no CG
% xyz_IC =[ 48 -5963	9908	5551
%           ];
% %       
%   % %line AA'-option1
% % x1= -12.5373
% % x2= -4.7761
% % y1= 10.6119
% % y2= -1.9254
% 
% x1= -13.3512;
% x2= -2.4632;
% y1= 6.9485;
% y2= 11.3320;
% 
% % %line BB'-option1
% % x1b= -9.7910
% % x2b= -3.8209
% % y1b= 11.6269
% % y2b= -1.6269
% 
% x1b= -11.9030;
% x2b= -4.2164;
% y1b= 19.5149;
% y2b= -2.6493;
% 
% 
% 
% hold all
%   
%   %for CG
% %   plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% %   text(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,num2str(xyz_CG(:,1)))
%   plot(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)
%   text(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,num2str(xyz_IC(:,1)))
%   plot([x1 x2],[y1 y2],'ks-')
%   plot([x1b x2b],[y1b y2b],'ks--')
%  addfm
% %  % Plot data
% % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%  plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
%  r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% %r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6)
% text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
%  
% % plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
% %  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
%  r_IC = convert2D([x1b y1b x2b y2b],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% % r_CG = convert2D([x1b y1b x2b y2b],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6)
% text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
% % plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
% %  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))    



%% Plot vertical radar data- Over the land storm 2
% This program is for generating the same line(AA') in different radar PPI , so
% tht vertical radar plane generated will be the same for each time period
% copy  and paste x'y'z point of the flash origination from excell sheeet
% per each 5 min period

%first using the plotter generate the PPI for the desired 'time' and 'angle' -we use angle 5(3.2 deg)here
 % User input
%    a = guidata(gcf);
%   fn = a.fileName;
%   hold all
%   %storm 2
%   xlim([-25 5]);
%   ylim([-15 20]);

% % % 
% % %  %***********IC/CG flash- time 17.00-17.05***************
% % % % xyz_CG = [29662.4	22075.9	4755.4]; no CG during this time
%   xyz_IC = [33 -16681	-7789	10138
%             34 -16499	-6473	9520
%             35 -17116	-5958	8082
%                         ];
% 
%    
%         
% % %cordinates for over the sea flash in '2011-07-22-Originationof flash project'
%  %line AA'
% 
%  x1 = -21.2313;
%  x2 = -11.828;
%  y1= -4.9179;
%  y2= -9.3060;
% % 
%  
% %line BB'
% 
% x1b= -19.6861;
% x2b= -14.9481;
% y1b= -14.4078;
% y2b= -1.5311;
% 
% % hold all
% %   %for CG
% % %   plot(xyz_CG(:,1)/1000,xyz_CG(:,2)/1000,'ko','markerfacecolor','b','markersize',6)
% %   %for IC
%    plot(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)
%    text(xyz_IC(:,2)/1000+0.2,xyz_IC(:,3)/1000+0.2,num2str(xyz_IC(:,1)))
%     plot([x1 x2],[y1 y2],'ks-');
% 
% %    
%   
%    
%    plot([x1b x2b],[y1b y2b],'ks--');
%  
% %  % Plot data
% % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%   plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6);
%  text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
%  
% plot_locations=0;
% %   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(1)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_IC = convert2D([x1b y1b x2b y2b],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6);
%  text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))

%  %***********IC/CG flash- time 17.05-17.10***************
% % xyz_IC = [29662.4	22075.9	4755.4]; no IC during this time
%   xyz_CG = [36 -16822.7	-5886.8	7546
%             37 -16817	-7198	6429
%             39 -16902	-6943	6708
%             40 -16885	-5193	6711
%             41 -16408	-6115	6804
%             42 -16438	-5686	6484
%             43 -15908	-5838	5782
%             44 -16650.8	-5338.6	6714
%             45 -16230	-5764	5966
%                       
%             ];
% 
%          
% % % %cordinates for over the sea flash in '2011-07-22-Originationof flash project'
%  %line AA'
% 
%  x1= -21.2313;
%  x2= -11.828;
%  y1= -4.9179;
%  y2= -9.3060;
% % 
%  %AA-opion 2
% % x1= -20.6032;
% %  x2= -11.4329;
% %  y1= -3.5180;
% %  y2= -8.3324;
% 
% %line BB'
% 
% x1b= -19.6861;
% x2b= -14.9481;
% y1b= -14.4078;
% y2b= -1.5311;
% % 
% % hold all
% %   %for CG
% % %   plot(xyz_CG(:,1)/1000,xyz_CG(:,2)/1000,'ko','markerfacecolor','b','markersize',6)
% %   %for IC
%    plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
%    text(xyz_CG(:,2)/1000+0.2,xyz_CG(:,3)/1000+0.2,num2str(xyz_CG(:,1)))
%     plot([x1 x2],[y1 y2],'ks-');
% 
% %      
%    plot([x1b x2b],[y1b y2b],'ks--');
% %  
% % %  % Plot data
% % % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%   plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6);
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_CG = convert2D([x1b y1b x2b y2b],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6);
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))


%  %***********IC/CG flash- time 17.10-17.15***************
% xyz_IC = [29662.4	22075.9	4755.4]; no IC during this time
%   xyz_CG = [46 -16456.5	-5454	6503.8
%             47 -16312.9	-5677.8	5052.7
%             ];
% 
%          
% % % %cordinates for over the sea flash in '2011-07-22-Originationof flash project'
%  %line AA'
%  
%  
% 
%  x1= -22.4754;
%  x2= -12.0060;
%  y1= -2.2571;
%  y2= -7.8357;
% % 
%  %AA-opion 2
% % x1= -20.6032;
% %  x2= -11.4329;
% %  y1= -3.5180;
% %  y2= -8.3324;
% 
% %line BB'
% 
% x1b= -19.6861;
% x2b= -14.9481;
% y1b= -14.4078;
% y2b= -1.5311;
% % 
% % hold all
% %   %for CG
% % %   plot(xyz_CG(:,1)/1000,xyz_CG(:,2)/1000,'ko','markerfacecolor','b','markersize',6)
% %   %for IC
%    plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
%    text(xyz_CG(:,2)/1000+0.2,xyz_CG(:,3)/1000+0.2,num2str(xyz_CG(:,1)))
%     plot([x1 x2],[y1 y2],'ks-');
% 
% %      
%    plot([x1b x2b],[y1b y2b],'ks--');
% %  
% % %  % Plot data
% % % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%   plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6);
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_CG = convert2D([x1b y1b x2b y2b],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6);
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))

 %  %***********IC/CG flash- time 17.15-17.20***************
% % xyz_IC = [29662.4	22075.9	4755.4]; no IC during this time
%   xyz_CG = [49	-16116.5 -3986.2	7966.2
% 	                    ];
         
 % %cordinates for over the sea flash in '2011-07-22-Originationof flash project'
%  %line AA'
% 
%  x1= -22.0672;
%  x2= -11.3060;
%  y1= -0.3209;
%  y2= -6.6940;
% % 
%  %AA-opion 2
% % x1= -20.6032;
% %  x2= -11.4329;
% %  y1= -3.5180;
% %  y2= -8.3324;
% 
% %line BB'
% 
% x1b= -19.6861;
% x2b= -14.9481;
% y1b= -14.4078;
% y2b= -1.5311;
% % 
% % hold all
% %   %for CG
% % %   plot(xyz_CG(:,1)/1000,xyz_CG(:,2)/1000,'ko','markerfacecolor','b','markersize',6)
% %   %for IC
%    plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
%    text(xyz_CG(:,2)/1000+0.2,xyz_CG(:,3)/1000+0.2,num2str(xyz_CG(:,1)))
%     plot([x1 x2],[y1 y2],'ks-');
% 
% %      
%    plot([x1b x2b],[y1b y2b],'ks--');
% %  
% % %  % Plot data
% % % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%   plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6);
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_CG = convert2D([x1b y1b x2b y2b],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6);
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))

 
 
 %  %***********IC/CG flash- time 17.20-17.25***************
% % xyz_IC = [29662.4	22075.9	4755.4]; no IC during this time
%   xyz_CG = [50  -15732	-3574	6651
% 	                    ];
  

% % % %cordinates for over the sea flash in '2011-07-22-Originationof flash project'
%  %line AA'
% 
%  x1= -22.9030;
%  x2= -11.5149;
%  y1= -0.0075;
%  y2= -5.6493;
% % 
%  %AA-opion 2
% % x1= -20.6032;
% %  x2= -11.4329;
% %  y1= -3.5180;
% %  y2= -8.3324;
% 
% %line BB'
% 
% x1b= -19.6861;
% x2b= -14.9481;
% y1b= -14.4078;
% y2b= -1.5311;
% % 
% % hold all
% %   %for CG
% % %   plot(xyz_CG(:,1)/1000,xyz_CG(:,2)/1000,'ko','markerfacecolor','b','markersize',6)
% %   %for IC
%    plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
%    text(xyz_CG(:,2)/1000+0.2,xyz_CG(:,3)/1000+0.2,num2str(xyz_CG(:,1)))
%     plot([x1 x2],[y1 y2],'ks-');
% 
% %      
%    plot([x1b x2b],[y1b y2b],'ks--');
% %  
% % %  % Plot data
% % % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%   plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6);
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% % r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_CG = convert2D([x1b y1b x2b y2b],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% hold all
% % plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6);
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))

  %  %***********IC/CG flash- time 17.25-17.30***************
%  xyz_IC = [53  -19150	-11471	7952
%            54  -18499	-10561	5158]; 
%   xyz_CG = [55  -18278.2	-10422.1	5576.4
%  	        56  -17964	    -11306.9	4639.4
% 	                    ];
% 
 
% % % %cordinates for over the sea flash in '2011-07-22-Originationof flash project'
%  %line AA'
% 
%  x1= -22.0672;
%  x2= -13.1866;
%  y1= -9.2015;
%  y2= -13.0672;
% % 
%  %AA-opion 2
% % x1= -20.6032;
% %  x2= -11.4329;
% %  y1= -3.5180;
% %  y2= -8.3324;
% 
% %line BB'
% 
% x1b= -19.6861;
% x2b= -14.9481;
% y1b= -14.4078;
% y2b= -1.5311;
% % 
% % hold all
% %   %for CG
%    plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
%    text(xyz_CG(:,2)/1000+0.2,xyz_CG(:,3)/1000+0.2,num2str(xyz_CG(:,1)))
% %   %for IC
%    plot(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)
%    text(xyz_IC(:,2)/1000+0.2,xyz_IC(:,3)/1000+0.2,num2str(xyz_IC(:,1)))
%     plot([x1 x2],[y1 y2],'ks-');
% 
% %      
%    plot([x1b x2b],[y1b y2b],'ks--');
% %  
% % %  % Plot data
% % % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%   plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
%  r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% hold all
%  plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6);
%  text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
%  r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% r_IC = convert2D([x1b y1b x2b y2b],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% hold all
%  plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6);
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
% 
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6);
%  text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
% 
%  
 
 %  %***********IC/CG flash- time 17.30-17.35***************
% %  xyz_IC = [53  -19150	-11471	7952 % no IC
% %            54  -18499	-10561	5158]; 
%   xyz_CG = [57  -18032.2	-9739.9 	4783
%  	        58  -17991	    -9445	    5065
% 	                    ];



        
% % %cordinates for over the sea flash in '2011-07-22-Originationof flash project'
 %line AA'

%  x1= -23.3160;
%  x2= -13.1905;
%  y1= -8.7145;
%  y2= -13.1086;
% 
 %AA-opion 2
%  x1= -22.9030;
%   x2= -12.4552;
%   y1= -7.2164;
%   y2= -12.0224;
% 
% %line BB'
% 
% x1b= -19.6861;
% x2b= -14.9481;
% y1b= -14.4078;
% y2b= -1.5311;
% % 
% % hold all
% %   %for CG
%    plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
%    text(xyz_CG(:,2)/1000+0.2,xyz_CG(:,3)/1000+0.2,num2str(xyz_CG(:,1)))
% %   %for IC
% %    plot(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)
% %    text(xyz_IC(:,2)/1000+0.2,xyz_IC(:,3)/1000+0.2,num2str(xyz_IC(:,1)))
%      plot([x1 x2],[y1 y2],'ks-');
% 
% %      
%    plot([x1b x2b],[y1b y2b],'ks--');
% %  
% %  % Plot data
% % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%   plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
%  r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% % r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% hold all
%  plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
% % plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6);
% %  text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
%  r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
% % r_IC = convert2D([x1b y1b x2b y2b],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% hold all
%  plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6);
%  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))

% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6);
%  text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))

 
%  %***********IC/CG flash- time 17.35-17.40***************
%   xyz_IC = [59  -17966	-7755	6345 
%             ]; 
% %   xyz_CG = [55  -18278.2	-10422.1	5576.4 % no CG
% %  	        56  -17964	    -11306.9	4639.4
% % 	                    ];
%        
% % % %cordinates for over the sea flash in '2011-07-22-Originationof flash project'
%  %line AA'
% 
%  x1= -24.6152;
%  x2= -13.0759;
%  y1= -5.2374;
%  y2= -9.5933;
% % 
%  %AA-opion 2
% % x1= -20.6032;
% %  x2= -11.4329;
% %  y1= -3.5180;
% %  y2= -8.3324;
% 
% %line BB'
% 
% x1b= -19.6861;
% x2b= -14.9481;
% y1b= -14.4078;
% y2b= -1.5311;
% % 
% % hold all
% %   %for CG
% %    plot(xyz_CG(:,2)/1000,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% %    text(xyz_CG(:,2)/1000+0.2,xyz_CG(:,3)/1000+0.2,num2str(xyz_CG(:,1)))
% %   %for IC
%    plot(xyz_IC(:,2)/1000,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)
%    text(xyz_IC(:,2)/1000+0.2,xyz_IC(:,3)/1000+0.2,num2str(xyz_IC(:,1)))
%      plot([x1 x2],[y1 y2],'ks-');
% 
% %      
%    plot([x1b x2b],[y1b y2b],'ks--');
% %  
% % %  % Plot data
% % % %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%   plot_locations = 0;
% %to generate RHIon AA
%  plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% %  r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
%  r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% hold all
% %  plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6)
% %  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6);
%  text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))
%  
% 
%   %to generate RHI on BB
% plot_vert_radar_plane2(fn,x1b,y1b,x2b,y2b,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)
% %  r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,2)/1000,xyz_CG(:,3)/1000);
%  r_IC = convert2D([x1b y1b x2b y2b],xyz_IC(:,2)/1000,xyz_IC(:,3)/1000);
% hold all
% %  plot(r_CG,xyz_CG(:,4)/1000,'ko','markerfacecolor','b','markersize',6);
% %  text(r_CG+0.1,xyz_CG(:,4)/1000+0.2,num2str(xyz_CG(:,1)))
% 
% plot(r_IC,xyz_IC(:,4)/1000,'ko','markerfacecolor','r','markersize',6);
%  text(r_IC+0.1,xyz_IC(:,4)/1000+0.2,num2str(xyz_IC(:,1)))


%% Plot vertical radar data- Over the sea storm
% This program is for generating the same line(AA') in different radar PPI , so
% tht vertical radar plane generated will be the same for each time period
% copy  and paste x'y'z point of the flash origination from excell sheeet
% per each 5 min period

%first using the plotter generate the PPI for the desired 'time' and 'angle' -we use angle 5(3.2 deg)here
 % User input
%    a = guidata(gcf);
%   fn = a.fileName;
%   hold all
%   xlim([5 45]);
%   ylim([5 45]);
% 
%  %**********IC/CG flash- time 21.45***************
% xyz_CG = [29662.4	22075.9	4755.4];
% xyz_IC = [26808.0	23119.0	7929.0];
% %cordinates for over the sea flash in '2011-07-22-Originationof flash project'
% %line AA'
%2nd cellat right
%  x1= 28.5955;
%  x2= 35.6067;
%  y1= 21.9888;
%  y2= 35.0225;


%1st cell(middle)
%  x1= 17.3559;
%  x2= 35.6844;
%  y1= 37.1037;
%  y2= 14.7406;
% 
%  
% hold all
%   %for CG
%   plot(xyz_CG(:,1)/1000,xyz_CG(:,2)/1000,'ko','markerfacecolor','b','markersize',6)
%   %for IC
%   plot(xyz_IC(:,1)/1000,xyz_IC(:,2)/1000,'ko','markerfacecolor','r','markersize',6)
%  plot([x1 x2],[y1 y2],'ks-')
%  
%  % Plot data
% %indicate wether you plot Ldar/CGlss,pbfa etc...,keep this 0 for now
%  plot_locations = 0;
% plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(1)
% delete(fID)
% 
%  %to generate RHI
% r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,1)/1000,xyz_IC(:,2)/1000);
% hold all
% plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_IC,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)


 %*********** CG/ IC flash- time 21.50-21.55 in the plotter***********
%   a = guidata(gcf);
%    fn = a.fileName;
%    plot_locations = 0;
   %CG
%    xyz_CG = [27586.2	24385.9	5972.0
%              27490.4	25179.2	6269.1
%              27109.6	25777.8	7240.6
%              ];
   %IC
%    xyz_IC = [25251	24660	7999];

 %line AA'
%  x1= 17.3559;
%  x2= 35.6844;
%  y1= 37.1037;
%  y2= 14.7406;
%  
% generate PPI
%    hold all
  % for CG
%   plot(xyz_CG(:,1)/1000,xyz_CG(:,2)/1000,'ko','markerfacecolor','b','markersize',6)
  % for IC
%   plot(xyz_IC(:,1)/1000,xyz_IC(:,2)/1000,'ko','markerfacecolor','r','markersize',6)
  %plot([x1 x2],[y1 y2],'ks-')


% % Plot data
% plot_locations =0;
% plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(1)
% delete(fID)


% r_CG = convert2D([x1 y1 x2 y2],xyz(:,1)/1000,xyz(:,2)/1000);
% r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,1)/1000,xyz_IC(:,2)/1000);
% hold all
% plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_IC,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)


 % ***********IC flash- time 22.01***************
%xyz_IC = [32464.7	26797.7	7647.9];
%xyz_CG= [NaN NaN NaN];% no CG flashes for this time 
% cordinates for over the flash in '2011-07-22-Originationof flash project'
%line AA'
%  x1= 17.3559;
%  x2= 35.6844;
%  y1= 37.1037;
%  y2= 14.7406;

%line BB' for flash8
% x1= 35.0000
% x2= 30.1601
% y1= 39.1993
% y2= 16.3523
% 
%   plot_locations = 0;
% 
%to generate PPI
%  hold all
% plot (xyz_CG(:,1),xyz_CG(:,2),'ro','makerfacecolor','r');
%  plot([x1 x2],[y1 y2],'ks-')
%  plot (xyz_IC(:,1)/1000,xyz_IC(:,2)/1000,'ro');
% 
% plot_vert_radar_plane2(fn,x1,y1,x2,y2,plot_locations)
% fID = gcf;
% radarTools(4)
% delete(fID)

%to generate RHI
%r_CG = convert2D([x1 y1 x2 y2],xyz_CG(:,1)/1000,xyz_CG(:,2)/1000);
% r_IC = convert2D([x1 y1 x2 y2],xyz_IC(:,1)/1000,xyz_IC(:,2)/1000);
% hold all
% %plot(r_CG,xyz_CG(:,3)/1000,'ko','markerfacecolor','b','markersize',6)
% plot(r_IC,xyz_IC(:,3)/1000,'ko','markerfacecolor','r','markersize',6)
%  
%% Field mill plot in xy (convert lat lon data to x, y , cordinates and plot in XPY plane)
% FM-07
% lat=28+38/60+32.5538/3600
% 
% lat =
% 
%    28.6424
% 
% lon=360-279 -15/60 -8.5739/3600
% 
% lon =
% 
%    80.7476
% 
% latlon2xy(lat,lon)
% 	x = -10266.9
% 	y = 11565.0
% 
% ans =
% 
%   -1.0267e+04
% 
% hold all
% plot(-10266.9/1000,11565.0/1000,'pk','markerfacecolor','b','markersize',6)
% text(-10266.9/1000,11565.0/1000,'FM-07')


 %FM-10
 %lat:28 37 26.8498
% FM22: lat: 28 30 24.9390	lon:279 18 23.0729
%% plot radar dBZ altitudes

% over the land, Storm 1
% ts = {'16:16:02', '16:19:28', '16:23:09', '16:26:49',  '16:31:00', ...
%     '16:34:41',  '16:38:51', '16:43:48', '16:47:43', '16:51:53', ...
%     '16:56:49', '17:01:14',  '17:05:54',  '17:10:20', '17:14:45'};
% 
% %t = nan(1,length(ts));
% for i =1:length(ts)
%     t(i) = hhmmss2sec(ts{i},0);
% end
% 
% 
% 
% 
% S1 = [ 1	6.66	6.558	6.525	NaN	    NaN
%        2	7.9	    7.783	7.298	7.065	NaN
%        3	8.213	7.94	7.9	    7.32	7.21
%        4	11.06	10.95	10.15	6.206	5.313
%        5	8.565	8.369	8.291	7.613	4.55
%        6	11.32	10.63	8.508	7.357	6.142
%        7	10.37	9.707	9.263	8.559	4.48
%        8	10.6	9.798	9.532	7.722	5.336
%        9	11.53	10.93	10.63	4.572	4.481
%        10	12.03	10.98	9.3874	4.687	NaN
%        11	11.76	11.37	6.961	4.529	NaN
%        12   11.8	11.22	6.396	5.39	NaN
%        13	11.66	10.18	6.395	4.6	    NaN
%        14	11.56	10.18	6.489	4.583	NaN
%        15	10.94	9.008	6.999	5.967	NaN];
%   figure
%   plot(t, [S1(:,2)], 'k*-');
%   hold all;
% 
%   plot(t, [S1(:,3)],'b*-');  
%   plot(t, [S1(:,4)],'g*-');  
%   plot(t, [S1(:,5)],'r*-');  
%   plot(t, [S1(:,6)], 'm*-');
%  % plot([59297 59297] ,[4 13], 'r--')  ;
%   plot(59297,7.922,'k*','MarkerFaceColor','k')  ;
%   text(59300,7.9, 'Origination point of the first IC flash');
%   title('Over the land : Thunderstorm 1');
%   xlabel('Time (UT)');
%   ylabel('Altitude (Km)');
%   legend('Cloud echo top ~ 10 dBZ',' 20 dBZ', '30 dBZ', '40 dBZ','50 dBZ');
%   set(gca,'xtick',58500:600:62100)
%   set(gca,'xticklabel',{'16:15','16:25','16:35','16:45','16:55','17:05','17:15'})
%   set(gca,'ytick',0:1:15);
%   
%  
%   
%  
% 
% % over the land, Storm 2 - Cell-1
% ts2 ={'16:56:49', '17:01:14','17:05:54','17:10:20','17:14:45',...
%     '17:19:10','17:23:35','17:28:01'}
% 
% for i=1:length(ts2)
%     t2(i)=hhmmss2sec(ts2{i},0);
% end
% S2= [	1   9.542	9.711	8.213	7.109	NaN
%         2   11.96	11.7	10.11	8.454	6.445
%         3   10.3	10.95	10.18	7.403	5.473
%         4   10.95	10.16	8.653	6.481	4.621
%         5   10.61	10.33	7.716	4.621	1.663
%         6   10.61	7.84	5.645	2.441	1.673
%         7   10.61	5.833	3.299	1.721	NaN
%         8   6.809	4.147	1.779	NaN	    NaN];
% 
% figure
%   plot(t2, [S2(:,2)], 'k*-');
%   hold all;
% 
%   plot(t2, [S2(:,3)],'b*-');  
%   plot(t2, [S2(:,4)],'g*-');  
%   plot(t2, [S2(:,5)],'r*-');  
%   plot(t2, [S2(:,6)], 'm*-');
%  % plot([59297 59297] ,[4 13], 'r--')  ;
%   plot(61383.4,10.1,'k*','MarkerFaceColor','k')  ;
%   text(61385, 9, 'Origination point of the first IC flash');
%   title('Over the land : Thunderstorm 1');
%   xlim([60600 63000]);
%   xlabel('Time (UT)');
%   ylabel('Altitude (Km)');
%   legend('Cloud echo top ~ 10 dBZ',' 20 dBZ', '30 dBZ', '40 dBZ','50 dBZ');
%   set(gca,'xtick',60600:600:63000)
%   set(gca,'xticklabel',{'16:50','17:00','17:10','17:20','17:30'})
%   %set(gca,'ytick',0:1:15);
% 
% 
% 
% % over the land, Storm 2 - Cell-2
% 
% ts3 ={'17:10:20','17:14:45','17:19:10','17:23:35','17:28:01','17:32:26',...
%     '17:36:51','17:41:16','17:45:41','17:50:37','17:55:02'};
% 
% for i=1:length(ts3)
%     t3(i)=hhmmss2sec(ts3{i},0);
% end
% 
% 
% 
% S3 =[1	7.027	7.027	7.248	6.242	NaN
%      2   7.806	7.344	6.676	4.777	4.042
%      3 	8.799	8.642	7.489	7.344	6.387
%      4  9.092	7.73	7.489	6.553	5.646
% 	 5  10.12	8.014	7.866	5.852	4.226
%      6	10.3	7.73	6.118	4.276	1.545
% 	 7  9.71	7.341	6.718	4.147	NaN
%      8	9.823	8.744	5.057	4.281	NaN
%   	 9  8.841	6.593	5.827	1.526	NaN
% 	10  7.537	6.469	5.936	5.059   NaN	
% 	11  7.863	6.801	5.21	4.869	NaN];
% 
% figure
% 
% plot(t3, [S3(:,2)], 'k*-');
%   hold all;
% 
%   plot(t3, [S3(:,3)],'b*-');  
%   plot(t3, [S3(:,4)],'g*-');  
%   plot(t3, [S3(:,5)],'r*-');  
%   plot(t3, [S3(:,6)],'m*-'); 
%   
%  % plot([59297 59297] ,[4 13], 'r--')  ;
%   plot(62793,7.95,'k*','MarkerFaceColor','k')  ;
%   text(62793, 7, 'Origination point of the first IC flash');
%   title('Over the land : Thunderstorm 2 (Cell 2)');
%   xlim([61500 64800]);
%   xlabel('Time (UT)');
%   ylabel('Altitude (Km)');
%   legend('Cloud echo top ~ 10 dBZ',' 20 dBZ', '30 dBZ', '40 dBZ','50 dBZ');
%   set(gca,'xtick',61500:600:64800)
%   set(gca,'xticklabel',{'17:05','17:15','17:25','17:35','17:45','17:55'})
% 
% 
% 
%   
% % over the sea storm -middle cell
% 
% ts4 ={'20:53:57','21:02:29','21:06:45','21:11:01','21:15:16','21:19:32',...
%     '21:23:48','21:28:03','21:36:35','21:40:50','21:45:06','21:49:22','21:53:38',...
%     '21:57:53','22:02:09','22:06:25'};
% 
% for i=1:length(ts4)
%     t4(i)=hhmmss2sec(ts4{i},0);
% %     t4_2(i)=hhmmss2sec(ts4{i},0)-hhmmss2sec(ts4{1},0);
% end
% 
% 
% 
% S4 =[   1 	5.881	4.728	2.972	NaN	    NaN
%         2   7.494	6.082	6.065	6.055   NaN	
%         3   7.61	7.61	7.559	3.899   NaN	
%         4   7.349	7.349	5.879	5.871   NaN	
%         5   9.126	7.344	7.366	5.886	4.686
%         6   9.251	9.157	7.466	7.429	5.947
%         7   9.326	9.317	7.52	6.008	4.818
%         8   11.46	9.253	5.977	2.997	2.045
%         9   9.79	9.716	7.893	7.811	6.283
%         10  12.49	12.38	9.985	9.951	6.457
%         11  12.61	12.37	10.19	8.192	6.542
%         12  12.53	12.2	9.986	6.611	4.137
%         13  12.66	10.3	8.286	5.346	3.346
%         14  12.61	9.984	6.731	4.193	1.37
%         15  12.45	10.29	5.359	3.359	2.347
%         16  12.36	10.29	5.483	3.435   NaN	
% ];

% figure
% plot(t4_2, [S4(:,2)], 'k*-');
%   hold all;
% 
%   plot(t4_2, [S4(:,3)],'b*-');  
%   plot(t4_2, [S4(:,4)],'g*-');  
%   plot(t4_2, [S4(:,5)],'r*-');  
%   plot(t4_2, [S4(:,6)],'m*-'); 
%   return

% figure
% plot(t4, [S4(:,2)], 'k*-');
%   hold all;
% 
%   plot(t4, [S4(:,3)],'b*-');  
%   plot(t4, [S4(:,4)],'g*-');  
%   plot(t4, [S4(:,5)],'r*-');  
%   plot(t4, [S4(:,6)],'m*-'); 
%   
%  % plot([59297 59297] ,[4 13], 'r--')  ;
%   plot(78304.5,4.75,'k*','MarkerFaceColor','k')  ;
%   text(78304, 4, 'Origination point of the first detectable CG flash (21:45,4.8');
%   title('Over the sea thunderstorm: middle Cell');
%   xlim([75000 79800]);
%   xlabel('Time (UT)');
%   ylabel('Altitude (Km)');
%   legend('Cloud echo top ~ 10 dBZ',' 20 dBZ', '30 dBZ', '40 dBZ','50 dBZ');
%   set(gca,'xtick',75000:600:79800)
%   set(gca,'xticklabel',{'20:50','21:00','21:10','21:20','21:40','21:50','22:00','22:10'})
% % 
% % % over the sea storm -Cell to the right
% 
% ts5 ={'21:36:35','21:40:50','21:45:06','21:49:22','21:53:38','21:57:53',...
%     '22:02:09','22:06:25','22:10:41','22:15:18','22:19:57'};
% 
% for i=1:length(ts5)
%     t5(i)=hhmmss2sec(ts5{i},0);
% end
% % 
% % 
% % 
% S5 =[  1 	12.28	6.191	2.145	NaN	    NaN
%         2   10.22	6.656	5.321	3.344	NaN
%         3   12.16	6.853	6.83	6.786	NaN
%         4   12.65	12.32	8.55	6.766	5.411
%         5   10.87	10.2	8.643	6.881	4.533
%         6   13.15	8.989	8.961	7.122	5.7
%         7   12.86	11.09	9.016	7.21	4.642
%         8   12.57	10.15	8.907	4.648	3.713
%         9   11.5	11.17	7.348	4.754	2.56
%         10  11.44	10.53	5.949	3.737	1.556
%         11  11.46	8.851	5.987	4.666	NaN
% 	
% ];
% % 
% % figure
% % 
% % plot(t5, [S5(:,2)], 'k*-');
% %   hold all;
% % 
% %   plot(t5, [S5(:,3)],'b*-');  
% %   plot(t5, [S5(:,4)],'g*-');  
% %   plot(t5, [S5(:,5)],'r*-');  
% %   plot(t5, [S5(:,6)],'m*-'); 
% %   
% %  % plot([59297 59297] ,[4 13], 'r--')  ;
% %   plot(79309.3,7.6,'k*','MarkerFaceColor','k')  ;
% %   text(79309, 7, 'Origination point of the first IC flash (22:01,7.6');
% %   title('Over the sea thunderstorm: Right Cell');
% %   %xlim([77700 80400]);
% %   xlabel('Time (UT)');
% %   ylabel('Altitude (Km)');
% %   legend('Cloud echo top ~ 10 dBZ',' 20 dBZ', '30 dBZ', '40 dBZ','50 dBZ');
% %  % set(gca,'xtick',777000:600:80400)
% %   set(gca,'xticklabel',{'21:35','21:45','21:55','22:05','22:15','22:25'})
% % 
% % %   
% % %   
% %   
% %   % over the sea storm -cell to the left
% % % 
% % ts6 ={'20:49:41','20:53:57','20:58:13','21:02:29','21:06:45','21:11:01',...
% %     '21:15:16','21:19:32','21:23:48','21:28:03','21:32:19','21:36:35','21:40:50',...
% %     '21:45:06','21:49:22','21:53:38','21:57:53','22:02:09'};
% % 
% % for i=1:length(ts6)
% %     t6(i)=hhmmss2sec(ts6{i},0);
% % end
% 
% %return
% % 
% % 
% % 
%  S6 =[  1   7.263	5.84	5.855	NaN     NaN    	
%         2   6.081	4.774	3.804	NaN	    NaN
%         3   6.108	4.817	3.827	NaN	    NaN
%         4   7.665	6.081	6.054	5.983   NaN	
%         5   7.592	6.188	6.108	6.027	NaN
%         6   9.671	7.889	6.376	6.242	5.033
%         7   9.914	8.021	7.823	5.055	5.142
%         8   10.04	9.76	7.987	5.163	3.255
%         9   9.992	8.063	8.062	6.52	1.755
%         10  9.671	7.889	6.376	6.242	5.033
%         11  10.33	6.754	6.727	4.287	1.056
%         12  10.44	8.319	8.253	4.261	1.067
%         13  10.48	8.419	5.403	4.45	NaN
%         14  10.36	8.253	4.579	1.118	NaN
%         15  8.885	8.32	4.609	NaN     NaN	
%         16  9.028	7.277	1.998	NaN     NaN	
%         17  12.12	9.683	4.892	NaN     NaN	
%         17  12.24	7.429	3.926	NaN	    NaN
% 
%  ];
% 
% figure
% 
% plot(t6, [S6(:,2)], 'k*-');
%   hold all;
% 
%     plot(t6, [S6(:,3)],'b*-');
%     plot(t6, [S6(:,4)],'g*-');
%     plot(t6, [S6(:,5)],'r*-');
%     plot(t6, [S6(:,6)],'m*-');
%   
%    % plot([59297 59297] ,[4 13], 'r--')  ;
%   %   plot(78304.5,4.75,'k*','MarkerFaceColor','k')  ;
%   %   text(78304, 4, 'Origination point of the first detectable CG flash (21:45,4.8');
%     title('Over the sea thunderstorm:  Cell to the left');
%   %   xlim([75000 79800]);
%     xlabel('Time (UT)');
%     ylabel('Altitude (Km)');
%     legend('Cloud echo top ~ 10 dBZ',' 20 dBZ', '30 dBZ', '40 dBZ','50 dBZ');
%    set(gca,'xtick',74700:600:79500)
%    xlim([74700 79500])
%   set(gca,'xticklabel',{'20:45','20:55','21:05','21:15','21:25','21:35','21:45','21:55','22:05'})
% 
% 
  
  
  %comparr over the sea cell as a bar chart
  
  
  
  
  
  
  
  













 
     function  addfm
 lats2 = [38 37 36 30];
 lats3 = [32.5538 26.8498 20.7371 24.9390];
 
 
 lons2 =[15 17 19 18];
 lons3=[8.5739 54.5539 34.0557 23.0729];
 
 lg = {'FM-07', 'FM-10', 'FM-11', 'FM-22'};

 
 for i = 1:length(lats2)
     lat=28+lats2(i)/60+lats3(i)/3600;
     lon=360-279 - lons2(i)/60 -lons3(i)/3600;
     
     [x,y] = latlon2xy(lat,lon);
     
     
     
     hold all
     plot(x/1000,y/1000,'pk','markerfacecolor','b','markersize',6)
     text(x/1000,y/1000,lg(i))
 end





%% Do not comment below this line
 function r = convert2D(vrbd,x,y)

x1 = vrbd(1);
x2 = vrbd(3);
y1 = vrbd(2);
y2 = vrbd(4);

% slope and intercept of the line (need later)
m = (y2 - y1)/(x2 - x1);
b = (x2*y1 - x1*y2)/(x2 - x1);

% Snap the location coordinates to line
if m == inf
    y = y;
    x = x1;
else
    x = (m*y+x-m*b)/(m*m + 1);
    y = (m*m*y+m*x+b)/(m*m + 1);
end

% Make them 2D
r = sqrt((x - x1).^2 + (y - y1).^2);

function setup_fig(rfn,xl,yl)


sns_x = [ -0.7524, 0.3372, -0.3149, -1.1658, -0.2701, -2.7640, -6.0091, 0.1825, -5.7394, -2.0637]*10;
sns_y = [  1.6555, 0.4446, -0.6838, -1.7020,  0.2631,  4.9254, -3.3983, -5.3008, 1.1923, 0.1569]*10;
sns_n = { 'K02', 'K14',    'K24',    'BCC',    'K17',    'EDW',    'STC',    'FLT',    'OVD',    'FFI'};


% setup

xlabel('East (km)')
ylabel('North (km)')

plot(sns_x,sns_y,'kp','markerfacecolor','k')
text(sns_x,sns_y,sns_n)
florida_map
box on

[rdn, rfn ,extn ] = fileparts(rfn);
title(sprintf('%s%s', rfn,extn),...
    'Interpreter','none')


xlim(xl)
ylim(yl)


daspect([1 1 1])



function ldar_plot_execute(handles)

if strcmp(handles.sen_set.autoBaseFolder,'on')
    handles = set_auto_b_folder(handles);
end

g = handles.g;
sen_set = handles.sen_set;

% set the base folder


if g.mm < 30
    ext=0;
else
    ext=30;
end

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);


linet_fn=sprintf('%s/LINET/%s/%s/%s/linet_%s%s%s_%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

% nldn_fn=sprintf('%s/LINET2/%s/%s/%s/linet_%s%s%s_%2.2i%2.2i.txt',...
%     sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
%     g.MM{:},g.DD{:},g.hh,ext);

nldn_fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

%Old PBFA fn
pbfaO_fn=sprintf('%s/PBFA_old/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    sen_set.loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);

nldn2_fn=sprintf('%s/NLDN2/%s/%s/NLDN2_%s%s%s.txt',...
        sen_set.loc_dir,g.YYYY{:},g.MM{:},g.YYYY{:},...
        g.MM{:},g.DD{:});
    
% Graph selection from the user
num=sen_set.ldar_graph_type;

% Do we need ldar time shift to be included?
if sen_set.ldar_tshiftOn==1
    x0=sen_set.x(sen_set.ldar_tshift_sn)-sen_set.x0;
    y0=sen_set.y(sen_set.ldar_tshift_sn)-sen_set.y0;
    z0=sen_set.z(sen_set.ldar_tshift_sn)-sen_set.z0;
else
    x0=0;
    y0=0;
    z0=0;
end

switch num
    case 1
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,sen_set,g)
    case 2
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,2,sen_set,g)
    case 3
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,3,sen_set,g)
    case 4
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,4,sen_set,g)
    case 5
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[1,0,0,0,0],sen_set,g);
    case 6
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[0,1,0,0,0],sen_set);
    case 7
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[0,0,1,0,0],sen_set,g);
    case 8
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[0,0,0,1,0],sen_set,g);
    case 9
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,1,[0,0,0,0,1],sen_set,g);
    case 10
        ldarColorTime2(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,0,[0,1,0,0,0],sen_set,g);
    case 11
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,5,sen_set,g)
    case 12
        LDAR_3D1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,nldn2_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,6,sen_set,g)
    case 13
        LDAR_Hist1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,sen_set)
    case 14
        LDAR_frequency(ldar_fn,linet_fn,pbfa_fn,nldn_fn,pbfaO_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),x0,y0,z0,sen_set)
        
    otherwise
        disp xc
        
end

function [radarFiles , radarFullFiles] = getRadarFiles(handles)

g = handles.g;
sen_set = handles.sen_set;

% Radar file folder
rdn = sprintf('%s/netCDF/%s/%s/%s/%s/%2.2i/',...
    sen_set.base_dir,sen_set.radarStations{sen_set.radarStationID},...
    g.YYYY{:},g.MM{:},g.DD{:},g.hh);

files = dir(rdn);

L = length(files);
times = nan(1,L-2);

radarFiles = {};
radarFullFiles = {};

if L > 2    
    for i = 3:L
        fn =  files(i).name;
        times(i-2) = str2double(fn(14:15))*3600+str2double(fn(16:17))*60+str2double(fn(18:19));
        radarFiles = [radarFiles; fn];
        radarFullFiles = [radarFullFiles; [rdn fn]];
    end
end

% if within the first 10 mins of the hour, let's load files from the
% previous hour
if g.mm < 10
    dNum1 = datenum(str2double(g.YYYY{:}),str2double(g.MM{:}),str2double(g.DD{:}),g.hh,0,0);
    
    % substract an hour
    dNum2 = dNum1 - 1/24;
    
    dvec = datevec(dNum2);
    
    % Radar file folder
    rdn = sprintf('%s/netCDF/%4.4i/%2.2i/%2.2i/%2.2i/',...
        sen_set.base_dir,dvec(1),dvec(2),dvec(3),dvec(4));
    
    files = dir(rdn);
    
    L = length(files);
    
    % We just need the last file of this folder
    if L > 2
        fn =  files(L).name;
        times =[times, str2double(fn(14:15))*3600+str2double(fn(16:17))*60+str2double(fn(18:19))];
        radarFiles = [radarFiles; fn];
        radarFullFiles = [radarFullFiles; [rdn fn]];
    end
end

% if within the first 10 mins of the hour, let's load files from the
% previous hour
if g.mm > 50
    dNum1 = datenum(str2double(g.YYYY{:}),str2double(g.MM{:}),str2double(g.DD{:}),g.hh,0,0);
    
    % add an hour
    dNum2 = dNum1 + 1/24;
    
    dvec = datevec(dNum2);
    
    % Radar file folder
    rdn = sprintf('%s/netCDF/%4.4i/%2.2i/%2.2i/%2.2i/',...
        sen_set.base_dir,dvec(1),dvec(2),dvec(3),dvec(4));
    
    files = dir(rdn);
    
    L = length(files);
    
    % We just need the first file of this folder
    if L > 2
        fn =  files(3).name;
        times =[times, str2double(fn(14:15))*3600+str2double(fn(16:17))*60+str2double(fn(18:19))];
        radarFiles = [radarFiles; fn];
        radarFullFiles = [radarFullFiles; [rdn fn]];
    end
end

if ~isempty(times)
    dts = abs(times - g.t1);
    [dts I] = sort(dts);
    radarFiles = radarFiles(I);
    radarFullFiles = radarFullFiles(I);
else
    radarFiles = '';
    radarFullFiles = '';
end
 