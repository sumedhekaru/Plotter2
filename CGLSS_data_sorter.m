function CGLSS_data_sorter

sen_set = open('sensor_setting.mat');

% Plot range
frm = '20110801';
too = '20110815';

% Date numbers
frmdn = datenum(frm,'yyyymmdd');
toodn = datenum(too,'yyyymmdd');

L = toodn-frmdn

rc = 30000;
t1 = 46800;
t2 = 7200;

d = cell(L,2);
xlab = cell(L,1);
cnt = 0;

for i=frmdn:toodn;
    fprintf('%s\n',datestr(i))
    cnt = cnt + 1;
    
    mm  = datestr(i,'mm');
    dd = datestr(i,'dd');
    fn1 = sprintf('%s/cglss/2011/%s/KSCCGLSS2011%s%s.dat',sen_set.base_dir,mm,mm,dd);
    
    data=CGLSS_extract(fn1,43200,86400,rc,0,0,0);
    tot1 = length(data);
    
    data=CGLSS_extract(fn1,t1,86400,rc,0,0,0);
    part1 = length(data);
    
    xlab{cnt,1} = datestr(i,'mmdd') ;
    
    
    mm  = datestr(i+1,'mm');
    dd = datestr(i+1,'dd');
    fn1 = sprintf('%s/cglss/2011/%s/KSCCGLSS2011%s%s.dat',sen_set.base_dir,mm,mm,dd);
    
    data=CGLSS_extract(fn1,0,43200,rc,0,0,0);
    tot2 = length(data);
    
    data=CGLSS_extract(fn1,0,t2,rc,0,0,0);
    part2 = length(data);
    
      
    %d{cnt,1} = cnt;    
    d{cnt,1} = tot1+tot2;
    d{cnt,2} = part1+part2;    
end



figure
d = cell2mat(d);

daylyAvg = mean(d(:,1))
TimeRangeAvg = mean(d(:,2))


bar(d,1.8,'hist');
xlim([0,L+2])
set(gca,'XTick',1:1:cnt)
set(gca,'XTickLabel',xlab)
rotateticklabel(gca,90);
title('CGLSS2 detection within 30km')
%xlabel('Date (mmdd)')
ylabel('Number of CGLSS2 locations')
legend('All','1300-0200UT')
