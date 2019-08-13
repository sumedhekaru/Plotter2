function convert_linet

wbh = waitbar(0,'Please Wait');

% This function convert the linet files to normal files we were using

% Input file name
fn = 'C:\Users\sumedhe\Downloads\110801.txt';

% Output file folder
out_dir='C:/Users/Sumedhe/Desktop/data/LINET-new/2010/08/01';
if exist(out_dir,'dir')==0    
    mkdir(out_dir)
end


waitbar(0,wbh,'Reading File');

% Open and read all data in the file
fid = fopen(fn);
data=textscan(fid,'%f %f:%f:%f %f %f %f %f %f %f','HeaderLines',0);
% Data contains the following
% 1 - date
% 2 - Hours
% 3 - minutes
% 4 - seconds
% 5 - Latitude
% 6 - Logitude
% 7 - Z
fclose(fid);

% Total number of points
tot1=length(data{1});

% find occuring time in seconds - store at 1st column
data1=data{2}*3600+data{3}*60+data{4};

% Convert data in to mat format
data=cell2mat(data);

%Date 
date = data(1,1);


%Lower limit
ll=0;
lln=1;
wbn=0;
%Write File names
for i = 0:47;   
    
    hh = floor(i/2);
    mm = rem(i,2)*30;
    
    wfn=sprintf('%s/linet_%.0f_%2.2i%2.2i.txt',out_dir,date,hh,mm)
    fid1 = fopen(wfn,'w');
    
    ul=ll+1800;    
    uln = sum(data1<=ul);
    
    % Writing the header 
    fprintf(fid1,'Linet\nh\tm\ts\tms\tmiks\tLon\tLat\tZ\tType\tIp\tAcc\n ');
    
    for j=lln:uln
        sec = floor(data(j,4));
        ms  = floor((data(j,4)-sec)*1000);
        us  = (data(j,4)-sec-ms/1000)*1e6;
        fprintf(fid1,'\n%2.2i\t%2.2i\t%2.2i\t%3.3i\t%7.3f\t%8.4f\t%7.4f\t%7.1f\t%i\t%6.1f\t%6.3f',...
            data(j,2),data(j,3),sec,ms,us,data(j,6),data(j,5),data(j,7),data(j,8),data(j,9),data(j,10));
        
    end
    
    fclose(fid1);
    
    if uln < lln
        delete(wfn)
    else
        lln=uln;
    end
    ll=ul;
    
    
    
    str=sprintf('Now writing...\nlinet-%.0f-%2.2i%2.2i.txt',date,hh,mm);
        
    try
        waitbar(i/48,wbh,str)
    catch
        fprintf('\nOperation Terminated by User\nV1.0 - Sumedhe \n')
        return
    end  
end


delete(wbh)




