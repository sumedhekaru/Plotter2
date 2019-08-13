function convert_linet2

wbh = waitbar(0,'Please Wait');

% This function convert the linet files to normal files we were using
% Note that this is the second version they send some time summer 2011.
% Make sure to give 3 inputs, input file name, output file name and the
% date.

% Input file name
fn = 'C:\Users\sumedhe\Downloads\TOM.txt';

% Output file folder
out_dir='C:/Users/Sumedhe/Desktop/data/LINET2/2010/07/28';

% date
date = 2010728;


if exist(out_dir,'dir')==0    
    mkdir(out_dir)
end


waitbar(0,wbh,'Reading File');

% Open and read all data in the file
fid = fopen(fn);
data=textscan(fid,'%f %f %f %f %f %f %f %f %f %f','HeaderLines',0);
% Data contains the following
% 1 - Hours
% 2 - minutes
% 3 - seconds
% 4 - mili seconds
% 5 - micro seconds
% 6 - Lon
% 7 - Lat
% 8 - Z
% 9 - unkown at the time this is writing
% 10 - unkown
fclose(fid);

% Total number of points
tot1=length(data{1});

% find occuring time in seconds 
data1=data{1}*3600+data{2}*60+data{3}+data{4}/1e3+data{5}/1e6;

% Convert data in to mat format
data=cell2mat(data);




%Lower limit
ll=0;
lln=1;
wbn=0;
%Write File names
for i = 0:47;   
    
    hh = floor(i/2);
    mm = rem(i,2)*30;
    
    wfn=sprintf('%s/linet_%.0f_%2.2i%2.2i.txt',out_dir,date,hh,mm);
    fid1 = fopen(wfn,'w');
    
    ul=ll+1800;    
    uln = sum(data1<=ul);
    
    % Writing the header 
    fprintf(fid1,'Linet\nh\tm\ts\tms\tmiks\tLon\tLat\tZ\n ');
    
    for j=lln:uln
        sec = data(j,3);
        ms  = data(j,4);
        us  = data(j,5);
        fprintf(fid1,'\n%2.2i\t%2.2i\t%2.2i\t%3.3i\t%.3f\t%.4f\t%.4f\t%.1f',...
            data(j,1),data(j,2),sec,ms,us,data(j,7),data(j,6),data(j,9));
        
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




