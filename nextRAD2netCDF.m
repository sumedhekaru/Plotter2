function nextRAD2netCDF

clc
%% User inputs
%bf = 'C:\data\2010-07-01 -- 2010-08-19\data\';
%bf = 'C:\data\2010-07-01 -- 2010-08-19\data\';
bf = 'C:\data\2011-08-05 -- 2011-08-16\data\';
%bf = 'C:\data\2011-07-17 -- 2011-08-04\data\';

date = '2011\08\08\';
station = 'MEL';

%% Start the program


wbh = waitbar(0,'Converting NEXRAD to netCDF');

bit = 0;
for i = 0:24 
    dn1 = sprintf('%sRADAR/%s/%s%2.2i/',bf,station,date,i);
    dn2 = sprintf('%snetCDF/%s/%s%2.2i/',bf,station,date,i);
    
    if ~exist(dn2,'dir')
        mkdir(dn2)
    end
    
    files = dir(dn1);
    
    L = length(files);
    
    
    for j = 3:L
        bit = bit + (1/24)*(2+j)/L;
        try
            waitbar(bit,wbh)
        catch
            return
        end
        
        fn = files(j).name;
        
        cmd = sprintf('"C:/Program Files/Java/jre7/bin/java" -classpath "C:/data/toolsUI-4.3.jar" ucar.nc2.FileWriter -in "%s%s" -out "%s%s.netcdf"',...
            dn1,fn,dn2,fn)
        system(cmd)
        %return
    end
end