function copyDataAtoB

%% User Inputs
% Copy to base dir
dn_bd = 'H:\data\';

% date and time
dates =      {'20110722','20110725', '20110726', '20110731', '20110812', '20110813', '20110814'}'
startTimes = [190000,    224000,     185500,     193500,      191500,    213000,     212500 ]'
endTimes =   [235500,    230000,     200500,     194000,      204500,    222000,     231000]'

dates = {'20100711', '20100712', '20100801', '20100815', '20100817'}'
startTimes = [190000, 180000, 180000, 004500, 155000]'
endTimes = [210000, 190000, 190000, 024500, 164000]'


dates = {'20110807'};
startTimes = [195000];
endTimes = [205000];

%% Start the program

for k = 1:length(startTimes)
    
    dateStr = dates{k};
    startT = startTimes(k);
    endT = endTimes(k);
    
    date = str2double(dateStr);
    
    prefix = 'C:';
    
    if date <= 20110627
        bd = [prefix '\data\2010-07-01 -- 2010-08-19\data\'];
    elseif date <= 20110716
        bd = [prefix '\data\2011-07-07 -- 2011-07-16\data\'];
    elseif date <= 20110804
        bd = [prefix '\data\2011-07-17 -- 2011-08-04\data\'];
    elseif date <= 20110816
        bd = [prefix '\data\2011-08-05 -- 2011-08-16\data\'];
    else
        % Do nothing
        return
    end
    
    % copy from dir
    cfdir = sprintf('%s%s/%s/%s/',bd,dateStr(1:4),dateStr(5:6),dateStr(7:8));
    
    % copy to dir
    ctdir = sprintf('%s%s/%s/%s/',dn_bd,dateStr(1:4),dateStr(5:6),dateStr(7:8));
    
    if ~exist(ctdir,'dir')
        mkdir(ctdir);
    end
    
    % get all the file names
    files = dir(cfdir);
    L = length(files);
    
    for i = 3:L
        
        fn = files(i).name;
        
        try
            t = str2double(fn(14:19));
            
            if t >= startT && t <=endT
            
            source = [cfdir fn];
            destin = [ctdir fn];
            
            % if destination exist compare size and copy if not equal the
            % file size
            if exist(destin,'file')
                % get the info of the file
                info1 = dir(source);
                info2 = dir(destin);
                
                % Small file size, copy.
                if info1.bytes > info2.bytes
                    disp(['Copying...  ' fn])
                    copyfile(source,destin)
                else
                    disp(['Skipped. Identical  ' fn])
                end
            else
                disp(['Copying...  ' fn])
                copyfile(source,destin)
            end
            end
        
        catch
            disp(['somthing wrong with file ' fn])
        end
    end    
end

    
    