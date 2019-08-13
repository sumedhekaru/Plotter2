function pulse_interval_histogram

%% User inputs
files = {'PulseIntervalData-20100712.txt'
         'PulseIntervalData-20100715.txt'         
         'PulseIntervalData-20100801.txt'
         'PulseIntervalData-20100815.txt'
         'PulseIntervalData-20100817.txt'
         'PulseIntervalData-20100818.txt'
         'PulseIntervalData-20110722.txt'
         'PulseIntervalData-20110814.txt'
    };

% Inter-flash time treshold in seconds (this number will automatically figure out the
% data belongs to next flash)
inter_flash_threshold = 5e-3;


%% Program
L = length(files);

dTs = [];

for i = 1:L
    
    fID = fopen(files{i});
    data = textscan(fID,'%s %f %f');
    fclose(fID);
    
    t = data{2};
    A = data{3};
    
    % Try to open the version two files   
    fID = fopen(['PulseIntervalData2-' files{i}(19:end)]);
    
    if fID ~= -1
        data = textscan(fID,'%s %f %f');
        fclose(fID);
        
        t = [t; data{2}];
        A = [A; data{3}];
    end
    
    % sort data
    [t inds] = sort(t);
    A = A(inds);
    
    % make groups
    dt = t(2:end)-t(1:end-1);   
    inds = find(dt > inter_flash_threshold);
    
    if isempty(inds)
        inds = length(dt);
    end
    % Number of groups
    L1 = length(inds);
    
    stInd = 1;
    enInd = inds(1);
        
    for i = 1:L1-1
        tsub = t(stInd:enInd);
        Asub = A(stInd:enInd);
        
        [mm, mind] = max(Asub);
        
        dTs = [dTs ;(tsub(mind)-tsub)];
        
        %max((tsub(mind)-tsub))
        
        stInd = inds(i)+1;
        enInd = inds(i+1);        
        
    end
end

figure
[N x] = hist(dTs*1000,[-5:0.5:10]);
ylabel('Number of pulses')
xlabel('dT (ms)')

hist(dTs*1000,30)

% print data
for i = 1:length(N)
    fprintf('%4.1f\t%4.i\n',x(i),N(i))
end