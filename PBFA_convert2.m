function PBFA_convert2
% There were 3 files generated for automatic PBFA for negative begining pulses,
% positive IBPs and return strokes. This program will take data in those files and
% converted in to the dirrectory structure that plotter 2 can understand.


% negative IBP file name for PBFA begining pulses
fn1 = 'C:\Users\Sumedhe\Desktop\PBFA_auto_TOAM5\PBFA_auto_test_20110814_2100-2200_IC.txt';

% file name positive PBFA pulses
fn2 = 'C:\Users\Sumedhe\Desktop\PBFA_auto_TOAM5\PBFA_auto_test_20110814_2100-2200.txt';

% file name for PBFA cg positions
fn3 = 'C:\Users\Sumedhe\Desktop\PBFA_auto_TOAM5\PBFA_CG_data_20110814_2100-2200.txt';






% Output file folder
out_dir='C:/Users/sumedhe/Desktop/PBFA_auto_TOAM5/PBFA2/2011/08/14/';

if exist(out_dir,'dir')==0    
    mkdir(out_dir)
end

%% Loading PBFA IBP data
% % Open and read all data in the first file

[data1 snsStr1 peaks1] = read_data(fn1,1);
[data2 snsStr2 peaks2] = read_data(fn2,2);
[data3 snsStr3 peaks3] = read_data(fn3,3);

size(data3)

data = [data1; data2; data3];
peaks = [peaks1; peaks2; peaks3];
snsStr = [snsStr1; snsStr2; snsStr3];

[data indx] = sortrows(data,2);
peaks = peaks(indx,:);
snsStr = snsStr(indx);


% Let's put into files
t1 = -1800;
t2 = 0;

for k=1:48
    
   
    t1 = t1+1800;
    t2 = t2+1800;
    
    lol = sum(data(:,2) < t1) + 1;
    ul  = sum(data(:,2) < t2);
    
    if ul > lol
        
        str=sec2hhmmss(t1);
        fn = [out_dir 'pbfa_20110814_' str(1:2) str(4:5) '.txt'];
        
        if exist(fn,'file')
            pbfaD = pbfaExtract(fn,t1,t2,1000000,0,0,0);
            [index ~] = size(pbfaD);
            index = index + 1;
            fID = fopen(fn,'a+');
        else
            index = 1;
            fID = fopen(fn,'a+');
            fprintf(fID,'\n\n');
        end
  
        
        
        for i=lol:ul
            
            % Is it a CG point? Just add it (filtereing has done by
            % PBFA_CG_finder.
            

            if data(i,8) < 3
                % Lets see Z coordinate is smaller than 2km
                if data(i,5) > 2000
                    fprintf('%i\t%0.8f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.2f\t%0.0f\t%0.0f\t%0.0f\t%0.1f\n',...
                        index,data(i,2:7),snsStr{i},data(i,9:13))
                    fprintf(fID,'%i\t%0.8f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.2f\t%0.0f\t%0.0f\t%0.0f\t%0.1f\n',...
                        index,data(i,2:7),snsStr{i},data(i,9:13));
                    index = index + 1;
                end
                
                
                
            elseif data(i,8) == 3
                fprintf('%i\t%0.8f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.2f\t%0.0f\t%0.0f\t%0.0f\t%0.1f\n',...
                        index,data(i,2:7),snsStr{i},data(i,9:13))
                fprintf(fID,'%i\t%0.8f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\t%0.2f\t%0.0f\t%0.0f\t%0.0f\t%0.1f\n',...
                    index,data(i,2:7),snsStr{i},data(i,9:13));
                index = index + 1;
            end
                
        end
        
        fclose(fID);
    end
end



function [data snsStr peaks] = read_data(fn,type)

% Discription 
% 1 - Index     2 - time    3 -5 - xyz
% 6 - N         7 - Ip      8 - SNS     9 - 12 dt,dx,dy,dz  13 - kisquired 
if type == 1 || type == 2
    % Reading files from PBFA_auto4 output   
    fid = fopen(fn);
    data = textscan(fid,'%f %f %f %f %f %f %f %s %f %f %f %f %f','HeaderLines',2);
    snsStr = data{8};
    data{8} = data{1} - data{1} + type;
    fclose(fid);
        
    fid = fopen([fn(1:end-4) '_peaks.txt']);
    peaks = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    fclose(fid);
       
elseif type == 3
    fid = fopen(fn);
    data = textscan(fid,'%f %f %f %f %f %f %f %s %f %f %f %f %f %*f','HeaderLines',2);
    snsStr = data{8};
    data{8} = data{1} - data{1} + type;
    data{13} = data{9}; % Column 13 should be ki squired
    data{9} = NaN(size(data{9}));
    data{10} = data{9};
    data{11} = data{9};
    data{12} = data{9};
    fclose(fid);
    
    fid = fopen([fn(1:end-4) '_peaks.txt']);
    peaks = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    fclose(fid);
end

data = cell2mat(data);
peaks = cell2mat(peaks);




