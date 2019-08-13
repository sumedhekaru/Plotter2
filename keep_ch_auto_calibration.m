function keep_ch_auto_calibration
% the 'auto_calibration.m' program will calibrate lp2 aginst field 
% mill data. However, not all the calibration plots are good. This
% program will help to keep good ones.
% Press 'y' to accept
% Press 'n' to reject
% Press 'e' to exit the program

% Folder name contain the files
%fn1 = 'C:\Users\sumedhe\Desktop\ch3_auto\';
fn1 = uigetdir('C:\Users\daqop\Desktop\Sumedhe\Calibration\Auto\K24\20110814_single_stroke\','Select the folder');

if fn1 == 0
    disp('User did''t select a folder')
    return
end
fn1 = [fn1 '\'];    

% Folder name to move files
fn2 = [fn1 'selected\'];
if ~exist(fn2,'dir')
    mkdir(fn2); 
end

% File to save calibration data
cal_fn = 'selected_ch1_cal_data.txt';
fID = fopen([fn2 cal_fn],'a+');


close all
list = dir(fn1);

[L ~] = size(list);

data = [];
%Find the file with all the data
for i = 3:L;    
    item = list(i);
    ext =  item.name(end-2:end);
    
    if strcmp(ext,'txt')
        dfID = fopen([fn1 item.name],'r');
        data=textscan(dfID,'%s %f %f %f');
        fclose(dfID);
        break
    end
end

if isempty(data)
    disp('Calibration data txt file couldn'' locate')
    return
end

% Sort raws according to peaks
[ix,ix] = sortrows(data{4},-1);

for i=1:4
    data{i} = data{i}(ix,:);
end

L = length(data{1});

name = data{1};
G    = data{2};

moved = 0;

for i = 1:L 
    
    try   
        fg = open([fn1 name{i}]);
        ylabel(sprintf('Moved %i',moved))
        
        gain = G(i);
        
        cc = 'o';
        
        while ~ismember(cc,{'e','y','n'})
            figure(fg)
            cc = get(fg,'CurrentCharacter');
            pause(0.1)
        end
        
        switch cc
            case 'e'
                close(fg)                
                fclose(fID);
                disp('User terminated the operation')
                return
            case 'y'
                moved = moved + 1;
                fprintf('moved %i\n',moved)
                movefile([fn1 name{i}],[fn2 name{i}])
                fprintf(fID,'%s\t%0.3f\n',name{i},gain);               
            case 'n'
                disp('skipped...')
        end        
        close(fg)
    catch
        fprintf('%s not found\n',name{i})
    end

end

fclose(fID);


