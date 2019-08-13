function keep_ch_auto_calibration
% The file is based on Keep_ch_auto_calibtation.
% For NBP data generated using NBP_finder3
% Press 'y' to accept
% Press 'n' to reject
% Press 'e' to exit the program

% Folder name contain the files
%fn1 = 'C:\Users\sumedhe\Desktop\ch3_auto\';
% fn1 = uigetdir('C:\Users\daqop\Desktop\Sumedhe\Calibration\Auto\K24\20110814_single_stroke\','Select the folder');
fn1 = 'C:\Users\Sumedhe\Desktop\NBP\20110814\20110814\FLT';


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
cal_fn = 'selected_NBP_data.txt';
fID = fopen([fn2 cal_fn],'a+');


close all
list = dir(fn1);

[L ~] = size(list);

%Find the file with all the data

dfID = fopen([fn1 '20110814-NBP_info.txt'],'r');
data=textscan(dfID,'%f %f %f %f');
fclose(dfID);


L = length(data{1});

name = data{1};
times = data{2};
FWHM = data{3};
SNR = data{4};

moved = 0;

i = 0;
while i <= L
    i = i + 1;
%for i = 1:L 
    
   try 
        fg = figure;
        cla
        [I,map] = imread(sprintf('%s%4.4i.png',fn1,i),'png');
        imshow(I,map);
        
        ylabel(sprintf('Moved %i',moved))
        
        time = times(i);
        
        cc = 'o';
        
        while ~ismember(cc,{'e','y','n','b'})
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
                sprintf('%s/%4.4f.png',fn1,name(i))
                copyfile(sprintf('%s/%4.4i.png',fn1,name(i)),sprintf('%s/%4.4i.png',fn2,name(i)))
                fprintf(fID,'%4.4i\t%13.6f\t%5.1f\t%5.1f\n',...
                    name(i),time,FWHM(i),SNR(i));   
            case 'b'
                % Go back one
                i = i - 2;
            case 'n'
                disp('skipped...')
        end        
        close(fg)
    catch
        fprintf(sprintf('%s/%4.4f.png not found.',fn2,name(i)))
    end

end

fclose(fID);


