function generateDataTxt()
% This function is writtent to generate text files containing CH data.
% Modification History:
%   2018-02-03 Created by Sumedhe Karunarathne
%
% Requirements/how to use
%   1. Run Plotter 2
%   2. Choose date/time range you want to generate txt data
%   3. Select only one of the CH channel you want to generate data
%   4. Run this file
%
%   
% Output
%   1. A text file will be generated in the current working directory
%       The name of the text file follows the following format.
%       TXT_data_SNS_YYYYMMDD_sssss.sss.txt 
%           SNS - Sensor ID
%           YYYYMMDD - Date
%           sssss.sss is the starting time of the time range
%   2. The full file path will be displayed on the command window
%   3. If multiple CH files were selected seperate files will be created
%   for all selected data.
%
% WARNING: If this file name already exsist, that file will be replaced   
                  
%% Start the programm
%% Get plotter2 data
try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

settings = h.sen_set;
settings.sen_IDs{1};
g = h.g;


% what are the selected ch channels
output = plot_all6(g);
close(gcf);

sns = str2num(output.trigrd);

% Save data in a txt file
for i = 1:length(output.ch_data)
    if ~isempty(output.ch_data(i).t)
        % Sensor ID
        SNS = settings.sen_IDs{ceil(i/3)};
        
        % Channel ID
        if mod(i,3) == 0
            channel = 3;
        else
            channel = mod(i,3);
        end
        
        fileName = sprintf('TXT_data_%s_CH%1.1i_%s%s%s_%8.3f.txt', ...
            SNS,channel,g.YYYY{:},g.MM{:},g.DD{:},g.t1);
        
        fprintf('File name: %s\n',[pwd '\' fileName])
        
        
        fileID = fopen(fileName,'w');
        fprintf(fileID,'Date: %s-%s-%s\n',g.YYYY{:},g.MM{:},g.DD{:});
        fprintf(fileID,'Station: %s\n',SNS);
        fprintf(fileID,'Channel: %1.1i\n\n',channel);
        fprintf(fileID,'%6s %12s\n','Time (s)','E (V/m)');
        fprintf(fileID,'%6.7f %12.8f\n',[output.ch_data(i).t;output.ch_data(i).y]);
        fclose(fileID);
              
    end
end


    


