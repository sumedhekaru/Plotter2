function PBFA_cleanup_backup_files
% Each time when you find a PBFA point, its going to make a backup file.
% This will build up over the time and running this file will cleanup those
% backup files.

% User inputs
bf = 'C:\data\2011-08-05 -- 2011-08-16\data\PBFA_old\';


% Program
for year = [2010 2011]
    for month = 1:12
        for day = 1:31
            bf2 = sprintf('%s/%4.4i/%2.2i/%2.2i/',bf,year,month,day);
            files = dir(bf2);
            
            L = length(files);
            
            if L > 0;
                for i = 1:L
                    answer = strfind(files(i).name,'backup');                    
                    if ~isempty(answer)
                        fprintf('Deleting %s ...\n',files(i).name);
                        delete([bf2 files(i).name])
                    end
                end
            end
        end
    end
end