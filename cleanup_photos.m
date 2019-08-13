function cleanup_photos

% 2013-06-20 I normally shoot with Cannon T4i and RAW+JPEG mode. When I delete
%            some photos in camera in computer, only the JPEG was delted and
%            RAW file will stay in the computer taking lot of space. This
%            function was written to do the cleanup process. 
%            Selected folder and its subfolders will be checked.



% Get the directory
dn = uigetdir('E:\Sumedhe\photos\');

if dn(1) == 0
    disp('User cancelled the cleanup process.')
    return
end

% delte files in the current folder
delete_files(dn)
   

function delete_files(dn)

% get all the raw files in the folder
raws = dir([dn '\*.CR2']);


L = length(raws);

if L > 0
    fprintf('\nWorking on directory %s\n',dn)
    fprintf('========================================\n')
end

% delete counter
del_count = 0;

for i = 1:L
    
    fn = raws(i).name;
    
    if exist([dn '\' fn(1:end-3) 'jpg'],'file')
        fprintf('deleting on %s\n',raws(i).name);
        delete([dn '\' fn])
        del_count = del_count + 1;
    end    
end

if L > 0
    fprintf('\n%i of %i (%0.1f%%) files were delted!\n\n',...
        del_count,L,del_count/L*100)
end

% Get subdirectories
dnn = dir(dn);
isub = [dnn(:).isdir];

nameFolds = {dnn(isub).name};


for i = 3:length(nameFolds)
   delete_files([dn '\' nameFolds{i}])
end



