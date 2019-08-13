function update_potter


URL = 'http://www.phy.olemiss.edu/~sumedhe/Plotter2.zip';
wd = pwd;
filename = [wd '/Plotter2.zip'];

wbh = waitbar(0.1,'Closing Plotter2...','name','Updating plotter2');

%% Close plotter2
try
    h=guidata(findall(0,'Tag','plotter2'));
    delete(h.plotter2)    
end

%% Download updates
waitbar(0.3,wbh,'Downloading updates ...','name','Updating plotter2');
[dispose status] = urlwrite(URL,filename);

if ~status
    delete(wbh)
    errordlg('Download failed. Please check the internet connection and retry.');
    return
end

%% Backup current folder.
waitbar(0.4,wbh,'Backing up current folder ...','name','Updating plotter2');
inds = find(wd == '/');

if isempty(inds)
    inds = find(wd == '\');
end

%% Take user aprove before backing up data
back_fol = [wd(1:inds(end)) 'Plotter2_backup_' datestr(now,'yyyymmdd-HHMMSS')];

msg = sprintf('Current Plotter2 files will be backed to the folder \n\n%s\n\nIs it ok?',back_fol);

button = questdlg(msg,'Backup','OK','Not cool!','Not cool!');

if ~strcmp(button,'OK')
    delete(wbh)
    errordlg('User did not aprove the update procedure. No harm done to your exsisting files')
    return
end
waitbar(0.5,wbh,'Backing up current folder ...','name','Updating plotter2');

success = movefile(wd,back_fol,'f');

if ~success
    msg = 'Auto backup failed. Please manually backup current Plotter2 files and click OK';
    button = questdlg(msg,'Backup','OK','Exit','OK');
    
    if ~strcmp(button,'OK')
        delete(wbh)
        errordlg('Update procedure failed. Please make sure your Plotter2 folder is ok.')
        return
    end        
end

%% Updating new files
waitbar(0.6,wbh,'Backing up current folder ...','name','Updating plotter2');

waitbar(0.7,wbh,'Updating Plotter2 ...','name','Updating plotter2');
% Put new files back to Plotter2 folder
unzip([back_fol '/Plotter2.zip'],wd(1:inds(end)));

% remove downloaded zip file
zipfile = [back_fol '/Plotter2.zip'];
delete(zipfile)



%% Restore last gui data
waitbar(0.8,wbh,'Restore Plotter2 settings ...','name','Updating plotter2');
wd = pwd;
g1 = open([wd '/plotter_last_gui_data.mat']);
g2 = open([back_fol '/plotter_last_gui_data.mat']);
newFields = fieldnames(g1);

for i = 1:length(newFields)
    if isfield(g2,newFields{i})
       g1.(newFields{i}) = g2.(newFields{i});
    end
end

save([wd '/plotter_last_gui_data.mat'],'-Struct','g1')

%% Restoring Setting files
waitbar(0.9,wbh,'Restore Plotter2 settings ...','name','Updating plotter2');
wd = pwd;
sen_set1 = open([wd '/sensor_setting.mat']);
sen_set2 = open([back_fol '/sensor_setting.mat']);
newFields = fieldnames(sen_set1);

for i = 1:length(newFields)
    if isfield(sen_set2,newFields{i})
       sen_set1.(newFields{i}) = sen_set2.(newFields{i});
    end
end

save([wd '/sensor_setting.mat'],'-Struct','sen_set1')

%% Reopenining plotter2
waitbar(1.0,wbh,'Done!.','name','Updating plotter2');
plotter2;
delete(wbh)
uiwait(msgbox('Plotter2 succesfully updated! Close this to contine.','Plotter2 Updator'));













