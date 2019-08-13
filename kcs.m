function kcs
% this will keep current settings

% Open the file
[fn,dir] = uigetfile('.mat','Choose the old setting file');

if fn ~= 0
    a = open([dir fn]);
else
    return
end

% verify it is the latest setting file suitable for this update
try 
    x = a.next_prev_dt;
catch
    msgbox('This is not the latest settings file - Please contact Sumedhe')
    return
end


try 
    x = a.cglssOn;
    msgbox('Your file is up to date. Nothing to do')
    return
catch
    % do nthing
end
    

a.cglssOn = 1;
save('sensor_setting.mat','-struct','a')

disp('Setting file successfully restored!')
