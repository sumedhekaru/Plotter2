function KSC_file_converting

time = '2330';

% Sending folder
%sfd = '\\SADAQ7\data\2010-07-01 -- 2010-08-19\data\FM\2010\08\01\';
sfd = 'C:\Users\sumedhe\Desktop\FMdata-conversion\15\FM';
%folder
dn = 'C:\Users\sumedhe\Desktop\FMdata-conversion\15\';



list = dir(dn);

L = length(list);


for i = 3:L;
    disp = list(i);
    fn = disp.name;

    ext = fn(1,end-2:end);
    sn  = str2double(fn(1,4:5));
    
    if strcmp(ext,'RAW')
        
            oldname = [dn fn];
            newname = [sfd fn(1,1:end-4) '_' time '.txt'];

            if exist(newname,'file')
                fprintf('File already exists!\n%s\n',newname)
            else
                movefile(oldname,newname);
            end      
    end
end
