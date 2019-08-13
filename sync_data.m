function sync_data
% This function was written to sync (oneway) all location data (Ex. LDAR2, 
% PBFA, etc) from the SADAQ7 surver to Dr. Marshall's computer. 
%
%   Requirements: Before use, the remote folder's should 
%                 be mounted in the mac system.
%
%   Usage : sync_data
%   
%   History:
%       2014-02-20 Created by Sumedhe Karunarathne

clc
rDir(1).dn = '/Volumes/data/2010-07-01 -- 2010-08-19/data';
rDir(2).dn = '/Volumes/data/2011-07-07 -- 2011-07-16/data';
rDir(3).dn = '/Volumes/data/2011-07-17 -- 2011-08-04/data';
rDir(4).dn = '/Volumes/data/2011-08-05 -- 2011-08-16/data';
rDir(5).dn = '/Volumes/data/newRaid2/data';


Ldir = '/Users/thomas/Documents/KSC_2011_data/LocData';

dType = {'cglss','ldar2','LINET','LINET2','NLDN2','PBFA','PBFA2','PBFA_old'};


for i = 1:5
    
    for j = 1:length(dType)
        % sync LDAR
        cmd = sprintf('rsync -rv --update "%s/%s" "%s/"',...
            rDir(i).dn,dType{j},Ldir);
        
        system(cmd);
    end    
end

