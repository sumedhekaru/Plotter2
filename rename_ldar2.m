function rename_ldar2

% Get directory name
dn = uigetdir;

if dn == 0
    disp('Folder not selected')
    return
end


l = ls(dn);

for i=3:length(l)
    % initial file name
    ifn = sprintf('%s/%s',dn,l(i,:));
    % final file name
    ffn = sprintf('%s/ldar2_%s',dn,l(i,:));
    
    movefile(ifn,ffn)
end

fprintf('%s Done\n', dn)