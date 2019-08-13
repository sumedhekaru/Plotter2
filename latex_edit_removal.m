function latex_edit_removal
fclose all;
clc

filename = 'C:\Users\Sumedhe\Desktop\NBP-paper-JGR\Second_submission\NBP_manuscript_R1_V04-no-tracking.tex';


fID1=fopen(filename, 'r');

%[~, fileName0, ~] = fileparts(filename);
fID2=fopen([filename(1:end-4) '-MLAB-untracked.tex'],'w');

tline = fgets(fID1);
%fprintf(fID2,'%s',tline);

k1 = strfind(tline, '\add{');

while ischar(tline)
    
    while ~isempty(k1)
        
        start = k1(1);
        k2 = strfind(tline, '}');
        k3 = strfind(tline, '{');
        
        k3(k3<start) = [];
        k2(k2<start) = [];
        
        k2 = [k2;zeros(size(k2))+1];
        k3 = [k3;zeros(size(k3))-1];
        
        ks = sortrows([k2'; k3'],1);
        
        [L1 L2] = size(ks);
        
        value = -1;
        
        for i=2:L1
            value = value + ks(i,2);
            
            if value == 0
                endind = ks(i,1);
                break
            end
        end
        tline
        
        
        tline(start:start+4) = [];
        tline(endind-5)      = [];
        
        k1 = strfind(tline, '\add{');
        
    end
    
    k1 = strfind(tline, '\remove{');
    
    while ~isempty(k1)
        
        start = k1(1);
        k2 = strfind(tline, '}');
        k3 = strfind(tline, '{');
        
        k3(k3<start) = [];
        k2(k2<start) = [];
        
        k2 = [k2;zeros(size(k2))+1];
        k3 = [k3;zeros(size(k3))-1];
        
        ks = sortrows([k2'; k3'],1);
        
        [L1 L2] = size(ks);
        
        value = -1;
        
        for i=2:L1
            value = value + ks(i,2);
            
            if value == 0
                endind = ks(i,1);
                break
            end
        end
        
        tline(start:endind) = [];
        
        k1 = strfind(tline, '\remove{');
        
    end
    
    
    % write that line to the new file
    fprintf(fID2,'%s',tline);
    
    
    tline = fgets(fID1);
    k1 = strfind(tline, '\add{');
    
end


fclose(fID1);
fclose(fID2);
