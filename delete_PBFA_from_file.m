function [fn1, backfile1, done1 ,fn2, backfile2, done2] = delete_PBFA_from_file(t,dataType,delType)
% This function is writtent to delete a PBFA point from a file.

% initialize return values
backfile1 = '';
backfile2 = '';
fn1 = '';
fn2 = '';
done1 = 0;
done2 = 0;

% Import plotter data
try; h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

settings = h.sen_set;
g = h.g;

if settings.ldar_tshiftOn
    x0 = settings.x(settings.ldar_tshift_sn);
    y0 = settings.y(settings.ldar_tshift_sn);
    z0 = settings.z(settings.ldar_tshift_sn);
else
    x0 = 0 ; y0 = 0; z0 = 0;
end
    

if g.mm < 30;     ext=0;
else   ext=30; end


switch dataType
    case 'PBFA-A'        
        pbfa_fn=sprintf('/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
            g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
            g.MM{:},g.DD{:},g.hh,ext);
    case 'PBFA'        
        pbfa_fn=sprintf('/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
            g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
            g.MM{:},g.DD{:},g.hh,ext);
    case 'PBFA-O'
        pbfa_fn=sprintf('/PBFA_old/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
            g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
            g.MM{:},g.DD{:},g.hh,ext);        
    otherwise
        disp('Point type coudn''t identified')
        return
end


if strcmp(delType,'Server') || strcmp(delType,'Local+Server')
    % Let's remove it from the local location
    [fn1, backfile1, done1] = delete_line([settings.base_dir pbfa_fn],t,x0,y0,z0);
    
    % If it is not done, let's remove the backup file
    if ~done1
        delete(backfile1)
    end
    
end

if strcmp(delType,'Local') || strcmp(delType,'Local+Server')
    % Let's remove it from the server location
    [fn2, backfile2, done2] = delete_line([settings.base_dir pbfa_fn],t,x0,y0,z0);
    
    % If it is not done, let's remove the backup file
    if ~done2
        delete(backfile2)
    end
    
end


function [fn, backfile, done] = delete_line(fn,t0,x0,y0,z0)

% Have I delete it
done = 0;

[pathstr, file, ext] = fileparts(fn);
backfile = sprintf('%s.backup.%4.4i%2.2i%2.2i_%2.2i%2.2i%2.2i',fn,round(clock));

%Create the backup the file
try 
    status = copyfile(fn,backfile);
    
    if ~status
        error('Coudn''t create the backup file. Point not deleted.')
    end        
catch
    error('Coudn''t create the backup file. Point not deleted.')
end
    
% Create a temporary file
tmpfile = [tempdir file ext];

fout = fopen(tmpfile,'w');

if fout < 1
    error('Coudn''t create the temorary file')
end

% Read the file and write the temporary file
fin = fopen(fn, 'r');

% read lines and write it in the out file

inline = fgets(fin);

while ischar(inline)
       
    % what is the the time
    vals = strsplit(inline);

    % Original time and location
    try 
        t = str2double(vals(2));
        x = str2double(vals(3));
        y = str2double(vals(4));
        z = str2double(vals(5));
        
        % arrival time
        t_ar = t + sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)/3e8;
        
        dt = abs(t_ar - t0);
        
        if dt < 0.1e-6
            disp('found the point and deleted.')
            done = 1;
        else
            fwrite(fout, inline);
        end
        
    catch
        % This line might be a blank. (let's write it as it is)
        fwrite(fout, inline);
    end
    
    % Read the next line
    inline = fgets(fin);
end

fclose(fin);
fclose(fout);

% Let's move the temporary file to the correct one.
if done 
    try
        status = copyfile(tmpfile,fn);
        
        if ~status
            error('Coudn''t update the file')
        end
    catch
        error('Coudn''t update the file.')
    end
else
    disp('The point not found')
end







    
    
    