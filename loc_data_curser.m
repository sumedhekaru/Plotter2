function output_txt = loc_data_curser(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');


%t=sprintf('%.8f',pos(1));
t = pos(1);
z = pos(2);

% Open CGLSS file
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

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
        -datenum(str2double(g.YYYY),0,0));


%% Assume it is a LDAR2 point
ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);


if exist(ldar_fn,'file')
    
    % Load PBFA data  
    [CG,CAL,ldar]= ldarExtract2(ldar_fn,t-1,t+1,str2double(settings.ldar_r),x0,y0,z0,0,0);
       
    ind = find(ldar(:,10)==t);
    ind2 = find(ldar(:,8)==z);

    indx = intersect(ind,ind2);
    
    if ~isempty(indx)
        ind = indx(1);
        output_txt = sprintf('T  = %0.7fs\nX  = %0.1fm \nY  = %0.1fm \nZ  = %0.1fm \nt = %0.7fs',...
           ldar(ind,1),ldar(ind,6),ldar(ind,7),ldar(ind,8),t);
       
       fprintf('%13.7f\t%8.1f\t%8.1f\t%8.1f\n',ldar(ind,1),ldar(ind,6),ldar(ind,7),ldar(ind,8))
       strClp = sprintf('%13.7f\t%8.1f\t%8.1f\t%8.1f',ldar(ind,1),ldar(ind,6),ldar(ind,7),ldar(ind,8));
       clipboard('copy', strClp)
       return
    end
end

%% Assume it is a PBFA Auto point
pbfa_fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);


if exist(pbfa_fn,'file')
    
    % Load PBFA data    
    pbfa = pbfaExtract(pbfa_fn,t-1,t+1,str2double(settings.ldar_r),x0,y0,z0,0);
    
    ind = find(pbfa(:,6)==t);
    ind2= find(pbfa(:,5)==z);

    indx = intersect(ind,ind2);

    if ~isempty(indx)
        ind = indx(1);
        output_txt = sprintf('T  = %0.7fs\nX  = %0.1fm \nY  = %0.1fm \nZ  = %0.1fm \nIp = %0.1fA\nN  = %i\nt = %0.7fs\nki-sqrd = %0.1f',...
            pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5),pbfa(ind,7),pbfa(ind,8),t,pbfa(ind,13));
        
        %fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
        %    pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5),pbfa(ind,8),pbfa(ind,7))
        
       fprintf('%13.7f\t%8.1f\t%8.1f\t%8.1f\n',pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5))
       strClp = sprintf('%13.7f\t%8.1f\t%8.1f\t%8.1f',pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5));
       clipboard('copy', strClp)
        
        return
    end    
end


%% Assume it is a PBFA point
pbfa_fn=sprintf('%s/PBFA_old/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);


if exist(pbfa_fn,'file')
    
    % Load PBFA data    
    pbfa = pbfaExtract(pbfa_fn,t-1,t+1,str2double(settings.ldar_r),x0,y0,z0,0);
    
    ind = find(pbfa(:,6)==t);
    ind2= find(pbfa(:,5)==z);

    indx = intersect(ind,ind2);

    if ~isempty(indx)
        ind = indx(1);
        output_txt = sprintf('T  = %0.7fs\nX  = %0.1fm \nY  = %0.1fm \nZ  = %0.1fm \nIp = %0.1fA\nN  = %i\nt = %0.7fs\nki-sqrd = %0.1f',...
            pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5),pbfa(ind,7),pbfa(ind,8),t,pbfa(ind,13));
        
        fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
            pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5),pbfa(ind,8),pbfa(ind,7))        
               
       strClp = sprintf('%13.7f\t%8.1f\t%8.1f\t%8.1f',pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5));
       clipboard('copy', strClp)
        return
    end    
end

%% Assume it is a PBFA 
pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);


if exist(pbfa_fn,'file')
    
    % Load PBFA data    
    pbfa = pbfaExtract(pbfa_fn,t-1,t+1,str2double(settings.ldar_r),x0,y0,z0,0);
    
    ind = find(pbfa(:,6)==t);
    ind2= find(pbfa(:,5)==z);

    indx = intersect(ind,ind2);

    if ~isempty(indx)
        ind = indx(1);
        output_txt = sprintf('T  = %0.7fs\nX  = (%0.1f+-%0.1f)m \nY  = (%0.1f+-%0.1f)m \nZ  = (%0.1f+-%0.1f)m \nIp = %0.1fA\nN  = %i\nt = %0.7fs\nki-sqrd = %0.1f',...
            pbfa(ind,2),pbfa(ind,3),pbfa(ind,10),pbfa(ind,4),pbfa(ind,11),pbfa(ind,5),pbfa(ind,12),pbfa(ind,7),pbfa(ind,8),t,pbfa(ind,13));
        
        fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
            pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5),pbfa(ind,8),pbfa(ind,7))        
               
       strClp = sprintf('%13.7f\t%8.1f\t%8.1f\t%8.1f',pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5));
       clipboard('copy', strClp)
        return
    end    
end

%% Assume it is a Linet Point
linet_fn=sprintf('%s/LINET/%s/%s/%s/LINET_%s%s%s_%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);


if exist(linet_fn,'file')
    
    % Load LINET data    
    linet = linetExtract2(linet_fn,t-1,t+1,str2double(settings.ldar_r),x0,y0,z0,0);
        
    % find closet LINET
    %x = pbfa(:,6)-t;
    ind = find(linet(:,3)==t);
    ind2 = find(linet(:,8)== z);

    indx = intersect(ind,ind2);

    if ~isempty(indx)
        ind = indx(1);
        output_txt = sprintf('T  = %0.7fs\nX  = %0.1fm \nY  = %0.1fm \nZ  = %0.1fm \nType = %i\nIp  = %0.1f\nt = %0.7fs\nerr = %0.3f',...
            linet(ind,1),linet(ind,6),linet(ind,7),linet(ind,8),linet(ind,9),linet(ind,10),t,linet(ind,11));
        
        str = sprintf('%13.7f\t%s\t%7.4f\t%8.4f\t%4.1f\t%i\t%6.1f\t%6.3f\t%0.3f\t%0.3f',...
            linet(ind,1),sec2hhmmss(linet(ind,1)),linet(ind,13),linet(ind,12),...
            linet(ind,8)/1000,linet(ind,9),linet(ind,10),linet(ind,11),linet(ind,6)/1000,linet(ind,7)/1000);

        fprintf('%s\n',str)
        clipboard('copy', str)
        return
    end

    %   data{1} - distance to the pulse from the given sensor 
    %   data{2} - Occuring time
    %   data{3} - x
    %   data{4} - y
    %   data{5} - z
    %   data{6} - Detection time at the sensor
    
end

%% Assume it is a CGLSS point

cglss_fn = sprintf('%s/cglss/%s/%s/KSCCGLSS%s%s%s.dat',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},g.MM{:},g.DD{:});

if exist(cglss_fn,'file')
    
    % Load CGLSS data
    cglss = CGLSS_extract(cglss_fn,t-1,t+1,str2double(settings.ldar_r),x0,y0,z0);
    
    % find closet CGLSS point
    %cglss(:,13)
    
    %x = cglss(:,13)-t;
    [m ind] = min(abs(cglss(:,13)-t));
    ind2 = find(cglss(:,12)==z);

    indx = intersect(ind,ind2);    
    
    % 1 : Occuring time
    % 2 : x
    % 3 : y
    % 12: z
    % 11: Current
    % 4 : distance to x0,y0,z0 given
    % 5,6,7 - Latitude
    % 8,9,10 - Logitude
    % 13 - detection time
    % 17 - Total number of sensors used
    
    if ~isempty(indx)
        ind = indx(1);
        output_txt = sprintf('T  = %0.7fs\nX  = %0.1fm\nY  = %0.1fm \nIp = %0.1fA\nN  = %i\nt = %0.7fs',...
            cglss(ind,1),cglss(ind,2),cglss(ind,3),cglss(ind,11),cglss(ind,17),t);
        
         fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
               cglss(ind,1),cglss(ind,2),cglss(ind,3),0,cglss(ind,17),cglss(ind,11))
        return    
    end
end

% If it come to this far, it means the program coudn't find the location.
output_txt = 'Couldn''t find the loction';



