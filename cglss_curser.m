function output_txt = cglss_curser(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');


%t=sprintf('%.8f',pos(1));
t = pos(1);

% Open CGLSS file
% Import plotter data
try; h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

settings = h.sen_set;
g = h.g;


cglss_fn = sprintf('%s/cglss/%s/%s/KSCCGLSS%s%s%s.dat',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},g.MM{:},g.DD{:});

if exist(cglss_fn,'file')
    
    if settings.ldar_tshiftOn
        x0 = settings.x(settings.ldar_tshift_sn);
        y0 = settings.y(settings.ldar_tshift_sn);
        z0 = settings.z(settings.ldar_tshift_sn);
    else
        x0 = 0 ; y0 = 0; z0 = 0;
    end
    
    % Load CGLSS data
    cglss = CGLSS_extract(cglss_fn,t-1,t+1,str2double(settings.ldar_r),x0,y0,z0);
    
    % find closet CGLSS point
    %cglss(:,13)
    
    %x = cglss(:,13)-t;
    [m ind] = min(abs(cglss(:,13)-t));
    
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
    
    if m < 3e-6
        output_txt = sprintf('T  = %0.7fs\nX  = %0.1fm\nY  = %0.1fm \nIp = %0.1fA\nN  = %i\nt = %0.7fs',...
            cglss(ind,1),cglss(ind,2),cglss(ind,3),cglss(ind,11),cglss(ind,17),t);
        
         fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
               cglss(ind,1),cglss(ind,2),cglss(ind,3),0,cglss(ind,17),cglss(ind,11))
    else
        output_txt = 'Coudn''t find a CGLSS point';
    end
    
else
    output_txt = 'CGLSS file not found!';
end


