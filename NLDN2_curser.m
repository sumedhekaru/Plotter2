function output_txt = NLDN2_curser(obj,event_obj)
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

nldn2_fn=sprintf('%s/NLDN2/%s/%s/NLDN2_%s%s%s.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},...
        g.MM{:},g.DD{:});

if exist(nldn2_fn,'file')
    
    if settings.ldar_tshiftOn
        x0 = settings.x(settings.ldar_tshift_sn);
        y0 = settings.y(settings.ldar_tshift_sn);
        z0 = settings.z(settings.ldar_tshift_sn);
    else
        x0 = 0 ; y0 = 0; z0 = 0;
    end
    
    % Load NLDN data
    [NLDN2c NLDN2g] = nldnExtract(nldn2_fn,0,0,t-3,t+3,...
        str2double(settings.ldar_r),x0,y0,z0,0);
    NLDN = [NLDN2c ; NLDN2g];
    
    NLDN = sortrows(NLDN,8);

    [m indx] = min(abs(NLDN(:,8)-t));
    %[m2 ind2] = min(abs(NLDN(:,4)-z))

    %indx = intersect(ind,ind2)

    if ~isempty(indx)
        ind = indx(1);
        output_txt = sprintf('T0  = %0.7fs\nX  = %0.1fm \nY  = %0.1fm \nZ  = %0.1fm \nIp = %0.1fkA\nt= %0.7f\n',...
            NLDN(ind,1),NLDN(ind,2),NLDN(ind,3),NLDN(ind,4),NLDN(ind,7),NLDN(ind,8));
        
        fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\n',...
            NLDN(ind,8),NLDN(ind,2),NLDN(ind,3),NLDN(ind,4),NLDN(ind,7))
    else
        output_txt = 'Coudn''t find a NLDN point';
    end

    %   data{1} - distance to the pulse from the given sensor 
    %   data{2} - Occuring time
    %   data{3} - x
    %   data{4} - y
    %   data{5} - z
    %   data{6} - Detection time at the sensor
    
%     if m < 3e-6
%         output_txt = sprintf('T  = %0.7fs\nX  = %0.1fm \nY  = %0.1fm \nZ  = %0.1fm \nIp = %0.1fA\nN  = %i\nt = %0.7fs\nki-sqrd = %0.1f',...
%             pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5),pbfa(ind,7),pbfa(ind,8),t,pbfa(ind,13));
%         
%         fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
%                pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5),pbfa(ind,8),pbfa(ind,7))
%     else
%         output_txt = 'Coudn''t find a PBFA point';
%     end
    
else
    output_txt = 'NLDN file not found!';
end