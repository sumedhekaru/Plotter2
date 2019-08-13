function output_txt = pbfaA_curser(obj,event_obj)
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

if g.mm < 30
    ext=0;
else
    ext=30;
end

pbfaA_fn=sprintf('%s/PBFA2/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);


if exist(pbfaA_fn,'file')
    
    if settings.ldar_tshiftOn
        x0 = settings.x(settings.ldar_tshift_sn);
        y0 = settings.y(settings.ldar_tshift_sn);
        z0 = settings.z(settings.ldar_tshift_sn);
    else
        x0 = 0 ; y0 = 0; z0 = 0;
    end
    
    % Load PBFA data    
    pbfa = pbfaExtract(pbfaA_fn,t-300,t+300,str2double(settings.ldar_r),x0,y0,z0,0);
    
    % find closet CGLSS point
    %cglss(:,13)
    %clc
    %x = pbfa(:,6)-t;
    %[m ind] = min(abs(pbfa(:,6)-t))
    %[m2 ind2] = min(abs(pbfa(:,5)-z))
    ind = find(pbfa(:,6) == t);
    ind2 = find(pbfa(:,5) == z);
    
    indx = intersect(ind,ind2);

    if ~isempty(indx)
        ind = indx(1);
        output_txt = sprintf('T  = %0.7fs\nX  = %0.1fm \nY  = %0.1fm \nZ  = %0.1fm \nIp = %0.1fkA\nN  = %i\nt = %0.7fs\nki-sqrd = %0.1f',...
            pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5),pbfa(ind,7),pbfa(ind,8),t,pbfa(ind,13));
        
        fprintf('%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
            pbfa(ind,2),pbfa(ind,3),pbfa(ind,4),pbfa(ind,5),pbfa(ind,8),pbfa(ind,7))
    else
        output_txt = 'Coudn''t find a PBFA point';
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
    output_txt = 'PBFA file not found!';
end