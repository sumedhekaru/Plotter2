function output_txt = LINET_curser(obj,event_obj)
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

linet_fn=sprintf('%s/LINET/%s/%s/%s/LINET_%s%s%s_%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
    g.MM{:},g.DD{:},g.hh,ext);


if exist(linet_fn,'file')
    
    if settings.ldar_tshiftOn
        x0 = settings.x(settings.ldar_tshift_sn);
        y0 = settings.y(settings.ldar_tshift_sn);
        z0 = settings.z(settings.ldar_tshift_sn);
    else
        x0 = 0 ; y0 = 0; z0 = 0;
    end
    
    % Load LINET data    
    linet = linetExtract2(linet_fn,t-300,t+300,str2double(settings.ldar_r),x0,y0,z0,0);
        
    % find closet LINET
    %x = pbfa(:,6)-t;
    [m ind] = min(abs(linet(:,3)-t));
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
    else
        output_txt = 'Coudn''t find a LINET point here';
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