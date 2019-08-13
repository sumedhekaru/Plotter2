function output_txt = ldar2_curser(obj,event_obj)
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

if g.mm < 30
    ext=0;
else
    ext=30;
end

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
        -datenum(str2double(g.YYYY),0,0));
    
ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);


if exist(ldar_fn,'file')
    
    if settings.ldar_tshiftOn
        x0 = settings.x(settings.ldar_tshift_sn);
        y0 = settings.y(settings.ldar_tshift_sn);
        z0 = settings.z(settings.ldar_tshift_sn);
    else
        x0 = 0 ; y0 = 0; z0 = 0;
    end
    
    % Load PBFA data  
    [CG,CAL,ldar]= ldarExtract2(ldar_fn,t-300,t+300,str2double(settings.ldar_r),x0,y0,z0,0,0);
    
   
    %x = pbfa(:,6)-t;
    [m ind] = min(abs(ldar(:,10)-t));
    
    if m < 3e-6
        output_txt = sprintf('T  = %0.7fs\nX  = %0.1fm \nY  = %0.1fm \nZ  = %0.1fm \nt = %0.7fs',...
           ldar(ind,1),ldar(ind,6),ldar(ind,7),ldar(ind,8),t);

    else
        output_txt = 'Coudn''t find a LDAR2 point';
    end
    
else
    output_txt = 'LDAR2 file not found!';
end