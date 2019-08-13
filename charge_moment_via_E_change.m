function charge_moment_via_E_change(g)

if nargin < 1
    h=guidata(findall(0,'Tag','plotter2'));
    g = h.g;
    g.x0 = -23511.8;
    g.y0 = 13372.3;
    g.z0 = 0;
    
end

sen_set = open('sensor_setting.mat');
bfolder=sprintf('%s/%s/%s/%s/', ...
    sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:});
tshift=sen_set.t_shift;


% Load data
for i = 1:60
    if g.chgraphs(i)
        ext=mod(i,3);
        if ext==0
            ext=3;
        end
        sn = ceil(i/3);
        % Finding the stattion ID
        sid=sen_set.sen_IDs{sn};        
        
        % finding the file name
        filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        
        hfilename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        
        
        % Let's get the data
        [t, y ] = FA_Extract1(filename,hfilename,g.t1,g.t2,tshift(i),sen_set,i);
        
        % High pass filter
        [t, y] = ch_high_pass(t,y,1200);
        %yOffset = mean(y(1:1000));
        %y = y - yOffset;
        
        figure
        hold all
        plot(t,y)
        y1 = cumtrapz(t,y);
        plot(t,y1*1e6)
        r = sqrt((sen_set.x(sn)-g.x0)^2 + (sen_set.z(sn)-g.z0)^2 + (sen_set.z(sn)-g.z0)^2);
        P = 4*pi*8.85e-12*(3.0e8)^2*r*trapz(t,y1)
    end
end