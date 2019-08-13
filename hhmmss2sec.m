function sec = hhmmss2sec(in,outOff,plot_data)
% in - input time in 'hh:mm:ss.sss' format or hhmmss.sss format
% if outOff (optional) == 0 no output will produce
% if plot_data = 1, data arround "in" will be plotted according to plotter2
if ischar(in)
    inx = find(in == ':');
    hh = str2double(in(1:inx(1)-1));
    mm = str2double(in(inx(1)+1:inx(2)-1));
    ss_tmp = str2double(in(inx(2)+1:end));
    ss  = floor(ss_tmp);
    ss_dec = ss_tmp - ss;
else
    in_sec = floor(in);
    str=sprintf('%6.6i',in_sec);
    hh = str2double(str(1:2));
    mm = str2double(str(3:4));
    ss = str2double(str(5:6));
    ss_dec = in - in_sec;
    
end

sec = hh*3600+mm*60+ ss + ss_dec;

if nargin < 2 || nargin > 1 && outOff == 1  
    clipboard('copy', sec)
    fprintf('%2.2i:%2.2i:%2.2i.%7.7i = %6.6fs \n',hh,mm,ss,round(ss_dec*1e7),sec)
end



if nargin > 2
    if plot_data
        h=guidata(findall(0,'Tag','plotter2'));
        
        g = h.g;
        g.hh = hh;
        g.mm = floor(mm/5)*5;
        g.ss = 0;
        
        
        g.t1=sec - 0.05;
        g.t2=sec + 0.10;
        
        
        
        plot_all4(g)
        
        
        
        sen_set = h.sen_set;
        
        if g.mm < 30
            ext=0;
        else
            ext=30;
        end
        
        
        dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
            -datenum(str2double(g.YYYY),0,0));
        
        ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
            sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);
        
        
        linet_fn=[];
        
        nldn_fn=sprintf('%s/LINET2/%s/%s/%s/linet_%s%s%s_%2.2i%2.2i.txt',...
            sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
            g.MM{:},g.DD{:},g.hh,ext);
        
        pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
            sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
            g.MM{:},g.DD{:},g.hh,ext);
        
        ldarColorTime1(ldar_fn,linet_fn,pbfa_fn,nldn_fn,g.t1,g.t2,...
            str2double(sen_set.ldar_r),0,0,0,1,[0,0,1,0,0],sen_set);

        
        
    end
end
    
