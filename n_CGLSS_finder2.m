function n_CGLSS_finder2(handles)

g = handles.g;
settings = handles.sen_set;



% Generating file name for LDAR data
if g.mm < 30
    ext=0;
else
    ext=30;
end

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

if handles.temp.cglss
    ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);
    
    % Store blank if the file does not exsist
    if ~exist(ldar_fn,'file')
        ldar_fn = '';
    end
else
    ldar_fn = '';
end

n = handles.temp.N;
delT = handles.temp.delT;


% Time shift for LDAR and Linet
if settings.ldar_tshiftOn==1
    sn=settings.ldar_tshift_sn;
    x0=settings.x(sn)-settings.x0;
    y0=settings.x(sn)-settings.y0;
    z0=settings.x(sn)-settings.z0;
    %name = settings.sen_IDs{sn};
    
else
    x0=0;
    y0=0;
    z0=0;
    %name = 'None';
end


if ~isempty(ldar_fn)
  
    % Load LDAR data
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
        x0,y0,z0,0);
    
    t=CG(:,1);                        % time
    x=CG(:,6);                         % East
    y=CG(:,7);                         % North
    z=CG(:,8);                         % Altitude

    k = 1;
    T = [];
    
    if isempty(t)
        fprintf('No CGLSS in the given time range\n')
    else
        for i=1:length(t)
            t1=t(i)-delT;
            t2=t(i)+delT;
            %lol=size(CG(:,10),1)- nnz(CG(:,10)>t1)+1;      % Index of the lower matrix element
            lol = nnz(CG(:,10) < t1)+1;
            ul=nnz(CG(:,10) <= t2);                           % Index of the upper matrix element
            slice=t(lol:ul);                             % CG data for given time interval
            
            if size(slice)== n
                T(k)=i;
                k=k+1;
            end
        end
    end

    if isempty(T)
        fprintf('\nNo (%d strock/s) CGLSS were found (Time threshold = %0.3fs)\n',n,delT)
    else
        
        % Plot all graphs
        for i = 1:n:k-1
            g.t1 = t(T(i))- delT;
            g.t2 = g.t1  + 2*delT;
            plot_all3(g);            
        end
        
        
        % print times on the command windows
        fprintf('\n%d (%d strock/s) CGLSS were found (%0.2f percent,dT= %0.2f s)\n',ceil((k-1)/n),n,(k-1)/length(t)*100,delT)
        fprintf('\n\tTime\t\t  Time(hh.mm.ss)\t  x(m) \ty(m) \tz(m) \n')
        
        for i=1:k-1
            fprintf('%.6f s\t %s \t %0.1f %0.1f %0.1f m\n',t(T(i)),sec2hhmmss(t(T(i))),x(T(i)),y(T(i)),z(T(i)))

            if mod(i,n)==0 && n>1
                fprintf('\n')
            end
        end
    end
    
end