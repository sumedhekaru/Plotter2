function n_CGLSS_finder3(handles)
%clc
g = handles.g;
settings = handles.sen_set;


if handles.temp.cglss
    cglss_fn=sprintf('%s/cglss/%s/%s/KSCCGLSS%s%s%s.dat',...
        settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},g.MM{:},g.DD{:});
    
    % Store blank if the file does not exsist
    if ~exist(cglss_fn,'file')
        cglss_fn = '';
    end
else
    cglss_fn = '';
end

n = handles.temp.N;
delT = handles.temp.delT;


% Time shift for CGLSS
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


if ~isempty(cglss_fn)
  
    % Load CGLSS data
    CG = CGLSS_extract(cglss_fn,g.t1,g.t2,str2double(settings.ldar_r),x0,y0,z0);
    
    t=CG(:,1);                       % time
    x=CG(:,2);                         % East
    y=CG(:,3);                         % North
    z=CG(:,12);                         % Altitude
    Ip = CG(:,11);                      % Peak Current
    %length(t)
    k = 1;
    T = [];
    
    if isempty(t)
        fprintf('No CGLSS in the given time range\n')
    else
        for i=1:length(t)
            t1=t(i)-delT;
            t2=t(i)+delT;
            %lol=size(CG(:,10),1)- nnz(CG(:,10)>t1)+1;      % Index of the lower matrix element
            lol = sum(t < t1)+1;
            ul = sum(t <= t2);                           % Index of the upper matrix element
            slice=t(lol:ul);                             % CG data for given time interval
            if length(slice)== n
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
        fprintf('\n\tTime(s)\t\t  Time(hh.mm.ss)\t  x(m) \t\t  y(m) \t\t  z(m) \t Ip(kA)\n')
        
        for i=1:k-1
            fprintf('%.6f \t %s \t %0.1f  \t %0.1f \t %0.1f \t %0.1f \n',...
                t(T(i)),sec2hhmmss(t(T(i))),x(T(i)),y(T(i)),z(T(i)),Ip(T(i)))

            if mod(i,n)==0 && n>1
                fprintf('\n')
            end
        end
    end
    
end