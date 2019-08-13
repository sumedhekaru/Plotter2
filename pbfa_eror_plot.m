function pbfa_eror_plot(arg)

% This program is written to find error finding of PBFA points. This will
% be a theoritical calculations rather than using actual Electric Field
% Change data.

arg.method      = 2;
arg.t_accuracy  = 0.0000005; % in seconds
arg.x_min       = -30000;    % in meters
arg.y_min       = -20000;
arg.x_max       = 10000;
arg.y_max       = 20000;
arg.z           = 5000;
arg.x_step      = 5000;
arg.y_step      = arg.x_step;
% arg.inner       = [ 1 2 3 6 11];      % Inner sensors
% arg.outer       = [ 5 7 8 9 10];      % Outer sensors
arg.inner       = [ 1 2 3 5 6 7];      % Inner sensors
arg.outer       = [8 9 10 11];
arg.N           = 10;   % N


switch arg.method
    case 1
    % When time of arrivals have accuracy, we can rount time to match that
    % accuracy. When we do that recalculated position for PBFA is differ
    % from the acctual value. This method will consider the absolute value
    % of that difference as the error.
    
    nx = (arg.x_max - arg.x_min)/arg.x_step;
    ny = (arg.y_max - arg.y_min)/arg.y_step;
    
    a1.z0 = arg.z;
    a1.inner = arg.inner;
    a1.outer = arg.outer;
    
    a2.inner = arg.inner;
    a2.outer = arg.outer;
    
    dx = zeros(nx,ny);
    dy = dx;
    dz = dx;
    x = dx;
    y = dx;
    z = dx;
    
    figure
   
    for i=1:nx
        a1.x0 = arg.x_min+(i-1)*arg.x_step;
        
        for j=1:ny
            
            a1.y0 = arg.y_min+(j-1)*arg.y_step;
            
            % find vertual times
            getT = pbfaVertualTimeGen(a1);
            

            a2.t_in  = round(getT.t_in/arg.t_accuracy)*arg.t_accuracy;
            a2.t_out = round(getT.t_out/arg.t_accuracy)*arg.t_accuracy;
            a2.method = 1;
            
            [x1,y1,z1,t1]=pbfa_finder(a2);
            
            dx(i,j) = abs(x1-a1.x0);   
            dy(i,j) = abs(y1-a1.y0);
            dz(i,j) = abs(z1-a1.z0);    
            %dt(i,j) = abs(t1);
            
            x(i,j) = a1.x0;     y(i,j) = a1.y0;     
            z(i,j) = a1.z0;     % t(i,j) = 0;            
            
             hold all
             eplot([x(i,j),dx(i,j)],[y(i,j),dy(i,j)],'K')
            
        end
    end
    xx=reshape(x,1,nx*ny);
    yy=reshape(y,1,nx*ny);
    dzz=reshape(dz,1,nx*ny);
        
    hold all
    scatter(xx,yy,25,dzz,'filled')
    colorbar
  
    sen_set = open('sensor_setting.mat');
    
    hold all
    for i=1:20
        if sen_set.x(i) ~= 0 && sen_set.y(i) ~= 0
            plot(sen_set.x(i),sen_set.y(i),'ro','MarkerFaceColor','r')
            text(sen_set.x(i),sen_set.y(i),sen_set.sen_IDs{i})
        end
    end
    
    daspect([1 1 1])
    box on
    %grid on
    xlabel('East (km)')
    ylabel('North (km)')
%     
%     fprintf('\n****************************\n')
%     fprintf('\tx = %0.1fm \n\ty = %0.1fm \n\tz = %0.1fm \n\tt = %0.7fs \n',x1,y1,z1,t1)
%     fprintf('****************************\n\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
    % Here we do not round the theoritical time of arrival obtained.
    % Instead we add random error to between +/- 1us (if 1us is the timing
    % accuracy). Then we calculate the position and find the error in the
    % position. We can do 1000s of calculation to obtain an error
    % distribution and then use the standard deviation as the error of the
    % each coordinate.
    
    nx = (arg.x_max - arg.x_min)/arg.x_step;
    ny = (arg.y_max - arg.y_min)/arg.y_step;
    
    a1.z0 = arg.z;
    a1.inner = arg.inner;
    a1.outer = arg.outer;
    
    a2.inner = arg.inner;
    a2.outer = arg.outer;
    
    dx = zeros(nx,ny);
    dy = dx;
    dz = dx;
    x = dx;
    y = dx;
    z = dx;
    
    figure
    % temporary store values
    dx_temp = NaN(1,arg.N);
    dy_temp = NaN(1,arg.N);
    dz_temp = NaN(1,arg.N);
    dt_temp = NaN(1,arg.N);
    
    wbh = waitbar(0,'Wait','name','Error Cal Progress');
   
    for i=1:nx
        a1.x0 = arg.x_min+(i-1)*arg.x_step;
        
        for j=1:ny
            
            a1.y0 = arg.y_min+(j-1)*arg.y_step;
            
            % find vertual times
            getT = pbfaVertualTimeGen(a1);
            
            n1 = length(getT.t_in); 
            n2 = length(getT.t_out);
            
            a2.method = 1;
            
            getT.t_in = round(getT.t_in/arg.t_accuracy)*arg.t_accuracy;
            getT.t_out = round(getT.t_out/arg.t_accuracy)*arg.t_accuracy;
            
            
            for k = 1:arg.N
                try
                    str = sprintf('Working on %u of %i (Iteration #: %i)',((i-1)*ny+j),nx*ny,k);
                    waitbar(((i-1)*ny+j)/(nx*ny),wbh,str)
                catch
                    return
                end

                    
                % Add randowm numbers between +/-(1/2)*time_accuracy to it
                r1 = -arg.t_accuracy + 2*arg.t_accuracy.*rand(1,n1);
                r2 = -arg.t_accuracy + 2*arg.t_accuracy.*rand(1,n2);
                
                a2.t_in  = getT.t_in + r1;
                a2.t_out = getT.t_out + r2;
                
                [x1,y1,z1,t1]=pbfa_finder(a2);
                
                dx_temp(k) = abs(x1-a1.x0);
                dy_temp(k) = abs(y1-a1.y0);
                dz_temp(k) = abs(z1-a1.z0);
                dt_temp(k) = abs(t1);
                
            end
            
            dx(i,j) = std(dx_temp);   
            dy(i,j) = std(dy_temp);
            dz(i,j) = std(dz_temp);   
            % dt(i,j) = std(dt_temp);
            
            x(i,j) = a1.x0;     y(i,j) = a1.y0;     
            z(i,j) = a1.z0;     %t(i,j) = 0;            
            
             hold all
             eplot([x(i,j),dx(i,j)],[y(i,j),dy(i,j)],'K')
            
        end
    end
    max(max(dx))
    max(max(dy))
    delete(wbh)
    
    xx=reshape(x,1,nx*ny);
    yy=reshape(y,1,nx*ny);
    dzz=reshape(dz,1,nx*ny);
        
    hold all
%     scatter(xx,yy,25,dzz,'filled')
%     colorbar

    figure
    contourf(x,y,dz)
    
  
    sen_set = open('sensor_setting.mat');
    
    hold all
    for i=1:20
        if sen_set.x(i) ~= 0 && sen_set.y(i) ~= 0
            plot(sen_set.x(i),sen_set.y(i),'ro','MarkerFaceColor','r')
            text(sen_set.x(i),sen_set.y(i),sen_set.sen_IDs{i})
        end
    end
    
    daspect([1 1 1])
    box on
    %grid on
    str = sprintf('PBFA error astimation for z0 =%0.1fm & dt =%fs\nUsing monte-carlo method sn = %s',...
        arg.z,arg.t_accuracy,num2str([arg.inner arg.outer]));
    title(str)
    xlabel('East (km)')
    ylabel('North (km)')
    
    
  
%     figure
%     surf(x,y,dz)
    
%     figure
%     contourf(x,y,dz)
    
%     
%     fprintf('\n****************************\n')
%     fprintf('\tx = %0.1fm \n\ty = %0.1fm \n\tz = %0.1fm \n\tt = %0.7fs \n',x1,y1,z1,t1)
%     fprintf('****************************\n\n')
%     
        
    
end