function pbfa_eror_finder(arg)
clc
% This program is written to find error finding of PBFA points. This will
% be a theoritical calculations rather than using actual Electric Field
% Change data.

arg.method      = 1;
arg.t_accuracy  = 0.0000002; % in seconds
arg.x_min       = -30000;    % in meters
arg.y_min       = -20000;
arg.x_max       = 10000;
arg.y_max       = 20000;
arg.z           = 5000;
arg.x_step      = 500;
arg.y_step      = arg.x_step;
arg.inner = [1 2 3 6 7 8 9 10];      % Inner sensors
arg.outer = [];      % Outer sensors

tic

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
    
    %figure
    a2.t_out = [];
    a2.method = 3;
    
    wbh = waitbar(0,'Calculating errors','Name','PBFA error plots');
    
    for i=1:nx
        a1.x0 = arg.x_min+(i-1)*arg.x_step;
        
        for j=1:ny
            
            a1.y0 = arg.y_min+(j-1)*arg.y_step;
            
            % find vertual times
            getT = pbfaVertualTimeGen(a1);
            
            %


            %a2.t_in  = round(getT.t_in/arg.t_accuracy)*arg.t_accuracy;
            %a2.t_out = round(getT.t_out/arg.t_accuracy)*arg.t_accuracy;
            %a2.method = 1;
            
            N = 100;
            x_vals = nan(1,N);
            y_vals = x_vals;
            z_vals = x_vals;
            t_vals = x_vals;
            
            for k=1:N
                
                try
                    %elapsed = toc;
                    val = ((i-1)*ny*N+(j-1)*N+k)/(nx*ny*N);
                    waitbar(val,wbh,sprintf('calculating errors.... (%0.3f%%)'...
                        ,val*100))
                catch
                    disp('User terminated the program!')
                    return
                end
                
                
                [r1 c1] = size(getT.t_in);
                %[r2 c2] = size(getT.t_out);
                a2.t_in = getT.t_in + normrnd(0,arg.t_accuracy,r1,c1);

            
                [x1,y1,z1,t1]=pbfa_finder(a2);
                
                x_vals(k) = x1;
                y_vals(k) = y1;
                z_vals(k) = z1;
                t_vals(k) = t1;
            end
            
            dx(i,j) = std(x_vals-a1.x0);   
            dy(i,j) = std(y_vals-a1.y0);
            dz(i,j) = std(z_vals-a1.z0);    
            dt(i,j) = std(t_vals);
            
            x(i,j) = a1.x0;     y(i,j) = a1.y0;     
            z(i,j) = a1.z0;     %t(i,j) = 0;            
            
             %hold all
             %eplot([x(i,j),dx(i,j)],[y(i,j),dy(i,j)],'K')
            
        end
    end
    
    delete(wbh)
    
    %xx=reshape(x,1,nx*ny);
    %yy=reshape(y,1,nx*ny);
    %dzz=reshape(dz,1,nx*ny);
        
%     hold all
%     scatter(xx,yy,25,dzz,'filled')
%     colorbar
  
    sen_set = open('sensor_setting.mat');
    
%     hold all
%     for i=1:20
%         if sen_set.x(i) ~= 0 && sen_set.y(i) ~= 0
%             plot(sen_set.x(i),sen_set.y(i),'ro','MarkerFaceColor','r')
%             text(sen_set.x(i),sen_set.y(i),sen_set.sen_IDs{i})
%         end
%     end
%     
%     daspect([1 1 1])
%     box on
%     %grid on
%     xlabel('East (km)')
%     ylabel('North (km)')
%     
    
    figure
    [ch ch]=contourf(x/1000,y/1000,dt*1e6);
    set(ch,'edgecolor','none')
    hold all
    florida_map
    daspect([1 1 1])
    c = colorbar;
    %sensor locations
    plot(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
        'ro','markerfacecolor','r')
    text(sen_set.x([1 2 3 5 6 7 8 9 10 11])/1000,sen_set.y([1 2 3 5 6 7 8 9 10 11])/1000,...
         {'K02','K14','K24','BCC','K17','EDW','STC','FLT','OVD','FFI'})
     xlabel('East (km)')
     ylabel('North (km)')
     %ylabel(c','Z-error (km)')
     ylabel(c','t-error (\mus)')
%     
%     fprintf('\n****************************\n')
%     fprintf('\tx = %0.1fm \n\ty = %0.1fm \n\tz = %0.1fm \n\tt = %0.7fs \n',x1,y1,z1,t1)
%     fprintf('****************************\n\n')
%     
   
    
end