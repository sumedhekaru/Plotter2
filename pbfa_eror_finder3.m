function pbfa_eror_finder3(arg)
clc
% This program is written to find error finding of PBFA points of suggested
% Olemiss array. This will
% be a theoritical calculations rather than using actual Electric Field
% Change data.

% Station	ShortName	Lat	Lon	x	y
% Olemiss	OLM	34.364732	-89.538443	0.0	0.0
% Field station of Olemiss	FST	34.432418	-89.390366	13606.8	7534.8
% UM Boonville Campus	BVL	34.667392	-88.564806	89467.4	33692.0
% UM Tupelo Campus	TPL	34.240198	-88.675344	79310.1	-13863.1
% UM Greneda Campus	GRD	33.763192	-89.818652	-25748.4	-66963.1
% UM Southheavan Campus	SHV	34.936553	-90.008993	-43238.8	63654.8





arg.method      = 1;
arg.t_accuracy  = 0.0000002; % in seconds
arg.x_min       = -50000;    % in meters
arg.y_min       = -50000;
arg.x_max       = 50000;
arg.y_max       = 50000;
arg.z           = 5000;
arg.x_step      = 1000;
arg.y_step      = arg.x_step;
arg.inner = [12:20];      % Inner sensors
arg.outer = [];      % Outer sensors

% Setup locations of current sensors
%arg.sen_set.x(15:20) = [0.0,13606.8,89467.4,79310.1,-25748.4,-43238.8];
%arg.sen_set.y(15:20) = [0.0,7534.8,33692.0,-13863.1,-66963.1,63654.8];
%arg.sen_set.z(15:20) = [0,0,0,0,0,0];
%station_names = {'OLM','FLS','BVL','TPL','GRD','SHV'};


% Setup locations
arg.sen_set.x(12:20) = [0.0,13.6068,30,  0, -30, -30, -30,   0,  30]*1000;
arg.sen_set.y(12:20) = [0.0,7.5348,  0, 30, 30,   0, -30, -30, -30]*1000;
arg.sen_set.z(15:20) = [0,0,0,0,0,0];
station_names = {'OLM','FLS','1','2','3','4','5','6'};

a2.sen_set = arg.sen_set;
a1.sen_set = arg.sen_set;

figure
plot(arg.sen_set.x,arg.sen_set.y,'ro')

% base folder and file name to save stuffs
bf = 'C:\Users\sumedhe\Desktop\Olemiss_array\';
fn = 'Hypothetical_30km_array';

sen_set = open('sensor_setting.mat');
v = 299792458/1.0003;
%a2.sen_set = sen_set;
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
    dy = dx;    dz = dx;    dt = dx;    ki_sqrd = dx;
    x = dx;    y = dx;    z = dx;
    x1 = dx;    y1 = dx;    z1 = dx;
    x2 = dx;    y2 = dx;    z2 = dx;
    x3 = dx;    y3 = dx;    z3 = dx;
    dx1 = dx;   dy1 = dx;    dz1 = dx;    dt1 = dx;    ki_sqrd1 = dx;
    dx2 = dx;   dy2 = dx;    dz2 = dx;    dt2 = dx;    ki_sqrd2 = dx;
    dx3 = dx;   dy3 = dx;    dz3 = dx;    dt3 = dx;    ki_sqrd3 = dx;

    
    %figure
    a2.t_out = [];
    
    
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
            
            N = 50;
            x_vals1 = nan(1,N);     y_vals1 = x_vals1;      
            z_vals1 = x_vals1;	    t_vals1 = x_vals1;
            
            x_vals2 = x_vals1;      y_vals2 = x_vals1;      
            z_vals2 = x_vals1;	    t_vals2 = x_vals1;
            
            x_vals3 = x_vals1;      y_vals3 = x_vals1;      
            z_vals3 = x_vals1;	    t_vals3 = x_vals1;
            
            
            val = ((i-1)*ny+j)/(nx*ny);
            waitbar(val,wbh,sprintf('calculating errors.... (%0.3f%%)'...
                ,val*100))

            
            for k=1:N                
%                 try
                    %elapsed = toc;
                    
%                 catch
%                     disp('User terminated the program!')
%                     return
%                 end
                
                
                [r1 c1] = size(getT.t_in);
                %[r2 c2] = size(getT.t_out);
                a2.t_in = getT.t_in + normrnd(0,arg.t_accuracy,r1,c1);

                a2.method = 3;
                [x1,y1,z1,t1]=pbfa_finder(a2);                
               
                
                x_vals1(k) = x1;
                y_vals1(k) = y1;
                z_vals1(k) = z1;
                t_vals1(k) = t1;
                
                a2.method = 4;
                [x1,y1,z1,t1]=pbfa_finder(a2);                
               
                
                x_vals2(k) = x1;
                y_vals2(k) = y1;
                z_vals2(k) = z1;
                t_vals2(k) = t1;
                
                a2.method = 5;
                [x1,y1,z1,t1]=pbfa_finder(a2);                
               
                
                x_vals3(k) = x1;
                y_vals3(k) = y1;
                z_vals3(k) = z1;
                t_vals3(k) = t1;
                
            end
            
                        
            dx1(i,j) = std(x_vals1-a1.x0);   
            dy1(i,j) = std(y_vals1-a1.y0);
            dz1(i,j) = std(z_vals1-a1.z0);    
            dt1(i,j) = std(t_vals1);
            ki_sqrd1(i,j) = (1/4)*sum(((getT.t_in - std(t_vals1) - sqrt((sen_set.x(arg.inner)-mean(x_vals1)).^2 + ...
                                                      (sen_set.y(arg.inner)-mean(y_vals1)).^2 +...
                                                      (sen_set.z(arg.inner)-mean(z_vals1)).^2)/v)/arg.t_accuracy).^2);
                     
            
                       
            dx2(i,j) = std(x_vals2-a1.x0);   
            dy2(i,j) = std(y_vals2-a1.y0);
            dz2(i,j) = std(z_vals2-a1.z0);    
            dt2(i,j) = std(t_vals2);
            ki_sqrd2(i,j) = (1/4)*sum(((getT.t_in - std(t_vals2) - sqrt((sen_set.x(arg.inner)-mean(x_vals2)).^2 + ...
                                                      (sen_set.y(arg.inner)-mean(y_vals2)).^2 +...
                                                      (sen_set.z(arg.inner)-mean(z_vals2)).^2)/v)/arg.t_accuracy).^2);
            
            dx2(i,j) = std(x_vals2-a1.x0);   
            dy3(i,j) = std(y_vals3-a1.y0);
            dz3(i,j) = std(z_vals3-a1.z0);    
            dt3(i,j) = std(t_vals3);
            ki_sqrd3(i,j) = (1/4)*sum(((getT.t_in - std(t_vals3) - sqrt((sen_set.x(arg.inner)-mean(x_vals3)).^2 + ...
                                                      (sen_set.y(arg.inner)-mean(y_vals3)).^2 +...
                                                      (sen_set.z(arg.inner)-mean(z_vals3)).^2)/v)/arg.t_accuracy).^2);
                                                  
            
%             if dz1(i,j) < dz2(i,j)
%                 dx(i,j) = dx1(i,j);
%                 dy(i,j) = dy1(i,j);
%                 dz(i,j) = dz1(i,j);
%                 dt(i,j) = dt1(i,j);
%                 ki_sqrd(i,j) = ki_sqrd1(i,j);
%                 method(i,j) = 3;
%             else
%                 dx(i,j) = dx2(i,j);
%                 dy(i,j) = dy2(i,j);
%                 dz(i,j) = dz2(i,j);
%                 dt(i,j) = dt2(i,j);
%                 ki_sqrd(i,j) = ki_sqrd2(i,j);
%                 method(i,j) = 4;
%             end
                                                  
                                                  
            
            x(i,j) = a1.x0;     y(i,j) = a1.y0;     
            z(i,j) = a1.z0;     t(i,j) = 0;
            
            x1(i,j) = mean(x_vals1);     y1(i,j) = mean(y_vals1);     
            z1(i,j) = mean(z_vals1);     t1(i,j) = mean(t_vals1);
            
            x2(i,j) = mean(x_vals2);     y2(i,j) = mean(y_vals2);     
            z2(i,j) = mean(z_vals2);     t2(i,j) = mean(t_vals2);
            
            x3(i,j) = mean(x_vals1);     y3(i,j) = mean(y_vals3);     
            z3(i,j) = mean(z_vals1);     t3(i,j) = mean(t_vals3);
            
            
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
    [ch ch]=contourf(x/1000,y/1000,dz3);
    set(ch,'edgecolor','none')
    hold all
    %florida_map
    daspect([1 1 1])
    c = colorbar;
    %sensor locations
    plot(a2.sen_set.x(arg.inner)/1000,a2.sen_set.y(arg.inner)/1000,...
        'ro','markerfacecolor','r')
    text(a2.sen_set.x(arg.inner)/1000,a2.sen_set.y(arg.inner)/1000,...
         station_names)
     xlabel('East (km)')
     ylabel('North (km)')
     ylabel(c','Z-error (km)')
     %ylabel(c','t-error (\mus)')
     
     saveas(gcf,[bf fn '.fig'])
     
     output.x = x;     output.y = y;     output.z = z;  output.t = t;
     output.x1 = x1;     output.y1 = y1;     output.z1 = z1;  output.t1 = t1;
     output.x2 = x2;     output.y2 = y2;     output.z2 = z2;  output.t2 = t2;
     output.x3 = x3;     output.y3 = y3;     output.z3 = z3;  output.t3 = t3;
     
     %output.dt = dt;   output.dx = dx;   output.dz = dz;    output.ki_sqrd = ki_sqrd;
     output.dt1 = dt1;   output.dx1 = dx1;   output.dz1 = dz1;    output.ki_sqrd1 = ki_sqrd1;
     output.dt2 = dt2;   output.dx2 = dx2;   output.dz2 = dz2;    output.ki_sqrd2 = ki_sqrd2;
     output.dt3 = dt3;   output.dx3 = dx3;   output.dz3 = dz3;    output.ki_sqrd3 = ki_sqrd3;
    
     output.a2 = a2;
     % output.method = method;
     [bf fn '.mat']
     save([bf fn '.mat'],'-Struct','output')
     
%     
%     fprintf('\n****************************\n')
%     fprintf('\tx = %0.1fm \n\ty = %0.1fm \n\tz = %0.1fm \n\tt = %0.7fs \n',x1,y1,z1,t1)
%     fprintf('****************************\n\n')
%     
   
    
end