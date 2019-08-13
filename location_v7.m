function location_v7()
tic

%This program can use to find location of lightning
% Error propagation is done using 'MonteCarlo' with their equation

% Open Sensor Settings
a=open('sensor_setting.mat');

%Speed of light in air
v=299702547;

%Take t from the graph

prompt = a.sen_IDs;
dlg_title = 'Absolute time for a pulse';
num_lines = 1;

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

answer = inputdlg(prompt,dlg_title,num_lines,{'','','','',''},options);

if ~isempty(answer)
    t=str2double(answer);
end

monticalo=1;

% t=   1.0e+04 *[
%   7.321503253125001
%    7.321503257177000
%    7.321503259714000
%    7.321503254791000
%    7.321503263002000
% ];

% t=    1.0e+04 *[
% 
%    6.569605973870000
%    6.569605969550000
%    6.569605970716000
%    6.569605971483000
%    6.569605974908000
% 
%    ];

% t= 1.0e-04 *[
% 
%    0.394624786569437
%    0.279942138166907
%    0.399780392516067
%    0.153610341350816
%    0.767565686327346];



% coordinate errors of sensors
dx=1.0e0;          % error in x in meters
dy=1.0e0;          % errors in y in meters
dz=1.0e0;          % erors in z in meters
dt=1e-6;        % errors in retrival time in seconds


t0=min(t);
ts=t-t0;


xxx=[0,0,0,0,0];
yyy=[0,0,0,0,0];
ttt=[0,0,0,0,0];
%t555=[0,0,0,0,0];
zzz=[0,0,0,0,0];

delta_ts=zeros(5,5);

% Uncertainties
% dx=[0.1 0.1 0.1 0.1 0.1]*1e-1;
% dy=[0.1 0.1 0.1 0.1 0.1]*1e-1;
% dz=[0.1 0.1 0.1 0.1 0.1]*1e-1;
% dt=[0.1 0.1 0.1 0.1 0.1]*1e-6 ;
    
    


for i=1:5
    
    %i
    
    switch i
        
        case 1
            
            %t0=min([t(2),t(3),t(4),t(5)]);
            %ts=t-t0;
            t1=ts(2); t2=ts(3); t3=ts(4); t4=ts(5);
            x1=a.x(2); x2=a.x(3); x3=a.x(4); x4=a.x(5);
            y1=a.y(2); y2=a.y(3); y3=a.y(4); y4=a.y(5);
            z1=a.z(2); z2=a.z(3); z3=a.z(4); z4=a.z(5);
            
%             dx1=dx(2); dx2=dx(3); dx3=dx(4); dx4=dx(5);
%             dy1=dy(2); dy2=dy(3); dy3=dy(4); dy4=dy(5);
%             dz1=dz(2); dz2=dz(3); dz3=dz(4); dz4=dz(5);
%             dt1=dt(2); dt2=dt(3); dt3=dt(4); dt4=dt(5);
            
            % used sensors
            index(1,1:4)=[2 , 3, 4, 5];
            
            %fprintf('%i should be 1\n',i)
        case 2
            %t0=min([t(1),t(3),t(4),t(5)]);
            %ts=t-t0;
            t1=ts(1); t2=ts(3); t3=ts(4); t4=ts(5);
            x1=a.x(1); x2=a.x(3); x3=a.x(4); x4=a.x(5);
            y1=a.y(1); y2=a.y(3); y3=a.y(4); y4=a.y(5);
            z1=a.z(1); z2=a.z(3); z3=a.z(4); z4=a.z(5);
            
%             dx1=dx(1); dx2=dx(3); dx3=dx(4); dx4=dx(5);
%             dy1=dy(1); dy2=dy(3); dy3=dy(4); dy4=dy(5);
%             dz1=dz(1); dz2=dz(3); dz3=dz(4); dz4=dz(5);
%             dt1=dt(1); dt2=dt(3); dt3=dt(4); dt4=dt(5);
            
            % used sensors
            index(2,1:4)=[1 , 2, 3, 4];
            
            %fprintf('%i should be 2\n',i)
            
        case 3
            %t0=min([t(1),t(2),t(4),t(5)]);
            %ts=t-t0;
            t1=ts(1); t2=ts(2); t3=ts(4); t4=ts(5);
            x1=a.x(1); x2=a.x(2); x3=a.x(4); x4=a.x(5);
            y1=a.y(1); y2=a.y(2); y3=a.y(4); y4=a.y(5);
            z1=a.z(1); z2=a.z(2); z3=a.z(4); z4=a.z(5);
            
%             dx1=dx(1); dx2=dx(2); dx3=dx(4); dx4=dx(5);
%             dy1=dy(1); dy2=dy(2); dy3=dy(4); dy4=dy(5);
%             dz1=dz(1); dz2=dz(2); dz3=dz(4); dz4=dz(5);
%             dt1=dt(1); dt2=dt(2); dt3=dt(4); dt4=dt(5);
            
            % used sensors
            index(3,1:4)=[1 , 2, 4, 5];
            
            %fprintf('%i should be 3\n',i)
            
        case 4
            %t0=min([t(1),t(2),t(3),t(5)]);
            %ts=t-t0;
            t1=ts(1); t2=ts(2); t3=ts(3); t4=ts(5);
            x1=a.x(1); x2=a.x(2); x3=a.x(3); x4=a.x(5);
            y1=a.y(1); y2=a.y(2); y3=a.y(3); y4=a.y(5);
            z1=a.z(1); z2=a.z(2); z3=a.z(3); z4=a.z(5);
            
%             dx1=dx(1); dx2=dx(2); dx3=dx(3); dx4=dx(5);
%             dy1=dy(1); dy2=dy(2); dy3=dy(3); dy4=dy(5);
%             dz1=dz(1); dz2=dz(2); dz3=dz(3); dz4=dz(5);
%             dt1=dt(1); dt2=dt(2); dt3=dt(3); dt4=dt(5);
%             
            % used sensors
            index(4,1:4)=[1 , 2, 3, 5];
            
            %fprintf('%i should be 4\n',i)
            
        case 5
            %t0=min([t(1),t(2),t(3),t(4)]);
            %ts=t-t0;
            t1=ts(1); t2=ts(2); t3=ts(3); t4=ts(4);
            x1=a.x(1); x2=a.x(2); x3=a.x(3); x4=a.x(4);
            y1=a.y(1); y2=a.y(2); y3=a.y(3); y4=a.y(4);
            z1=a.z(1); z2=a.z(2); z3=a.z(3); z4=a.z(4);
            
%             dx1=dx(1); dx2=dx(2); dx3=dx(3); dx4=dx(4);
%             dy1=dy(1); dy2=dy(2); dy3=dy(3); dy4=dy(4);
%             dz1=dz(1); dz2=dz(2); dz3=dz(3); dz4=dz(4);
%             dt1=dt(1); dt2=dt(2); dt3=dt(3); dt4=dt(4);
            
            % used sensors
            index(5,1:4)=[1 , 2, 3, 4];
            
            %fprintf('%i should be 5\n',i)
            
           
        otherwise
            disp('Error Occured')
        return
    end
        
    
    K=[ x1-x2    y1-y2   t1- t2
        x2-x3    y2-y3   t2- t3
        x3-x4    y3-y4   t3- t4
        ];
    % Error in K
%     E=[ sqrt(dx1^2+dx2^2)    sqrt(dy1^2+dy2^2)  sqrt(dt1^2+dt2^2)
%         sqrt(dx2^2+dx3^2)    sqrt(dy2^2+dy3^2)  sqrt(dt1^2+dt2^2)
%         sqrt(dx3^2+dx4^2)    sqrt(dy3^2+dy4^2)  sqrt(dt1^2+dt2^2)
%         ];
% 
%     E=[ dx1-dx2  dy1-dy2     dt1-dt2
%         dx2-dx3  dy2-dy3     dt2-dt3
%         dx3-dx4  dy3-dy4     dt3-dt4];


    gg=[ x1^2+y1^2+z1^2-x2^2-y2^2-z2^2-v^2*(t1^2-t2^2)
        x2^2+y2^2+z2^2-x3^2-y3^2-z3^2-v^2*(t2^2-t3^2)
        x3^2+y3^2+z3^2-x4^2-y4^2-z4^2-v^2*(t3^2-t4^2)
        ]/2;
    
%     % Error in gg
%     delta =  [ (x1*dx1+y1*dy1+z1*dz1) - (x2*dx2+y2*dy2+z2*dz2) ...
%                     + 1/2*(sqrt(dx1^2+dy1^2+dz1^2)-sqrt(dx2^2+dy2^2+dz2^2)) ...
%                     - v^2*(sqrt(dt1^2+dt2^2)/2+(t1-t2)*sqrt(dt1^2+dt2^2))
%                     
%                (x2*dx2+y2*dy2+z2*dz2) - (x3*dx3+y3*dy3+z3*dz3) ...
%                     + 1/2*(sqrt(dx2^2+dy2^2+dz2^2)-sqrt(dx3^2+dy3^2+dz3^2)) ...
%                     - v^2*(sqrt(dt2^2+dt3^2)/2+(t2-t3)*sqrt(dt2^2+dt3^2))
%                     
%                (x3*dx3+y3*dy3+z3*dz3) - (x4*dx4+y4*dy4+z4*dz4) ...
%                     + 1/2*(sqrt(dx3^2+dy3^2+dz3^2)-sqrt(dx4^2+dy4^2+dz4^2)) ...
%                     - v^2*(sqrt(dt3^2+dt4^2)/2+(t3-t4)*sqrt(dt3^2+dt4^2))
%         
%                 ];
            
        
    %Find x,y and T ---- Sol(1) is x, Sol(2) is y and Sol(3) is related to T
    Sol=K\gg;
    f(1:3,i)=Sol;
    
    xxx(i)=Sol(1);
    yyy(i)=Sol(2);    
    ttt(i)=-Sol(3)/v^2 ;
    %t555(i)=t5;
    if i<5
        zzz(i)=abs(real(sqrt(v^2*(ts(i+1)-ttt(i))^2-(a.x(i+1)-Sol(1))^2-(a.y(i+1)-Sol(2))^2)));
    else
        zzz(i)=abs(real(sqrt(v^2*(ts(1)-ttt(i))^2-(a.x(1)-Sol(1))^2-(a.y(1)-Sol(2))^2)));
    end
                
    % Time difference from calculated and measured at each sensor in us
    delta_ts(i,:)=abs(ttt(i)+(sqrt((a.x-xxx(i)).^2+(a.y-yyy(i)).^2+(a.z-zzz(i)).^2))./v-ts')*1e6;   
    
    
    ttt(i)=ttt(i)+t0;
    
    
    % Finding Error
%     error1=(eye(3)+K\E)\(K\(E*Sol-delta));
%     
%     inv(K)*delta
%     
%     dxxx(i)=error1(1);
%     dyyy(i)=error1(2);
%     dttt(i)=-error1(3)/v^2;
    
    
end
% 
% dxxx=dxxx'
% dyyy=dyyy'
% dttt=dttt'


[best_dt best_index]=min(mean(delta_ts'));

%% Printing Results So far
% fprintf('Times\n')
fprintf('\n===============================================================')
Times=t

fprintf('Arrival time difference at each sensor (in us)\n')
fprintf('Excluded\t%s\t\t%s\t\t%s\t\t%s\t\t%s\n',...
    a.sen_IDs{1},a.sen_IDs{2},a.sen_IDs{3},a.sen_IDs{4},a.sen_IDs{5})


for i=1:5
    fprintf('%s\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n',...
        a.sen_IDs{i},delta_ts(i,1),delta_ts(i,2),delta_ts(i,3),...
        delta_ts(i,4),delta_ts(i,5))
end

fprintf('\nExcluded\tt (s)\t\tx (m)\t\ty (m)\t\tz (m)\n')

for i=1:5
        fprintf('%s\t\t%.8f\t%.1f\t\t%.1f\t\t%.1f\n',...
            a.sen_IDs{i},ttt(i),xxx(i),yyy(i),zzz(i))    
end

fprintf('Average\t\t%.8f\t%.1f\t\t%.1f\t\t%.1f\n\n',mean(ttt),mean(xxx),mean(yyy),mean(zzz))

%% Generating file name for LDAR data

global g

settings=open('sensor_setting.mat');

if g.mm < 30
    ext=0;
else
    ext=30;
end

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);



%Loading Ldar data
[CG,CAL,DLS]=ldarExtract(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
    0,0,0,0);

t1=[DLS(:,9);CG(:,9)];
x1=[DLS(:,6);CG(:,6)];
y1=[DLS(:,7);CG(:,7)];
z1=[DLS(:,8);CG(:,8)];


% print ldar info
if isempty(x1)==0
   fprintf('Near by LDAR points\n') 
    for i=1:length(x1)
        fprintf('\t\t%.6f\t%.1f\t\t%.1f\t\t%.1f\n',t1(i),x1(i),y1(i),z1(i))
    end
end

fprintf('\nBest Estimation has %.8f uS diff with calculation # %i\n',...
    best_dt*5,best_index) 


%% MonitCalo Error findings
if monticalo == 1
    
    n_steps = 10000;   % Number of steps
    
    
    % whitch sensor do you wish to omit
    cal_index=best_index;       % Use this enable if you wish do error propagation using best cal
    %cal_index=5;                 % Other wise use cal_index = 1 trough 5 
    
    %best_index=5;
    
%     % coordinate errors of sensors
%     dx=1.0e0;          % error in x in meters
%     dy=1.0e0;          % errors in y in meters
%     dz=1.0e0;          % erors in z in meters
%     dt=0.1e-6;        % errors in retrival time in micro seconds
    
    % Best astimations
    
    x_b=xxx(cal_index);    
    y_b=yyy(cal_index);      
    z_b=zzz(cal_index);    
    t_b=ttt(cal_index);
    
    % New astimations sores in the following variables
    x_a=zeros(1,n_steps);   y_a=x_a;
    z_a=x_a;                t_a=x_a;
    
    wbh=waitbar(0,'Moticarlo Error Anylisis...- 0%');
    
    % Good solution counter
    counter=0;
    
    for i=1:n_steps
        try
         mssg=sprintf('Monte-Carlo Error Anylisis...- %.1f%%',i/n_steps*100);
          waitbar(i/n_steps,wbh,mssg)
        catch
            fprintf('Analysis stopped by the user!\n')
            return
        end
        
         %Generate values from the uniform distribution on the interval [a, b].
         % r = a + (b-a).*rand(100,1);
         
         % coordintes uncertainties of the sensors
         dxns=(-dx)+2*dx.*rand(1,4);
         dyns=(-dy)+2*dy.*rand(1,4);
         dzns=(-dz)+2*dz.*rand(1,4);
         dtns=(-dt)+2*dt.*rand(1,4);
         
         % Error Matrix
         E=[ dxns(1)-dxns(2)  dyns(1)-dyns(2)     dtns(1)-dtns(2)
             dxns(2)-dxns(3)  dyns(2)-dyns(3)     dtns(2)-dtns(3)
             dxns(3)-dxns(4)  dyns(3)-dyns(4)     dtns(3)-dtns(4)
             ];
         
         x1=a.x(index(cal_index,1));       y1=a.y(index(cal_index,1)); z1=a.z(index(cal_index,1));   t1=t(index(cal_index,1));
         x2=a.x(index(cal_index,2));       y2=a.y(index(cal_index,2)); z2=a.z(index(cal_index,2));   t2=t(index(cal_index,2));
         x3=a.x(index(cal_index,3));       y3=a.y(index(cal_index,3)); z3=a.z(index(cal_index,3));   t3=t(index(cal_index,3));
         x4=a.x(index(cal_index,4));       y4=a.y(index(cal_index,4)); z4=a.z(index(cal_index,4));   t4=t(index(cal_index,4));
         
         
         
         % Find back the position of the pulese
         tns=[t1 t2 t3 t4];
         t0=min(tns);
         
         
         t1=t1-t0;
         t2=t2-t0;
         t3=t3-t0;
         t4=t4-t0;
         
         
         % Error in g column matrix           
         delta = [x1*dxns(1)+y1*dyns(1)+z1*dzns(1)-x2*dxns(2)-y2*dyns(2)-z2*dzns(2) + ...
                        1/2*(dxns(1)^2+dyns(1)^2+dzns(1)^2-dxns(2)^2-dyns(2)^2-dzns(2)^2) - ...
                        v^2*(t1*dtns(1)-t2*dtns(2)) + (v^2)/2*(dtns(1)^2-dtns(2)^2);
                        
                  x2*dxns(2)+y2*dyns(2)+z2*dzns(2)-x3*dxns(3)-y3*dyns(3)-z3*dzns(3) + ...
                        1/2*(dxns(2)^2+dyns(2)^2+dzns(2)^2-dxns(3)^2-dyns(3)^2-dzns(3)^2) - ...
                        v^2*(t2*dtns(2)-t3*dtns(3))+(v^2)/2*(dtns(2)^2-dtns(3)^2);
                        
                   x3*dxns(3)+y3*dyns(3)+z3*dzns(3)-x3*dxns(3)-y4*dyns(4)-z4*dzns(4) + ...
                        1/2*(dxns(3)^2+dyns(3)^2+dzns(3)^2-dxns(4)^2-dyns(4)^2-dzns(4)^2) - ...
                        v^2*(t3*dtns(3)-t4*dtns(4))+(v^2)/2*(dtns(3)^2-dtns(4)^2);                 
                    ];
             
                  
        

         
         K=[ x1-x2    y1-y2   t1- t2
             x2-x3    y2-y3   t2- t3
             x3-x4    y3-y4   t3- t4
             ];
         
         % Kin=inv(K);
         
%          gg=[ x1^2+y1^2+z1^2-x2^2-y2^2-z2^2-v^2*(t1^2-t2^2)
%              x2^2+y2^2+z2^2-x3^2-y3^2-z3^2-v^2*(t2^2-t3^2)
%              x3^2+y3^2+z3^2-x4^2-y4^2-z4^2-v^2*(t3^2-t4^2)
%              ]/2;
         
         %epsilon=inv(eye+Kin*E)*(Kin*E*f(:,cal_index)-Kin*delta);
         epsilon = K\(E*f(:,best_index)-delta);
         
        
         % Check whether the solution is real
         xcal=epsilon(1);
         ycal=epsilon(2);
         tcal=-epsilon(3)/v^2;
         zcal=(sqrt(v^2*(t2+t0-ttt(cal_index)-tcal)^2-(x2-xxx(cal_index)-xcal)^2-(y2-yyy(cal_index)-ycal)^2)+z2);
         %zcal=(sqrt(v^2*(t3-tcal)^2-(x3-xcal)^2-(y3-ycal)^2)+z2);
         %zcal=(sqrt(v^2*(t4-tcal)^2-(x4-xcal)^2-(y4-ycal)^2)+z2);
         %zcal=(sqrt(v^2*(t1-tcal)^2-(x1-xcal)^2-(y1-ycal)^2)+z2);
         
         %zcal=abs(real(sqrt(v^2*(t2-tcal)^2-(x2-xcal)^2-(y2-ycal)^2)+z2))
         
         if isreal(zcal)
             
             counter=counter+1;
             
             x_a(counter)=xcal;
             y_a(counter)=ycal;
             z_a(counter)=zcal;
             t_a(counter)=tcal;
             
         end
         
         %delta_ts=abs(tcal+(sqrt((xns-xcal).^2+(yns-ycal).^2+(zns-zcal).^2))./v-t')*1e6;
         
         
         
%          x_a(i)=Sol(1);
%          y_a(i)=Sol(2);
%          t_a(i)=-Sol(3)/v^2;
%          
%          z_a(i)=abs(real(sqrt(v^2*(t1-t_a(i))^2-(x1-x_a(i))^2-(y1-y_a(i))^2)));
%          
%          t_a(i)=t_a(i)+t0;    
%                          
    end
    
    % Statistics
    mean_dx=mean(x_a(1:counter));   std_dx=std(x_a(1:counter));   var_dx=var(x_a(1:counter));    min_dx=min(abs(x_a(1:counter)));    max_dx=max(abs(x_a(1:counter)));
    mean_dy=mean(y_a(1:counter));   std_dy=std(y_a(1:counter));   var_dy=var(y_a(1:counter));    min_dy=min(abs(y_a(1:counter)));    max_dy=max(abs(y_a(1:counter)));
    mean_dz=mean(z_a(1:counter)-zzz(best_index));   std_dz=std(z_a(1:counter)-zzz(best_index));   var_dz=var(z_a(1:counter)-zzz(best_index));    min_dz=min(abs(z_a(1:counter)-zzz(best_index)));    max_dz=max(abs(z_a(1:counter)-zzz(best_index)));
    mean_dt=mean(t_a(1:counter));   std_dt=std(t_a(1:counter));   var_dt=var(t_a(1:counter));    min_dt=min(abs(t_a(1:counter)));    max_dt=max(abs(t_a(1:counter)));

    
    delete(wbh)
end






%%%% printing statistical info   %%%%
if monticalo == 1
    fprintf('\nError Analysis Using Monti Carlo (%i/%i steps - %s Excluded)',counter,n_steps,a.sen_IDs{cal_index})
    fprintf('\n\t\tdx = %0.2fm ; \tdy = %0.2fm ; \tdz = %0.2fm \tdt = % 0.2fus ', dx,dy,dz,dt*1e6)
    fprintf('\n\nCoord \t\t mean \t\t STD \t\tVAR \t\t min \t\t max')
    fprintf('\ndx(m)\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f',mean_dx,std_dx,var_dx,min_dx,max_dx)
    fprintf('\ndy(m)\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f',mean_dy,std_dy,var_dy,min_dy,max_dy)
    fprintf('\ndz(m)\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f',mean_dz,std_dz,var_dz,min_dz,max_dz)
    fprintf('\ndt(us)\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n\n',mean_dt*1e6,std_dt*1e6,var_dt*1e12,min_dt*1e6,max_dt*1e6)
end











% figure
% plot3(xxx(1),yyy(1),zzz(1),'*');
% 
% hold all
% plot3(xxx(2),yyy(2),zzz(2),'*');
% plot3(xxx(3),yyy(3),zzz(3),'*');
% plot3(xxx(4),yyy(4),zzz(4),'*');
% plot3(xxx(5),yyy(5),zzz(5),'*');
% plot3(mean(xxx),mean(yyy),mean(zzz),'*');
% 
% if isempty(x1)==0
%     plot3(x1,y1,z1,'p')
%     legend('1','2','3','4','5','M','LDAR')
% else
%     legend('1','2','3','4','5','M')
% end
%     
% hold off
% 
% 
% box on
% grid on
% toc
fprintf('\n===============================================================\n')