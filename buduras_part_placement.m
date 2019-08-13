function buduras_part_placement
clc

%% Middle big circle
% ang = 250:-250/15:0;
% ang = ang-35;
% 
% % First raw
% rs = [0.9 1.2 1.5 1.8];
% LEDs = [1:56 377:384]; 
% for i = 1:length(ang)    
%     for j = 1:4
%         r = rs(j);    
%         x = r*cos(ang(i)*pi/180);
%         y = r*sin(ang(i)*pi/180);
%         fprintf('MOVE D%i (R %0.3f %0.3f)\n',LEDs((i-1)*4+j),x,y)
%     end
% end

%% Circles
% % Left most circle
% %x0 = -3.0311;
% %y0 = -1.75;
% %Ds = 57:120;
% 
% 
% % Right Most Circle
% %x0 = 3.0311;
% %y0 = -1.75;
% %Ds = 313:376;
% 
% % Left second circle
% %x0 = -3.5*cos(pi/6);
% %y0 = 3.5*sin(pi/6);
% %Ds = 121:184;
% 
% % Right second circle
% %x0 = 3.5*cos(pi/6);
% %y0 = 3.5*sin(pi/6);
% %Ds = 249:312;
% 
% % Top circle
% x0 = 0;
% y0 = 3.5;
% Ds = 185:248;
% 
% 
% rs = [0.5 0.9 1.3];
% Ns = [12 20 32];
% 
% k = 0;
% 
% for i=1:3
%     for j = 1:Ns(i)
%         k = k+1;
%         ang = 2*pi:-2*pi/Ns(i):0;
%         ang = ang+pi/2;
%         
%         x = x0 + rs(i)*cos(ang(j));
%         y = y0 + rs(i)*sin(ang(j));
%         
%         fprintf('MOVE D%i (R %8.3f %8.3f)\n',Ds(k),x,y);
%         
%     end
% end

%% Little curve segment on the perimeter
% 
% % Left most
% %ang0 = 180;
% %Ds = 481:488;
% 
% % Right Most
% %ang0 = 0;
% %Ds = 505:512;
% 
% % Left 2nd
% %ang0 = 120;
% %Ds = 489:496;
% 
% % Right 2nd
% ang0 = 60;
% Ds = 497:504;
% 
% thetas = 17.5:-35/8:-17.5;
% ang = thetas + ang0;
% r = 4.83;
% 
% k = 1;
% for i = 1:9
%     
%     x = r*cos(ang(i)*pi/180);
%     y = r*sin(ang(i)*pi/180);
%     
%     if ang(i) ~= ang0
%         fprintf('MOVE D%i (R %8.3f %8.3f)\n',Ds(k),x,y);
%         k = k+1;
%     end
% end

%% Little Stright strips
% 
% % Left most
% %ang = 180;
% %Ds = 449:459;
% 
% % Left 2nd
% %ang = 120;
% %Ds = 457:464;
% 
% % Right 2nd
% %ang = 60;
% %Ds = 465:472;
% 
% % Right most
% ang = 0;
% Ds = 473:480;
% 
% dL = 0.15;
% 
% rs = [2.1 2.35 2.6 2.85];
% 
% for i = 1:4
%     
%     x0 = rs(i)*cos(ang*pi/180);
%     y0 = rs(i)*sin(ang*pi/180);
%     
%     x = x0 + dL*cos((ang+90)*pi/180);
%     y = y0 + dL*sin((ang+90)*pi/180);
%     
%     fprintf('MOVE D%i (R %8.3f %8.3f)\n',Ds(i),x,y);
%     
%     x = x0 + dL*cos((ang-90)*pi/180);
%     y = y0 + dL*sin((ang-90)*pi/180);
%     
%     fprintf('MOVE D%i (R %8.3f %8.3f)\n',Ds(i+4),x,y);
% end

%% Diamond sections
dh = 0.25;
dv = 0.25;
ang = -26;

Ds = 385:400;

figure(100)
cla
hold all
%axis equal

xs = 0:dh:4*dh;
ys = 0:dv:4*dv;

for i = 1:4
    
    for j = 1:4
        x = xs(i)+(j-1)*dv;
        y = ys(j) ;
        
        r = sqrt(x^2+y^2);
        
        if r == 0
            ang0 = 0;
        else            
            ang0 = asin(y/r);
        end
        
        % new xy after rotation
        x = r*cos(ang0+ang*pi/180)
        y = r*sin(ang0+ang*pi/180);
        
        fprintf('MOVE D%i (R %8.3f %8.3f)\n',Ds((i-1)*4+j),x,y);
        plot(x,y,'ro','markerfacecolor','r')
    end
end


    
    
    
    


    