function [x1,y1,z1,t1]=pbfa_finder(arg)
% clc

% Speed of light in air
v = 299792458/1.0003;

% Open sensor settings
sen_set = arg.sen_set; % open('sensor_setting.mat');

%%% Arguments %%%
% arg.inner = [ 1 2 3 6 11];      % Inner sensors
% arg.outer = [ 5 7 8 9 10];      % Outer sensors
% 
% arg.t_in  = getT.t_in;
% arg.t_out = getT.t_out;
% 
% arg.method = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of outer and inner sensors
n_o = numel(arg.outer);
n_i = numel(arg.inner);

%tic
% if arg.method == 3
%     % do nothing
% end

%toc
%tic

%method = arg.method;


switch arg.method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% METHOD 1 - 
%    Using over ditermined solution and considering innner and outer
%   pools together to find x and y. But only use the closest sensor to find
%   z coordinate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creating K and g arrays
all_sen = [arg.outer arg.inner];
t       = [arg.t_out arg.t_in ];

tmin    = min(t);
t = t - tmin;

x       = sen_set.x(all_sen);
y       = sen_set.y(all_sen);
z       = sen_set.z(all_sen);


for i=1:n_o+n_i-1;
    K(i,1) = x(i+1) - x(i);
    K(i,2) = y(i+1) - y(i);
    K(i,3) = z(i+1) - z(i);
    K(i,4) = ( t(i+1) - t(i) );
    
    g(i,1) = 0.5*(x(i+1)^2 + y(i+1)^2 + z(i+1)^2 ....
                - x(i)^2 - y(i)^2 - z(i)^2 ...
                - v^2*( t(i+1)^2-t(i)^2 ));
end



% f=K\g;
% f = pinv(K)*g;
[f,flag,relres]= lsqr(K,g,eps*2,100);

%flag

%relres

x1 = f(1);
y1 = f(2);
z1 = f(3);
t1 = -f(4)/v^2; 

% Find the closest sensor for the pulse
r = sqrt((x -x1).^2+ (y-y1).^2);
[rmin,ind] = min(r);

% Closest sensor number
csn = all_sen(ind);

z1 = sqrt(v^2*(t(ind)-t1)^2 - (x(ind)-x1)^2 - (y(ind) - y1)^2) - z(ind);

t1=t1+tmin;
% fprintf('\n****************************\n')
% fprintf('\tx = %0.1fm \n\ty = %0.1fm \n\tz = %0.1fm \n\tt = %0.7fs \n',x1,y1,z1,t1)
% fprintf('****************************\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% METHOD 2 - 
%    Using over ditermined solution and considering innner and outer
%   pools together to find x and y. But only use the closest TWO sensors to find
%   z coordinate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creating K and g arrays
all_sen = [arg.outer arg.inner];
t       = [arg.t_out arg.t_in ];

tmin    = min(t);
t = t - tmin;

x       = sen_set.x(all_sen);
y       = sen_set.y(all_sen);
z       = sen_set.z(all_sen);

for i=1:n_o+n_i-1;
    K(i,1) = x(i+1) - x(i);
    K(i,2) = y(i+1) - y(i);
    K(i,3) = z(i+1) - z(i);
    K(i,4) = ( t(i+1) - t(i) );
    
    g(i,1) = 0.5*(x(i+1)^2 + y(i+1)^2 + z(i+1)^2 ....
                - x(i)^2 - y(i)^2 - z(i)^2 ...
                - v^2*( t(i+1)^2-t(i)^2 ));
end


% f=K\g;
% f = pinv(K)*g;
[f,flag,relres]= lsqr(K,g,eps*2,100);

%flag

%relres

x1 = f(1);
y1 = f(2);
z1 = f(3);
t1 = -f(4)/v^2; 

% Find the closest sensor for the pulse
r = sqrt((x -x1).^2+ (y-y1).^2);
[rmin,ind] = min(r);

% Closest sensor number
%csn = all_sen(ind);

z11 = sqrt(v^2*(t(ind)-t1)^2 - (x(ind)-x1)^2 - (y(ind) - y1)^2) - z(ind);


% find 2nd closest sensor for the pulse
r(ind) = 1e10;
[rmin,ind] = min(r);
z12 = sqrt(v^2*(t(ind)-t1)^2 - (x(ind)-x1)^2 - (y(ind) - y1)^2) - z(ind);

% % find 3rd closest sensor for the pulse
% r(ind) = 1e10;
% [rmin,ind] = min(r);
% z13 = sqrt(v^2*(t(ind)-t1)^2 - (x(ind)-x1)^2 - (y(ind) - y1)^2) - z(ind);




z1 = (z11+z12)/2;
t1=t1+tmin;
% fprintf('\n****************************\n')
% fprintf('\tx = %0.1fm \n\ty = %0.1fm \n\tz = %0.1fm \n\tt = %0.7fs \n',x1,y1,z1,t1)
% fprintf('****************************\n\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   case 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% METHOD 3 - 
%    Using over ditermined solution and considering innner and outer
%   pools together to find x and y. But only use the closest available
%   sensor to find z coordinate. We will ignore sensors less than 6km
%   radius because closer sensor may have time shift problems due
%   induction term.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating K and g arrays
tic
all_sen = [arg.outer arg.inner];
t       = [arg.t_out arg.t_in ];

tmin    = min(t);
t = t - tmin;

x       = sen_set.x(all_sen);
y       = sen_set.y(all_sen);
z       = sen_set.z(all_sen);

K = nan(n_o+n_i-1,4);

for i=1:n_o+n_i-1;
    K(i,1) = x(i+1) - x(i);
    K(i,2) = y(i+1) - y(i);
    K(i,3) = z(i+1) - z(i);
    K(i,4) = ( t(i+1) - t(i) );
    
    g(i,1) = 0.5*(x(i+1)^2 + y(i+1)^2 + z(i+1)^2 ....
                - x(i)^2 - y(i)^2 - z(i)^2 ...
                - v^2*( t(i+1)^2-t(i)^2 ));
end

%f=K\g;
% f = pinv(K)*g;
[f,flag,relres]= lsqr(K,g,eps*2,100);

%flag

%relres

x1 = f(1);
y1 = f(2);
%z1 = f(3);
t1 = -f(4)/v^2; 

% Find the closest sensor for the pulse
r = sqrt((x -x1).^2+ (y-y1).^2);
%return

% Let's find all the sensors within 6-30 km
temp = [r' , all_sen' ,(1:1:length(r))'];
temp = sortrows(temp,1);

lol = sum(temp(:,1) < 6000) + 1;
ul  = sum(temp(:,1) < 30000);
%disp('working1')

if ul < lol
    %csn = temp(lol,2);
    ind = temp(lol,3);
else
    %csn = temp(lol:ul,2)';
    ind = temp(lol:ul,3)';
end

z1 = mean((v^2*(t(ind)-t1).^2 - (x(ind)-x1).^2 - (y(ind) - y1).^2).^0.5 - z(ind));

t1=t1+tmin;
% fprintf('\n****************************\n')
% fprintf('\tx = %0.1fm \n\ty = %0.1fm \n\tz = %0.1fm \n\tt = %0.7fs \n',x1,y1,z1,t1)
% fprintf('****************************\n\n')

%toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     case 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% METHOD 3 - 
%    Using over ditermined solution and considering innner and outer
%   pools together to find x and y. But only use the closest available
%   sensor to find z coordinate. We will ignore sensors less than 6km
%   radius because closer sensor may have time shift problems due
%   induction term.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creating K and g arrays
all_sen = [arg.outer arg.inner];
t       = [arg.t_out arg.t_in ];

tmin    = min(t);
t = t - tmin;

x       = sen_set.x(all_sen);
y       = sen_set.y(all_sen);
z       = sen_set.z(all_sen);

K = nan(n_o+n_i-1,4);

for i=1:n_o+n_i-1;
    K(i,1) = x(i+1) - x(i);
    K(i,2) = y(i+1) - y(i);
    K(i,3) = z(i+1) - z(i);
    K(i,4) = ( t(i+1) - t(i) );
    
    g(i,1) = 0.5*(x(i+1)^2 + y(i+1)^2 + z(i+1)^2 ....
                - x(i)^2 - y(i)^2 - z(i)^2 ...
                - v^2*( t(i+1)^2-t(i)^2 ));
end


% f=K\g;
% f = pinv(K)*g;
[f,flag,relres]= lsqr(K,g,eps*2,100);

%flag

%relres

x1 = f(1);
y1 = f(2);
%z1 = f(3);
t1 = -f(4)/v^2; 

% Find the closest sensor for the pulse
r = sqrt((x -x1).^2+ (y-y1).^2);
%return

% Let's find all the sensors more than 6- km
temp = [r' , all_sen' ,(1:1:length(r))'];
temp = sortrows(temp,1);

lol = sum(temp(:,1) < 6000) + 1;
% ul  = sum(temp(:,1) < 30000);
ul = length(temp);
%disp('working2')

if ul < lol
    %csn = temp(lol,2);
    ind = temp(lol,3);
else
    %csn = temp(lol:ul,2)';
    ind = temp(lol:ul,3)';
end

z1 = mean((v^2*(t(ind)-t1).^2 - (x(ind)-x1).^2 - (y(ind) - y1).^2).^0.5 - z(ind));

t1=t1+tmin;
% fprintf('\n****************************\n')
% fprintf('\tx = %0.1fm \n\ty = %0.1fm \n\tz = %0.1fm \n\tt = %0.7fs \n',x1,y1,z1,t1)
% fprintf('****************************\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% METHOD 5 - First estimate the x,y,z,t using Koshak's method and then use
% Levenberg-marquardt method to find the optimum solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arg.method = 3;
    [x0(2),x0(3),x0(4),x0(1)]=pbfa_finder(arg);
    
    if ~isreal(x0(4))
        x0(4) = 0;
    end
    
    
    sns = [arg.outer arg.inner];
    t       = [arg.t_out arg.t_in ];
    
    mm = min(t);
    t = (t - mm);
    x0(1) = x0(1) - mm;
        
   
    options=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.001},'TolFun',1e-20,'TolX',1e-20);
    f = @(x) myfun(x,t,sns,sen_set);
    %[x temp exitFlag] = fsolve(f,x0,options);
    x = fsolve(f,x0,options);
    
    t1 = x(1)+mm;    x1 = x(2);  y1 = x(3);  z1 = x(4);
   

    case 6
    % Especially decigned for to find RS locations. Here z will be ignored
    % (assume z = 0) and then calculate xyt using Koshak's method and then
    % use Levenburg method to optimize the answer
    
    % Creating K and g arrays
    all_sen = [arg.outer arg.inner];
    t       = [arg.t_out arg.t_in ];

    tmin    = min(t);
    t = t - tmin;

    x       = sen_set.x(all_sen);
    y       = sen_set.y(all_sen);
    z       = sen_set.z(all_sen);

    K = nan(n_o+n_i-1,3);

    for i=1:n_o+n_i-1;
        K(i,1) = x(i+1) - x(i);
        K(i,2) = y(i+1) - y(i);
        %K(i,3) = z(i+1) - z(i);
        K(i,3) = ( t(i+1) - t(i) );

        g(i,1) = 0.5*(x(i+1)^2 + y(i+1)^2 + z(i+1)^2 ....
                    - x(i)^2 - y(i)^2 - z(i)^2 ...
                    - v^2*( t(i+1)^2-t(i)^2 ));
    end


    f=K\g;
    % f = pinv(K)*g;
    % [f,flag,relres]= lsqr(K,g,eps*2,100);

    %flag

    %relres

    x0(2) = f(1);
    x0(3) = f(2);
    x0(1) = -f(3)/v^2; 
          
   
    options=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.001},...
        'TolFun',1e-20,'TolX',1e-20);
    f = @(x) myfunCG(x,t,all_sen,sen_set);
    %[x temp exitFlag] = fsolve(f,x0,options);
    try 
        x = fsolve(f,x0,options);
    catch
        x = [0 0 0];
    end
    
    t1 = x(1)+tmin;    x1 = x(2);  y1 = x(3);  z1 = 0;
    
end


function F = myfun(x,t,sns,sen_set)

c = 299702547;
%global t sns sen_set

F = (c*(x(1)-t)).^2 - ((x(2) - sen_set.x(sns)).^2 + ...
                      (x(3) - sen_set.y(sns)).^2 + ...
                      (x(4) - sen_set.z(sns)).^2);
                  
function F = myfunCG(x,t,sns,sen_set)

c = 299702547;
%global t sns sen_set

F = (c*(x(1)-t)).^2 - ((x(2) - sen_set.x(sns)).^2 + ...
                      (x(3) - sen_set.y(sns)).^2);