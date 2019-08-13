function dQ_finder4(a)
% Charge tranfer due to electrostatic change

if nargin == 0
    %% User inputs
    % Sensor numbers
    a.sns = [2 3 6];
    
    % measured dEs
    a.dEs = [0.9091 2.0104 0.4156];
    
    % Higher position [x y z]
    a.p1 = [1188.7 -5293.7 6200];
    
    % Lower position [x y z]
    a.p2 = [1188.7 -5293.7 5597.2];
    
    % lower bound tolerence of p2 (dx,dy,dz,dq)
    a.lb = [2000 2000 2000 10];
    
    % Upper bound tolerence of p2 (dx,dy,dz,dq)
    a.ub = [2000 2000 2000 10];
    
end


%% Start program

a.k = 1 / 2 / pi / 8.85418782e-12;

% Sensor positions
a.x = 1.0e+004*[ -0.752420000000000
    0.337160000000000
    -0.314850000000000
    -0.426630000000000
    -1.165790000000000
    -0.270095000000000
    -2.764040000000000
    -6.009070000000000
    0.182483000000000
    -5.739400000000000
    -2.063660000000000];


a.y = 1.0e+004 * [1.655450000000000
    0.444560000000000
    -0.683760000000000
    0.954490000000000
    -1.701960000000000
    0.263128000000000
    4.925430000000000
    -3.398320000000000
    -5.300820000000000
    1.192280000000000
    0.156929000000000];

%sensor IDS
sen_IDs = {'K02', 'K14', 'K24', 'WSB', 'BCC', 'K17' ,'EDW','STC','FLT','OVD'};
a.x = a.x(a.sns);
a.y = a.y(a.sns);

%% Optimize just the charge using Levenburg
[dQ ki dEc] = optimize_q(a);

fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('Calculated charge  \t\t= %0.4f C\n',dQ);
fprintf('Reduced Chi-squired \t= %0.2f\n',ki);
fprintf('Upper location \t\t\t= [%0.1f\t%0.1f\t%0.1f] m\n',a.p1);
fprintf('Lower location \t\t\t= [%0.1f\t%0.1f\t%0.1f] m\n',a.p2);

H = abs(a.p1(3)-a.p2(3))/1000;
fprintf('Charge Moment \t\t\t= %0.04f C km (Verticle Height = %0.2f km)\n',...
    2*dQ*H,H)
D = sqrt(sum((a.p1 - a.p2).^2))/1000;
fprintf('Charge Moment \t\t\t= %0.04f C km (Distance = %0.2f km)\n',...
    2*dQ*D,D)

%fprintf('Number of iterations\t= %i\n',N);
fprintf('_____________________________________________________\n')
fprintf('Sensor\t\tMeasured\t\tCalculated\t\t%%diff\n')
fprintf('=====================================================\n')
for i = 1:length(a.sns)
    fprintf('%s\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
        sen_IDs{a.sns(i)},a.dEs(i),dEc(i),abs(a.dEs(i)-dEc(i))/abs(a.dEs(i)+dEc(i))*200);
end


%% Levenberg-Marquardt for lower charge and its location

a.dQ = dQ;
[P2 ki2 dEc2] = optimize_p2(a);

fprintf('\nOptimized lower charge location using Levenberg-Marquardt\n')
fprintf('\nCalculated charge  \t\t= %0.4f C\n',P2(4));
fprintf('Reduced Chi-squired \t= %0.2f\n',ki2);
fprintf('Lower location \t\t\t= [%0.1f\t%0.1f\t%0.1f] m\n',P2(1:3));

H2 = abs(a.p1(3)-P2(3))/1000;
fprintf('Charge Moment \t\t\t= %0.04f C km (Verticle Height = %0.2f km)\n',...
    2*P2(4)*H2,H2)
D2 = sqrt(sum((a.p1 - P2(1:3)).^2))/1000;
fprintf('Charge Moment \t\t\t= %0.04f C km (Distance = %0.2f km)\n',...
    2*P2(4)*D2,D2)

%fprintf('Number of iterations\t= %i\n',N);
fprintf('_____________________________________________________\n')
fprintf('Sensor\t\tMeasured\t\tCalculated\t\t%%diff\n')
fprintf('=====================================================\n')
for i = 1:length(a.sns)
    fprintf('%s\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
        sen_IDs{a.sns(i)},a.dEs(i),dEc2(i),abs(a.dEs(i)-dEc2(i))/abs(a.dEs(i)+dEc2(i))*200);
end

%% Trust region reflective optimization for lower charge and its location

if length(a.sns) > 3
    [P2 ki2 dEc2] = optimize_p2_2(a);
else
    fprintf('\nCanot complete Trust-Region-Reflective algorithm as it requires at least as many \nequations as variables; aborting.\n')
    return
end

fprintf('\nOptimized lower charge location using Trust Region Reflective\n')
fprintf('\nCalculated charge  \t\t= %0.4f C\n',P2(4));
fprintf('Reduced Chi-squired \t= %0.2f\n',ki2);
fprintf('Lower location (P2) \t= [%0.1f\t%0.1f\t%0.1f] m\n',P2(1:3));
fprintf('P2 Lower tolence\t\t= [-%0.1f\t-%0.1f\t-%0.1f\t-%0.1f]\n',a.lb)
fprintf('P2 Upper tolence\t\t= [+%0.1f\t+%0.1f\t+%0.1f\t+%0.1f]\n',a.ub)

H2 = abs(a.p1(3)-P2(3))/1000;
fprintf('Charge Moment \t\t\t= %0.04f C km (Verticle Height = %0.2f km)\n',...
    2*P2(4)*H2,H2)
D2 = sqrt(sum((a.p1 - P2(1:3)).^2))/1000;
fprintf('Charge Moment \t\t\t= %0.04f C km (Distance = %0.2f km)\n',...
    2*P2(4)*D2,D2)

%fprintf('Number of iterations\t= %i\n',N);
fprintf('_____________________________________________________\n')
fprintf('Sensor\t\tMeasured\t\tCalculated\t\t%%diff\n')
fprintf('=====================================================\n')
for i = 1:length(a.sns)
    fprintf('%s\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
        sen_IDs{a.sns(i)},a.dEs(i),dEc2(i),abs(a.dEs(i)-dEc2(i))/abs(a.dEs(i)+dEc2(i))*200);
end

fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')



% outputs
% argout.dQ = dQ;
% argout.ki = ki;
% argout.dEc = dEc;
% argout.dQ2 = P2(4);
% argout.LP2 = P2(1:3);
% argout.ki2 = ki2;
% argout.dEc2 = dEc2;








function F = myfun(x,a)

% c = 299702547;
% %global t sns sen_set
%
% F = (c*(x(1)-t)).^2 - ((x(2) - sen_set.x(sns)).^2 + ...
%                       (x(3) - sen_set.y(sns)).^2 + ...
%                       (x(4) - sen_set.z(sns)).^2);

F = a.dEs' + a.k*x(4)*a.p1(3)./((a.x-a.p1(1)).^2+(a.y-a.p1(2)).^2+a.p1(3).^2).^1.5 - ...
    a.k*x(4)*x(3)./((a.x-x(1)).^2+(a.y-x(2)).^2+x(3).^2).^1.5;




function [P2 ki_sqrd dEc] = optimize_p2(a)
% Optimize the lower point using Levenburg-Marquite algorythm

x0 = [a.p2 a.dQ];

options=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.001},'TolFun',1e-20,'TolX',1e-20);
f = @(x) myfun(x,a);
%[x temp exitFlag] = fsolve(f,x0,options);
x = fsolve(f,x0,options);

L0 = length(a.sns);
if L0 > 5; np = L0 - 4;
else       np = 1;
end

% Measurement errors
mErr = [0.32 0.31 0.31 NaN 0.18 0.29 0.53 0.16 0.13 0.17 NaN];

dEc = -a.k*x(4)*a.p1(3)./((a.x-a.p1(1)).^2+(a.y-a.p1(2)).^2+a.p1(3).^2).^1.5 + ...
    a.k*x(4)*x(3)./((a.x-x(1)).^2+(a.y-x(2)).^2+x(3).^2).^1.5;

ki_sqrd = sum(((dEc' - a.dEs)./mErr(a.sns)).^2)/np;

P2 = x;

function [P2 ki_sqrd dEc] = optimize_p2_2(a)
% Optimize the lower point using Levenburg-Marquite algorythm

x0 = [a.p2 a.dQ];

options=optimset('Display','off','TolFun',1e-20,'TolX',1e-20);
f = @(x) myfun(x,a);
%[x temp exitFlag] = fsolve(f,x0,options);
%x = fsolve(f,x0,options);


x = lsqnonlin(f,x0,x0-a.lb,x0+a.ub, options);

L0 = length(a.sns);
if L0 > 5; np = L0 - 4;
else       np = 1;
end

% Measurement errors
mErr = [0.32 0.31 0.31 NaN 0.18 0.29 0.53 0.16 0.13 0.17 NaN];

dEc = -a.k*x(4)*a.p1(3)./((a.x-a.p1(1)).^2+(a.y-a.p1(2)).^2+a.p1(3).^2).^1.5 + ...
    a.k*x(4)*x(3)./((a.x-x(1)).^2+(a.y-x(2)).^2+x(3).^2).^1.5;

ki_sqrd = sum(((dEc' - a.dEs)./mErr(a.sns)).^2)/np;

P2 = x;

function F = myfun2(x,a)

F = a.dEs' + a.k*x*a.p1(3)./((a.x-a.p1(1)).^2+(a.y-a.p1(2)).^2+a.p1(3).^2).^1.5 - ...
             a.k*x*a.p2(3)./((a.x-a.p2(1)).^2+(a.y-a.p2(2)).^2+a.p2(3).^2).^1.5;
         
function [P2 ki_sqrd dEc] = optimize_q(a)
% Optimize the lower point using Levenburg-Marquite algorythm

x0 = 5.02;

options=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.001},'TolFun',1e-20,'TolX',1e-20);
f = @(x) myfun2(x,a);
%[x temp exitFlag] = fsolve(f,x0,options);
x = fsolve(f,x0,options);

L0 = length(a.sns);
if L0 > 1; np = L0 - 1;
else       np = 1;
end

% Measurement errors
mErr = [0.32 0.31 0.31 NaN 0.18 0.29 0.53 0.16 0.13 0.17 NaN];

dEc = -a.k*x*a.p1(3)./((a.x-a.p1(1)).^2+(a.y-a.p1(2)).^2+a.p1(3).^2).^1.5 + ...
       a.k*x*a.p2(3)./((a.x-a.p2(1)).^2+(a.y-a.p2(2)).^2+a.p2(3).^2).^1.5;

ki_sqrd = sum(((dEc' - a.dEs)./mErr(a.sns)).^2)/np;

P2 = x; 
