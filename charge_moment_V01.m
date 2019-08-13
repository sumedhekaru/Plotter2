function charge_moment_V01(a)
% This function was writtent to find charge moment calculations due static
% changes of IC flashes. This will take electrostatic E-changes as input
% parameters and upper and lower charge locations.
%
%
% Modification History
%   2014-06-23 Created this function


clc
if nargin == 0
    %% User inputs
    % Sensor numbers
    a.sns = [1 2 3  6];
    
    % measured dEs
    a.dEs = [5.1771 8.8319 2.9915 5.1282];
    
    % Higher position [x y z]
    a.p1 = [23418.7 23295.2 9613];
    
    % Lower position [x y z]
    a.p2 = [23509.8 23353.6 7212.4];
    
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

[P2 ki_sqrd dEc] = linear_q(a)

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

function F = myfun2(x,a)

F = a.dEs' + a.k*x*a.p1(3)./((a.x-a.p1(1)).^2+(a.y-a.p1(2)).^2+a.p1(3).^2).^1.5 - ...
             a.k*x*a.p2(3)./((a.x-a.p2(1)).^2+(a.y-a.p2(2)).^2+a.p2(3).^2).^1.5;
         

function F = myfun_linear(x,a)

% Number of descrete poings
N = 100;
p2sx = a.p1(1):(a.p2(1)-a.p1(1))/N:a.p2(1);
p2sy = a.p1(2):(a.p2(2)-a.p1(2))/N:a.p2(2);
p2sz = a.p1(3):(a.p2(3)-a.p1(3))/N:a.p2(3);
qs = 1:-1/100:0;
factor = x /sum(qs);
qs = qs*factor;
L = length(a.x);
F = zeros(1,L);
for i = 1:L
    
    F(i) = a.dEs(i) + a.k*x*a.p1(3)./((a.x(i)-a.p1(1)).^2+(a.y(i)-a.p1(2)).^2+a.p1(3).^2).^1.5 - ...
        sum(a.k*qs.*p2sz./((a.x(i)-p2sx).^2+(a.y(i)-p2sy).^2+p2sz.^2).^1.5);
end
         
         
function [P2 ki_sqrd dEc] = linear_q(a)  

x0 = 1;

options=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.001},'TolFun',1e-20,'TolX',1e-20);
f = @(x) myfun_linear(x,a);
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
