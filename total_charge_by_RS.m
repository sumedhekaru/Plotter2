function total_charge_by_RS(a)
% This function was writtent to find calculate charge/charge moment calculations
% dyue static changes of RSs.

% Modification History
%   2014-08-12 Created this function

if nargin == 0
    %% User inputs
    % Sensor numbers
    %a.sns = [1 2 3 5 6 10];
	a.sns = [1 2 3 6];
    
    % measured dEs
    %a.dEs = [-16.25 8.8319 2.9915 5.1282];
    %a.dEs = [-16.2857 -5.2571 -5.3571 -3.4714 -7.8000 -3.0286 ];
	%a.dEs = [-16.2857 -5.2571 -5.3571 -7.8000 ];
    a.dEs = [-15.9000 -33.4000 -30.3000  -85.9000];
    
    % Trial charge location and charge [x y z q] (First PBFA/LDAR2 location
    % may be apropriate here)
    %a.p = [-25328.9 11707.8 6047.0 -2];
    a.p = [-12444.2	  2072.6	  6194.6 1];
    
    % lower bound tolerence of p (dx,dy,dz,dq)
    a.lb = [5000 2000 5000 10];
    
    % Upper bound tolerence of p (dx,dy,dz,dq)
    a.ub = [7000 2000 5000 10];
    
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
sen_IDs = {'K02', 'K14', 'K24', 'WSB', 'BCC', 'K17' ,'EDW','STC','FLT','OVD','FFI'};
a.x = a.x(a.sns);
a.y = a.y(a.sns);

%% Optimize just the charge using Levenburg


[qH, ki, dEc] = optimize_q(a);

fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('Optimized using LM algorithm and forcing xy locations \nprovided and allowing q and H to change\n')
fprintf('\nCalculated charge  \t\t= %0.4f C\n',-qH(1));
fprintf('Charge Altititude  \t\t= %0.1f m\n',qH(2));

fprintf('Reduced Chi-squired \t= %0.2f\n',ki);
fprintf('Upper location \t\t\t= [%0.1f\t%0.1f\t%0.1f] m\n',a.p(1:2),qH(2));

H = abs(qH(2))/1000;
fprintf('Charge Moment \t\t\t= %0.04f C km (Vertical Height = %0.2f km)\n',...
    -2*qH(1)*H,H)


%fprintf('Number of iterations\t= %i\n',N);
fprintf('_____________________________________________________\n')
fprintf('Sensor\t\tMeasured\t\tCalculated\t\t%%diff\n')
fprintf('=====================================================\n')
for i = 1:length(a.sns)
    fprintf('%s\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
        sen_IDs{a.sns(i)},a.dEs(i),dEc(i),abs(a.dEs(i)-dEc(i))/abs(a.dEs(i)+dEc(i))*200);
end

[x, ki, dEc] = optimize_q2(a);

fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('Optimized using TR algorithm allowing xyzq to change\nin the given range.\n')
fprintf('\nCalculated charge  \t\t= %0.4f C\n',-x(4));
fprintf('Reduced Chi-squired \t= %0.2f\n',ki);
fprintf('Upper location \t\t\t= [%0.1f\t%0.1f\t%0.1f] m\n',x(1:3));

H = abs(x(3))/1000;
fprintf('Charge Moment \t\t\t= %0.04f C km (Vertical Height = %0.2f km)\n',...
    -2*x(4)*H,H)


%fprintf('Number of iterations\t= %i\n',N);
fprintf('_____________________________________________________\n')
fprintf('Sensor\t\tMeasured\t\tCalculated\t\t%%diff\n')
fprintf('=====================================================\n')
for i = 1:length(a.sns)
    fprintf('%s\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
        sen_IDs{a.sns(i)},a.dEs(i),dEc(i),abs(a.dEs(i)-dEc(i))/abs(a.dEs(i)+dEc(i))*200);
end

%% LM optimization
function [x, ki_sqrd, dEc] = optimize_q(a)
% Optimize the lower point using Levenburg-Marquite algorythm

x0 = [1 5000];

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

dEc = a.k*x(1)*x(2)./((a.x-a.p(1)).^2+(a.y-a.p(2)).^2+x(2).^2).^1.5;

ki_sqrd = sum(((dEc' - a.dEs)./mErr(a.sns)).^2)/np;


function F = myfun2(x,a)

F = ((a.dEs' - a.k*x(1)*x(2)./((a.x-a.p(1)).^2+(a.y-a.p(2)).^2+x(2).^2).^1.5)./a.dEs');



%% Trust Reagion Reflecive algorythm
function [x, ki_sqrd, dEc] = optimize_q2(a)

options=optimset('Display','off','TolFun',1e-20,'TolX',1e-20);

f = @(x) myfun3(x,a);
x = lsqnonlin(f,a.p,a.p - a.lb,a.p + a.ub, options);

L0 = length(a.sns);
if L0 > 1; np = L0 - 1;
else       np = 1;
end

% Measurement errors
mErr = [0.32 0.31 0.31 NaN 0.18 0.29 0.53 0.16 0.13 0.17 NaN];

dEc = a.k*x(4)*x(3)./((a.x-x(1)).^2+(a.y-x(2)).^2+x(3).^2).^1.5;

ki_sqrd = sum(((dEc' - a.dEs)./mErr(a.sns)).^2)/np;



function F = myfun3(x,a)

F = ((a.dEs' - a.k*x(4)*x(3)./((a.x-x(1)).^2+(a.y-x(2)).^2+x(3).^2).^1.5)./a.dEs');


