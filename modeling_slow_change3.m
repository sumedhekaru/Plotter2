function modeling_slow_change3
% This is an attempt to model slow change in IBP stage
% In this versionjust use a single dipole growing and changing.

%% User inputs
sns = [1 2 3 6 8 9 10];
bf = 'C:\Users\Sumedhe\Desktop\IBP modeling TCM\20110814_01828\';
fn = '20110814_01828s_809ms_';

p1 = [-17611.0	-42082.3	8684.8];    % Lower location
p2 = [-17873.6	-42599.9	9651.9];    % Upper location
Nh = 100;                               % Number of height steps to consider
t1 = 1828.81376819;
t2 = 1828.82576824;



%% Load data
sen_set = open('sensor_setting.mat');

Ns = length(sns);
% figure;
% hold all

for i = 1:Ns
    sid = sns(i);
    data(i) = open(sprintf('%s%s%s.mat',bf,fn,sen_set.sen_IDs{sid}));
    
    % remove propagation time (roughly)
    tshift = sqrt((sen_set.x(sid) - mean([p1(1),p2(1)]))^2 + ...
                  (sen_set.y(sid) - mean([p1(2),p2(2)]))^2 + ...
                  (sen_set.z(sid) - mean([p1(3),p2(3)]))^2)/ 3.0e8;
    data(i).t = data(i).t - tshift;
    
    %plot(data(sid).t,data(sid).v)
end


% Point array
dx = (p2(1)-p1(1))/Nh;
dy = (p2(2)-p1(2))/Nh;
dz = -(p1(3)-p2(3))/Nh;
dt = (t2 - t1)/Nh;

xs = p1(1):dx:p2(1);
ys = p1(2):dy:p2(2);
zs = p1(3):dz:p2(3);
ts = t1:dt:t2;

if isempty(xs); disp('change sign of dx and re run'); return; end;
if isempty(ys); disp('change sign of dy and re run'); return; end;
if isempty(zs); disp('change sign of dz and re run'); return; end;

a.x = sen_set.x(sns);
a.y = sen_set.y(sns);
a.z = sen_set.z(sns);
a.k = 1 / 2 / pi / 8.85418782e-12;
a.p1 = p1;
a.sns = sns;

wbh = waitbar(1/Nh,'Working on...','name','Slow change modeling');

for i = 2:Nh
    
    try
        waitbar(i/Nh)
    catch
        disp('User stopped modeling slow change')
        return
    end
    
    dEs = nan(1,Ns);
    for j = 1:Ns
        E2 = interp1(data(j).t,data(j).v,ts(i),'pchip');
        E1 = interp1(data(j).t,data(j).v,ts(1),'pchip');
        dEs(j) = E2 - E1;
        
        if i == 2
            data(j).t2(1) = ts(i-1);
            data(j).v2(1) = E1;
            data(j).Q(1) = 0;
            data(j).x(1) = p1(1);
            data(j).y(1) = p1(2);
            data(j).z(1) = p1(3);
        end
    end
    
        
    dEs = dEs(~isnan(dEs));
    a.dEs = dEs;
    a.p2 = [xs(i) ys(i) zs(i)]; 
    [dQ ki dEc] = optimize_q(a);
    a.dQ = dQ;
    a.lb = [5000 5000 5000 100];
    a.ub = [5000 5000 5000 100];
    [P2 ki dEc] = optimize_p2_2(a);
    dQ = P2(4);
    
    for j = 1:Ns
        data(j).t2(i)  = ts(i);
        data(j).v2(i)  = dEc(j) + data(j).v2(1);
        data(j).kis(i) = ki;
        data(j).dQ(i)  = dQ;
        data(j).Q(i)   = dQ;
        data(j).x(i) = P2(1);
        data(j).y(i) = P2(2);
        data(j).z(i) = P2(3);
    end   
    
end

% Display results
for i = 1:Ns
    figure
    plot(data(i).t,data(i).v)
    hold all
    plot(data(i).t2,data(i).v2)
    legend([sen_set.sen_IDs{sns(i)} ':Real'],...
           [sen_set.sen_IDs{sns(i)} ':Modeled'])
end

figure
plot(data(1).t2,data(1).kis)
title('\chi^2')

figure
plot(data(1).t2,data(1).Q)
title('Cumulative Charge')

figure
plot(data(1).t2,data(1).z)
title('Altitude')

figure
plot3(data(1).x,data(1).y,data(1).z)
lims = zlim;
hold all
plot3(data(1).x,data(1).y,zeros(size(data(1).z)) + lims(1))
title('Upper charge location')
daspect([1 1 1])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
box on


delete(wbh)

 function F = myfun(x,a)

F = a.dEs + a.k*x(4)*a.p1(3)./((a.x-a.p1(1)).^2+(a.y-a.p1(2)).^2+a.p1(3).^2).^1.5 - ...
    a.k*x(4)*x(3)./((a.x-x(1)).^2+(a.y-x(2)).^2+x(3).^2).^1.5;

function F = myfun2(x,a)

F = a.dEs + a.k*x*a.p1(3)./((a.x-a.p1(1)).^2+(a.y-a.p1(2)).^2+a.p1(3).^2).^1.5 - ...
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

ki_sqrd = sum(((dEc - a.dEs)./mErr(a.sns)).^2)/np;

P2 = x;

function [P2 ki_sqrd dEc] = optimize_p2_2(a)
% Optimize the lower point using Trust region reflective

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

ki_sqrd = sum(((dEc - a.dEs)./mErr(a.sns)).^2)/np;

P2 = x;
