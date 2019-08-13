function leaderModelingPostProcessing2
% F18 1st RS
fn = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F02\20140509-K24\\ModelDataOut_.mat';
clc

d = open(fn)
L = length(d.tuni);
b = d.b;
sen_set = open('sensor_setting.mat');

eps0 = 8.85418782e-12;

[PATH,NAME] = fileparts(fn);


% %figure
% xl = xlim;
% hold all
% [AX,H1,H2]=plotyy(nan,nan,nan,nan);
% hold(AX(2), 'on')
% %figure
% plot(AX(2),d.tuni(2:end),cumtrapz(d.Qs))
% set(AX,'xlim',xl)
% linkaxes(AX,'x')
% set(AX(1),'ylimmode','auto')
% % mmm=ceil(max([y1;y2])/1000)*1000;
% 
% 
% set(AX(2),'xtick',[])
% set(AX, 'YColor', [0 0 0])


d

%% FieldChange should produced from RS
sns = [1 2 3 5 6 7 8 9 10 11];
% F18 measured RS field changes
%MRS = [-67.9219 -73.76 -101.38 NaN NaN -154.87 NaN NaN NaN NaN];

% for F02
MRS = [-117.7358 -126.1456 -216.9811 NaN NaN -436.6577 NaN  NaN -2.1397 NaN];

fprintf('\nField changes at each sensor\n')
fprintf('Sensor\tMeasured\tCalculated\t%%Diff\n')

%1596
factor = 1.00;
lol = 1;
ul  = 1596;


x = b.x(:,1);
y = b.y(:,1);
z = b.z(:,1);
q = b.qs(:,1);

data = [x  y z q];
data = sortrows(data,3);

N = length(q);

% It is beleive that all return stroke does is nutralizing charges on the
% leader. We will check that here. We will make a dipole by placing
% negative charge on the ground (doesn't make any vertical E change)
% and equal and positive charge on each PBHSV
% locations. See what is the field change that produce.
E1 = MRS - MRS;

for i = sns
    E1(i) = sum(1/(2*pi*eps0)* q.* z ./ ...
        ((x-sen_set.x(i)).^2+(y-sen_set.y(i)).^2+z.^2).^1.5);
end
E1
% Now we have matched E1 amount of the E field, now we need to match the
% rest.
E_rest = MRS -E1


% We are going to match E_rest by placing positive charges on the leader
% path assuming it is in a constant E field

return


%% Finding the charge distribution.

a.data = data(:,1:3);
a.sen_set = sen_set;
a.Em = E_rest;
a.sns = sns;
a.E = -1.5e5;        % External electric field

x0 = 0.01;

% Add PBFA data
a.data = [a.data ; b.pbfa(:,2:4)];
a.data = unique(a.data,'rows');
data = a.data;

x = data(:,1);
y = data(:,2);
z = data(:,3);
N = length(x);


% Potential data
% Voltage data
vData = open('Stolz_and_Marsh2008_4th_baloon_flight_data.mat');
%figure
%hold all
%plot(vData.V,vData.z)

vVals = nan(1,N);
for i = 1:N
    ind = nnz(vData.z < z(i)/1000);
    if ind == 0
        vVals(i) = 0;
    else
        vVals(i) = vData.V(ind) + (vData.V(ind+1)-vData.V(ind))*(z(i)/1000-vData.z(ind))/(vData.z(ind+1)-vData.z(ind));
    end
end

%plot(vVals,z/1000,'ro')
a.vVals = -vVals'*1e6*4;

%% Levenberg-Marquardt
%options=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.01},'TolFun',1e-5,'TolX',1e-5);
%f = @(x) myfun(x,a);
%x = fsolve(f,x0,options);

%% Trust Region reflective

options=optimset('Display','iter','Algorithm','trust-region-reflective','TolFun',1e-10,'TolX',1e-10);
f = @(x) myfun(x,a);
r = lsqnonlin(f,x0,.01,3,options);
%x = fsolve(f,x0,options);

r = abs(r); % Radius of the return stroke channel
%r = 0.24;



Zmn = zeros(N,N); % For constant V

for m = 1:N
    for n = 1:N
        
        % Distance to self charge and its image
        if n == m
            Zmn(m,n) = 1/r - ...
                1/sqrt((x(m)-x(n))^2+(y(m)-y(n))^2+(z(m)+z(n)-r)^2);
            
            % Fix -Inf problem when z(n) == 0
            if z(n) == 0
                Zmn(m,n) = 1/r;
            end
            
            % Distances due to other charges and its image
        else
            Zmn(m,n) = 1/sqrt((x(m)-x(n))^2+(y(m)-y(n))^2+(z(m)-z(n))^2) - ...
                1/sqrt((x(m)-x(n))^2+(y(m)-y(n))^2+(z(m)+z(n))^2);
        end
    end
end


% Follow potential curve from vData
Vm = 4*pi*8.85418782e-12.*(a.vVals);

% Constant E
%Vm = 4*pi*8.85418782e-12.*(z*-a.E);

% Getting solutions for constatnt V
q = Zmn\Vm;


figure
scatter3(a.data(:,1),a.data(:,2),a.data(:,3),5,q,'fill')
colorbar
daspect([1 1 1])
view([pi/2 0 0])
box on
title(sprintf('External E = %0.2e V/m Q_{tot} = %0.3f C   r = %0.1f cm',a.E,sum(q),r*100))



% Now we have the charge distribution, calculate the field change that it
% produce at each station.


for i = sns
    E2(i) = -sum(1/(2*pi*eps0)* q.* z ./ ...
        ((x-sen_set.x(i)).^2+(y-sen_set.y(i)).^2+z.^2).^1.5) ;
end


% Print results

for i = sns
    fprintf('%1.1i\t\t%8.1f\t\t%8.1f\t\t%8.2f\n',i,a.Em(i),E2(i),(a.Em(i)-E2(i))/a.Em(i)*100)
end

figure
hold all
plot(vData.V*1e6,vData.z*1000)
plot(-a.vVals,z,'ro')
plot(z*a.E,z)

legend('Original','Interpolated','Previous')

% % Let's calculate the potential at a point to see it is what we should be
% % getting.
% for k = 1:N    
%     vk = 0;
%     for i = 1:N
%         if k == i
%             if z(i) == 0
%                 vk = vk + 1/(4*pi*8.85418782e-12)*q(i)/r;
%             else
%                 vk = vk + 1/(4*pi*8.85418782e-12)*q(i)*(1/r - ...
%                 1/sqrt((x(m)-x(n))^2+(y(m)-y(n))^2+(z(m)+z(n)-r)^2));
%             end
%         else
%             vk = vk + 1/(4*pi*8.85418782e-12)*q(i)*(1/sqrt((x(k)-x(i))^2+(y(k)-y(i))^2+(z(k)-z(i))^2) - ...
%                 1/sqrt((x(k)-x(i))^2+(y(k)-y(i))^2+(z(k)+z(i))^2));
%         end
%     end
%    fprintf('%0.5e\t\t%0.5e\n',vk,z(k)*a.E)
% end
%     
%     

    



function F = myfun(r,a)

data = a.data;
sen_set= a.sen_set;
data = data(:,1:3);
sns = a.sns;

eps0 = 8.85418782e-12;

data = unique(data,'rows');
x = data(:,1);
y = data(:,2);
z = data(:,3);
N = length(x);

Zmn = zeros(N,N); % For constant V

for m = 1:N
    for n = 1:N
        
        % Distance to self charge and its image
        if n == m
            Zmn(m,n) = 1/r - ...
                1/sqrt((x(m)-x(n))^2+(y(m)-y(n))^2+(z(m)+z(n)-r)^2);
            
            % Fix -Inf problem when z(n) == 0
            if z(n) == 0
                Zmn(m,n) = 1/r;
            end
            
            % Distances due to other charges and its image
        else
            Zmn(m,n) = 1/sqrt((x(m)-x(n))^2+(y(m)-y(n))^2+(z(m)-z(n))^2) - ...
                1/sqrt((x(m)-x(n))^2+(y(m)-y(n))^2+(z(m)+z(n))^2);
        end
    end
end

Vm = 4*pi*8.85418782e-12.*(a.vVals);

%Vm = 4*pi*8.85418782e-12.*(z*-a.E);

% Getting solutions for constatnt V
q = Zmn\Vm;

E2 = nan(1,10);
% Now we have the charge distribution, calculate the field change that it
% produce at each station.
for i = sns
    E2(i) = -sum(1/(2*pi*eps0)* q.* z ./ ...
        ((x-sen_set.x(i)).^2+(y-sen_set.y(i)).^2+z.^2).^1.5) ;
end

F = nanmean(((E2-a.Em)).^2);



