function k_change_charge_moment
% This an attemt to find charge moment of a K change.
tic


sen_set = open('sensor_setting.mat');

x = sen_set.x(1:10);
y = sen_set.y(1:10);
z = sen_set.z(1:10);

% Let's calculate field changes at each sensor location using
% charge Q and in positions H1 and H2

Q0  = 4.0;

H1 = 5500;
H2 = 6200;
x0 = 0;
y0 = 0;
x02 = 100;
y02 = 100;
epsi0 = 8.85418782e-12;



H = (H1+H2)/2;

E = -(1/2/pi/epsi0)*Q0*...
     (H1./(H1^2+(x-x0).^2+(y-y0).^2).^1.5...
    - H2./(H2^2+(x-x02).^2+(y-y02).^2).^1.5);

        
%E = round(E*10)/10;

a.h1 = [4000:100:7500];
a.h2 = a.h1;
a.E = E;
a.sen_set = sen_set;
a.x0 = x0;
a.y0 = y0;

Q = 1:.1:10;
xi_m = size(Q);
H11 = Q-Q;
H22 = H11;

for i = 1:length(Q)
    a.Q = Q(i);
    [H11(i), H22(i), m] = find_H1H2(a);
    
    xi_m(i) = m;
end

% Calculating charge moment
[peakLoc peakMag]= peakfinder(-xi_m);

fprintf('\nQ(C)\tH1(m)\tH2(m)\tdH(m) \tQdH(C.m)\n')
for i=1:length(peakLoc)
    ind = peakLoc(i);
    fprintf('%0.1f\t%i\t%i\t%i\t%.1f\n', ...
        Q(ind),H11(ind),H22(ind),(H22(ind)-H11(ind)),Q(ind)*(H22(ind)-H11(ind)))
end

% Average charge moment
avg_q_mom = mean(Q(peakLoc).*(H22(peakLoc)-H11(peakLoc)));

fprintf('\nAverage chare moment = % 0.1f C.m\n\n', avg_q_mom)


figure
tit = sprintf('dH = %.1fm \ndQ = %0.3fC   \nH1 = %0.1f   \nH2 = %0.1f  \nQ = %0.1f',...
    (a.h1(2)-a.h1(1)),Q(2)-Q(1),H1,H2,Q0);
plot(Q,xi_m)
hold all
plot(Q(peakLoc),xi_m(peakLoc),'or','MarkerFaceColor','r')
ylabel('Minimum \chi^2')
xlabel('Charge (C)')
annotation('textbox',[.7 .3 .1 .1],'String',tit)


figure
plot(Q,H11);
hold all
plot(Q,H22)
legend('H1','H2')
xlabel('Charge (C)')
ylabel('Altitude (m)')
annotation('textbox',[.7 .3 .1 .1],'String',tit)




    
    

% [H1, H2] = find_H1H2(a)
% 
% a.H1 = H1;
% a.H2 = H2;
% %
% find_Q(a)

toc


function [H1,H2,m]= find_H1H2(a)

wbh = waitbar(0,'Please wait....',....
    'Name', 'Charge moment Finder');

x = a.sen_set.x(1:10);
y = a.sen_set.y(1:10);
z = a.sen_set.z(1:10);

x0 = a.x0;
y0 = a.y0;

epsi0 = 8.85418782e-12;

L = length(a.h1);
ki_sqrd = NaN(L,L);

for i = 1:L
    H1 = a.h1(i);
    for j = 1:L;
        
        try
            temp = ((i-1)*L+j)/L/L;
            msg = sprintf('Please wait.... %.1f%% [%0.1fC]',temp*100,a.Q);
            waitbar(temp,wbh,msg)
        catch
            fprintf('User stopped the process!\n')
            return
        end
        
        H2 = a.h2(j);
        
        E_cal = -(1/2/pi/epsi0)*a.Q*(H1./(H1^2+(x-x0).^2+(y-y0).^2).^1.5...
               - H2./(H2^2+(x-x0).^2+(y-y0).^2).^1.5);
           
        % diff      
         ki_sqrd(i,j) = sum((E_cal-a.E).^2);         
    end
end

m = min(min(ki_sqrd));

[i,j] = find(ki_sqrd == m);

H1 = a.h1(i);
H2 = a.h2(j);

delete(wbh) 

function Q = find_Q(a)

x = a.sen_set.x(1:10);
y = a.sen_set.y(1:10);
z = a.sen_set.z(1:10);

x0 = a.x0;
y0 = a.y0;

H1 = a.H1;
H2 = a.H2;

epsi0 = 8.85418782e-12;

   Q = - a.E ./ (1/2/pi/epsi0) ./ ...
                (H1./(H1^2+(x-x0).^2+(y-y0).^2).^1.5...
               - H2./(H2^2+(x-x0).^2+(y-y0).^2).^1.5);

 