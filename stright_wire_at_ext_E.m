function stright_wire_at_ext_E

% This program is written to find out the charge distribution (and then
% E-field at a point) of a finite stright conductor. The method following
% can be found in "Advanced Engineering Electromagnetics" by Constantine A.
% Balanis in page numbers 670-675.


%% User inputs
L = 300;          % Length of the conductor in meters
N = 101;          % Number of segments of the conductor to consider
V = -0e6;          % Voltage of the conductor
r = 0.005;      % Radius of the conductor
d = 500;          % calculate E-field at d away from the end of the conductor
E = 1e5;        % External electric field

%% Finding the charge distribution.
dy = L/N;
ys = -L/2+dy/2:dy:L/2;

Zmn = zeros(N,N); % For constant V


% Ex = (ys-ys(N/2+1))./((ys-ys(N/2+1)).^2+a^2).^1.5;
% 
% figure
% plot(Ex)
% 
% return

for m = 1:N
    for n = 1:N

        % y' (for numerical intergration)
        %yps = (n-1)*dy:dy/1000:n*dy;
        yps = ys(n)-dy/2:dy/1000:ys(n)+dy/2;

     % For constant V     
        % function evaluation at y' s
        fps = 1./sqrt((ys(m)-yps).^2+r^2);
        
        % Intergration
        Zmn(m,n) = trapz(yps,fps);        
    end
end

Vm = 4*pi*8.85418782e-12.*(V+ys'*E);

% Getting solutions for constatnt V
a = Zmn\Vm;
a = a';

%% Plotting and post calculations

% Plotting charge vs length diagram


figure
[xst yst] = stairs([ys-dy/2 ys(end)+dy/2],[a a(end)]);
plot(yst,xst);
xlabel('Charge density (C/m)')
ylabel('Altitude (m)')

% Total Charge
q = dy*a;
Qtot = sum(q);

title(sprintf('E = %0.1e V/m      V = %0.1e V     Q = %0.2e C', E, V, Qtot))
tools2fig


% E field along the axis
y2s = [-d:dy:-L/2-dy ys  L/2+dy:dy:d];
N1 = length(y2s);
Ex = zeros(1,N1);
for i = 1:N1
    Ex(i) = sum(1/(4*pi*8.85418782e-12).*q.*(y2s(i) - ys)./(r^2+(y2s(i) - ys).^2).^1.5) +E;
end

figure
hold all
plot(Ex,y2s,'ro-')
ylabel('Distance (m)')
xlabel('E (V/m)')
title(sprintf('E = %0.1e V/m      V = %0.1e V     Q = %0.2e C', E, V, Qtot))
plot(xlim,[0 0],'k')
plot([0 0],ylim,'k')
tools2fig
box on

%% Calculating the E field at a point
x0 = -d;
y0 = 0;
Ex = zeros(1,N);
Ey = zeros(1,N);

for i = 1:N
    x = x0 - ys(i);
    [Ex(i) Ey(i)] = E_due_to_a_ring(x,y0,r,q(i));
end

Ex = sum(Ex);
Ey = sum(Ey);


%% Find off axis E-filed of a charged ring
function [Ex Ey] = E_due_to_a_ring(x,y,r,q)

if x == 0 && y == 0
    Ex = 0;
    Ey = 0;
else
    
    % Number of parts to consider in the ring
    N = 10;
    dpi = 2*pi/N;
    pis = 0:dpi:2*pi;
    dq = q/(2*pi)*dpi;
    
    dys = (y-r*sin(pis));
    
    Es = (1/(4*pi*8.85418782e-12))*dq./(x^2+dys.^2+(r*cos(pis)).^2);
    
    Exs = Es.*x./sqrt(x^2+dys.^2);
    Eys = Es.*dys./sqrt(x^2+dys.^2);
    
    Ex = sum(Exs);
    Ey = sum(Eys);
end







