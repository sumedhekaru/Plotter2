function stright_wire_of_constant_V

% This program is written to find out the charge distribution (and then
% E-field at a point) of a finite stright conductor. The method following
% can be found in "Advanced Engineering Electromagnetics" by Constantine A.
% Balanis in page numbers 670-675.


%% User inputs
L = 1;          % Length of the conductor in meters
N = 5;          % Number of segments of the conductor to consider
V = 1;          % Voltage of the conductor
a = 0.001;      % Radius of the conductor
d = 1;          % calculate E-field at d away from the end of the conductor


%% Finding the charge distribution.
dy = L/N;
ys = dy/2:dy:L;

Zmn = zeros(N,N);

for m = 1:N
    for n = 1:N
        % y' (for numerical intergration)
        yps = (n-1)*dy:dy/1000:n*dy;
        
        % function evaluation at y' s
        fps = 1./sqrt((ys(m)-yps).^2+a^2);
        
        % Intergration
        Zmn(m,n) = trapz(yps,fps);        
    end
end

Vm = zeros(N,1) + 4*pi*8.85418782e-12*V;

a = inv(Zmn)*Vm;
a = a';

%% Plotting and post calculations

% Plotting charge vs length diagram
figure
stairs([ys-dy/2 ys(end)+dy/2],[a a(end)])
xlabel('Length (m)')
ylabel('Charge density (C/m)')

% Total Charge
q = dy*a;
Qtot = sum(q);

% E-field at distance d;
E = sum(1/(4*pi*8.85418782e-12).*q./((L+d - ys).^2));

% print info to the title
title(sprintf('Total charge = %0.2e C    E. field at %0.1f m away = %0.2e V/m',Qtot,d,E))






