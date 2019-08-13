function leader_simple
% This function deposite charge along lines and just to see electric fields
% calculated at different distance. This program does not intend to match
% or model real e field recorded. But this will help to understant how
% electric field look at different distances when leader go horizontally.

test_case  = 1;

switch test_case
    case 1
        % Leader going vertically up and then towards sensors
        y = -5000:100:5000;
        x = zeros(size(y));
        %z = [0:100:5000 zeros(1,50)+5000];
        z = zeros(size(x))+5000;
       
        % The position to deposite positive charge
        x0 = 0;
        y0 = -5000;
        z0 = 5000;
        
               
    otherwise
        % do nothing
end

% sensor positions (assuming all sensors are in x axis and y=0 and z=0)
x_s= -50000;
k  = 1; % Coulomb law constant (half the normal value becuse earth's conductivity)
q  = 1;  % charge depsited at each step


E_old = 0;
E_cal = x-x;

for i = 1:length(x)
   E_cal(i) = k.*( -q*z(i)./((x(i)-x_s)^2+y(i)^2+z(i)^2)) + ...
              k.*( q*z0./((x0-x_s)^2+y0^2+z0^2)) + ...
              E_old;
   E_old = E_cal(i);
end

figure(gcf)
hold all
plot(E_cal)