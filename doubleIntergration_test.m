function doubleIntergration_test
% clc
% t = 0:1e-6:0.001;
% E = zeros(size(t));
% L = length(t);
% wbh = waitbar(0,'working on...');
% for i=1:L
%    waitbar(i/L,wbh)
%    E(i)=dblquad(@F,0,t(i),0,4000);
% end
% delete(wbh)
% 
% figure
% plot(t,E)
% 
% function Y = F(T,H) 
%     
%   D = 2500;
%   c = 3e8;
%   Ip=10000;
%   L = length(T);
%   Y = size(T);
%   
%   for i = 1:L
%       t = T(i);
%       if t < 0
%           Y(i) = 0;
%       elseif t < 2.5e-6
%           Y(i) = D^2./(c^2*(D^2+H.^2).^1.5)*Ip*t/2.5e-6;
%       elseif t < 25e-6
%           Y(i) = D^2./(c^2*(D^2+H.^2).^1.5)*(Ip-Ip*t/22.5e-6);
%       else
%           Y(i) = 0;
%       end
%   end
%       
%     
 
clc
t = 0:1e-7:0.00001;
E = zeros(size(t));
I = E;
E2 = E;
L = length(t);
wbh = waitbar(0,'working on...');
for i=1:L
   waitbar(i/L,wbh)   
  I(i) = myfun(t(i));
%    if i > 1     
%     E(i) = trapz(t(1:i),I(1:i));
%    end
   
   %E2(i) = integral(@myfun,0,t(i),'AbsTol',1e-15);
   
   %E2(i)=quad(@myfun,0,t(i));
end
 
% add point
ind = find(t <= 2.5e-6);
ind = ind(end);
Itemp = I(ind+1:end);
ttemp = t(ind+1:end);

t(ind+1) = t(ind);
t(ind+2:end+1) = ttemp;
I(ind+2:end+1) = Itemp;

E = cumtrapz(t,I);





delete(wbh)

figure
hold all
plot(t,E*1000000)
plot(t,I)
%plot(t,E2*1000000)

    
 function Y = myfun(T) 
  
  
  Ip=10000;
  L = length(T);
  Y = zeros(size(T));

  for i = 1:L
      t = T(i);
      if t < 0
          % do nothing
      elseif t <= 2.5e-6
          Y(i) = 0.10;
      elseif t <= 5.0e-6
          Y(i) = -0.10;
      else
          % do nothing
      end
  end


    