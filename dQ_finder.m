function argout = dQ_finder(a)
% Charge tranfer due to electrostatic change

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
    
    a.printOut = 1;
    
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
[dQ N ki dEc] = find_dQ2(a,-100,100,1);

if a.printOut
    fprintf('\nCalculated charge  \t\t= %0.4f C\n',dQ);
    fprintf('Reduced Chi-squired \t= %0.2f\n',ki);
    fprintf('Upper location \t\t\t= [%0.1f\t%0.1f\t%0.1f] m\n',a.p1);
    fprintf('Lower location \t\t\t= [%0.1f\t%0.1f\t%0.1f] m\n',a.p2);
    H = abs(a.p1(3)-a.p2(3))/1000;
    fprintf('Charge Moment \t\t\t= %0.04f C km (Verticle Height = %0.2f km)\n',...
        2*dQ*H,H)
    D = sqrt(sum((a.p1 - a.p2).^2))/1000;
    fprintf('Charge Moment \t\t\t= %0.04f C km (Distance = %0.2f km)\n',...
        2*dQ*D,D)
    
    fprintf('Number of iterations\t= %i\n',N);
    fprintf('_____________________________________________________\n')
    fprintf('Sensor\t\tMeasured\t\tCalculated\t\t%%diff\n')
    fprintf('=====================================================\n')
    for i = 1:length(a.sns)
        fprintf('%s\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
            sen_IDs{a.sns(i)},a.dEs(i),dEc(i),abs(a.dEs(i)-dEc(i))/abs(a.dEs(i)+dEc(i))*200);
    end
end

% outputs
argout.dQ = dQ;
argout.ki = ki;
argout.dEc = dEc;



 function  [dQ N ki dEc] = find_dQ2(a,q1,q2,N)
    
    N = N + 1;
    
    dq = (q2 - q1)/4;
    qs = q1:dq:q2;
    
    L0 = length(a.sns);
    
    % Measurement errors
    mErr = [0.32 0.31 0.31 NaN 0.18 0.29 0.53 0.16 0.13 0.17 NaN];
    
    % number of free parameters
    if L0 > 2; np = L0 - 2;
    else       np = 1;
    end
    
    kis = zeros(1,4);    
    
    for i = 1:4
        
        dEc = -a.k*qs(i)*a.p1(3)./((a.x-a.p1(1)).^2+(a.y-a.p1(2)).^2+a.p1(3).^2).^1.5 + ...
              a.k*qs(i)*a.p2(3)./((a.x-a.p2(1)).^2+(a.y-a.p2(2)).^2+a.p2(3).^2).^1.5;
        
        % D21  = (a.x-a.p2(1)).^2+(a.y-a.p2(2)).^2;
        % D22  = (a.x-a.p2(1)).^2+(a.y-a.p2(2)).^2;
                
        % dEc = -a.k*qs(i).*(2*a.p1(3)^2-D21)./(D21+a.p1(3).^2).^2.5 + ...
        %       a.k*qs(i).*(2*a.p2(3)^2-D22)./(D22+a.p2(3).^2).^2.5;
          
        % kis(i) = sqrt(sum((dEc' - a.dEs).^2)/4);
        kis(i) = sum(((dEc' - a.dEs)./mErr(a.sns)).^2)/np;
    end

    
    [ki indx] = min(kis);  
    
    if N > 100 || (q2 - q1) < 10000*eps
        dQ = (q1+q2)/2;      
        return
    end
    
    [dQ N ki dEc] = find_dQ2(a,qs(indx)-dq,qs(indx)+dq,N);
   