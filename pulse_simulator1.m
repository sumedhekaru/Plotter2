function data = pulse_simulator1

%% User Inputs
a=struct( ...
    'maxA' , {50000}       , ...
    't1' , {2.5000e-006} , ...
    't2' , {4.500e-006} , ...
    'v' , {1e8}   , ...
    'H1' , {5000}        , ...
    'H2' , {7000}        , ...
    'x0' , {0}, ...
    'y0' , {0} , ...
    'T1' , {0.1500e-004}           , ...
    'T2' , {0.550e-004} , ...
    'lamda' , {200}         , ...
    't_step' , {1.0000e-007} , ...
    'x' , {1000}       , ...
    'y' , {0}         ...
    );

%% Other parameters
a.dh = 10;   % Height increments


%%

% lamda serves a purpose similar to the characteristic length in RB
a.k       = (a.t2-a.t1)/a.t1;
a.alpha   = 10/a.t2;

eps0    = 8.8542e-12;               % Permitivity of free space
c       = 299792458/1.0003;         % Speed of light in air


% get the distance from the pulse to sensors
a.x = a.x - a.x0;
a.y = a.y - a.y0;
D=sqrt(a.x^2+a.y^2);

% Counter
i=0;


ln=length(a.T1:a.t_step:a.T2);

Z = (a.H1:a.dh:a.H2)';
E_ind=zeros(1,ln);
E_rad=E_ind;
E_stat=E_ind;

wbh = waitbar(0,'Please wait...','Name','NBP Simulator');

for t=a.T1:a.t_step:a.T2
    
    i=i+1;
    try
        waitbar(i/ln,wbh,sprintf('Please Wait... %.0f%%',i*100/ln))
    catch
        disp('NBP Program Stopped by User!')
        return
    end
    
    % Current and derivation of current
    [ii,didt,intI]=current(Z,t-(Z-a.H1)/a.v-((D^2+Z.^2).^0.5)/c,a);
   
    
    % 1/R radiation term
    %Y = (-1/2/pi/eps0)*D^2./(c^2*(D^2+Z.^2).^1.5).*didt;
    
    % 1/R^1.3 radiation term
    Y = (-1/2/pi/eps0)*D^2./(c^2*(D^2+Z.^2).^1.5).*didt;    
    E_rad(i) = trapz(Z,Y);
%     Y = didt;
%     if max(Y) > 1e11;
%         figure
%         plot(Y)
%         return
%     end
%     E_rad(i) = trapz(Z,Y);
    
    % Induction Term
    Y = (1/2/pi/eps0)*((2*Z.^2-D^2)./(c*((Z.^2+D^2).^2))).*ii;
    E_ind(i) = trapz(Z,Y);
        
%     % Electro Static Term
    Y=(1/2/pi/eps0)*(2*Z.^2 - D^2)./((D^2+Z.^2).^2.5).*int_i(Z,t-(Z-a.H1)/a.v-(D^2+Z.^2).^0.5/c,a);
    % Y=(1/2/pi/eps0)*(2*Z.^2 - D^2)./((D^2+Z.^2).^2.5).*intI;
     E_stat(i) = trapz(Z,Y);
    
    % Dipole moment
%     if i==ln
%         Y=int_i(Z,t-(Z-a.H1)/a.v-(D^2+Z.^2).^0.5/c);
%         dipole=trapz(Z,Y);
%     end
end

t=(a.T1:a.t_step:a.T2);
delete(wbh)
figure
plot(t,E_rad)
 hold all

plot(t,E_ind)
plot(t,E_stat)
plot(t,E_rad+E_ind+E_stat)

legend('rad','ind','stat','tot')
tools2fig


function [ii didt intI]=current(z,t,a)

    ii=zeros(size(z));
    didt=zeros(size(z));
    
    ind = sum(t <= a.t1);

    % Current
    %ii(1:ind) = a.maxA*exp(-(z(1:ind)-a.H1)/a.lamda).* ...
    %    exp(-(a.alpha*(t(1:ind)-a.t1)).^2);
    
    %ii(ind+1:end) = a.maxA*exp(-(z(ind+1:end)-a.H1)/a.lamda).* ...
    %    exp(-(a.alpha/a.k*(t(ind+1:end)-a.t1)).^2);
    
    % Derivative of current
%     didt(1:ind) = a.maxA*exp(-(z(1:ind)-a.H1)/a.lamda).* ...
%                 (-2*a.alpha^2.*(t(1:ind)-a.t1)).* ....
%                 exp(-(a.alpha*(t(1:ind)-a.t1)).^2);
%     
%     didt(ind+1:end) = a.maxA*exp(-(z(ind+1:end)-a.H1)/a.lamda).* ...
%                 (-2*(a.alpha/a.k)^2 .*(t(ind+1:end)-a.t1)).* ....
%                 exp(-(a.alpha/a.k *(t(ind+1:end)-a.t1)).^2);

%             % MTLE
%             ii(1:ind)  = a.maxA*exp(-(z(1:ind)-a.H1)/a.lamda).*exp(-a.alpha^2*(t(1:ind)-a.t1).^2);
%             didt(1:ind)= exp(-(z(1:ind)-a.H1)/a.lamda).*(-2)*a.maxA*a.alpha.^2.*(t(1:ind)-a.t1).*exp(-a.alpha^2*(t(1:ind)-a.t1).^2);
%             
%             % MTLE
%             ii(ind+1:end)=a.maxA*exp(-(z(ind+1:end)-a.H1)/a.lamda).*exp(-(a.alpha/a.k)^2*(t(ind+1:end)-a.t1).^2);
%             didt(ind+1:end)=exp(-(z(ind+1:end)-a.H1)/a.lamda).*(-2)*a.maxA*(a.alpha/a.k)^2.*(t(ind+1:end)-a.t1).*exp(-(a.alpha/a.k)^2*(t(ind+1:end)-a.t1).^2);

           

    for j=1:length(t)
        
        if t(j) <= a.t1
            % MTLEI
            %                 ii(j)=a.maxA*exp((z(j)-a.H1)/a.lamda)*exp(-a.alpha^2*(t(j)-a.t1)^2);
            %                 didt(j)=exp((z(j)-a.H1)/a.lamda)*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);

                             % TL
                            % ii(j)=a.maxA*exp(-a.alpha^2*(t(j)-a.t1)^2);
                            % didt(j)=(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);

            % MTLL
            %                 ii(j)=a.maxA*(1 - (z(j)-a.H1)/a.H1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
            %                 didt(j)=(1 - (z(j)-a.H1)/a.H1)*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);

            % MTLE
            ii(j)  = a.maxA*exp(-(z(j)-a.H1)/a.lamda)*exp(-a.alpha^2*(t(j)-a.t1)^2);
            didt(j)= exp(-(z(j)-a.H1)/a.lamda)*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
            
            


        else
            %                 % MTLEI
            %                ii(j)=a.maxA*exp((z(j)-a.H1)/a.lamda)*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            %                 didt(j)=exp((z(j)-a.H1)/a.lamda)*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            %
            %                 % TL
                           % ii(j)=a.maxA*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                            % didt(j)=(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            %
            %                 % MTLL
            %                 ii(j)=a.maxA*(1 - (z(j)-a.H1)/a.H1)*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            %                 didt(j)=(1 - (z(j)-a.H1)/a.H1)*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);

            % MTLE
            ii(j)=a.maxA*exp(-(z(j)-a.H1)/a.lamda)*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            didt(j)=exp(-(z(j)-a.H1)/a.lamda)*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);

        end
        
    end
    
    if length(t) > 1
        intI = cumtrapz(t,ii);
    end
    
    
function inti=int_i(z,t,a)

    tau=t-a.t_step:a.t_step/10:t;

    inti=zeros(size(z));

    for j=1:length(z)
        % t=t-(z-a.H1)/a.v;
        cr=zeros(size(tau));
        for kk=1:length(tau)
           cr(kk)= current(z(j),tau(kk),a);
        end

        %length(tau)
        %length(cr)

        inti(j)=trapz(tau,cr);

    end



