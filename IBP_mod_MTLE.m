function [t,E_stat,E_ind,E_rad,P] = IBP_mod_MTLE(maxA,t1,t2,v,H1,H2,x0,y0,T1,...
       T2,t_step,lamda,x,y,dh,alpha)
% User Inputs
if nargin == 0
    a=struct( ...
        'maxA' , {-420}        , ...
        't1' , {10e-006}    , ...
        't2' , {32.0e-006}     , ...
        'v' ,  {4.6e7}             , ...
        'H1' , {7600}           , ...
        'H2' , {8400}          , ...
        'x0' , {0}              , ...
        'y0' , {0}              , ...
        'T1' , {-0.100e-004}    , ...
        'T2' , {3.5e-004}     , ...
        't_step' , {1.0000e-006}, ...
        'lamda', {113}, ...
        'alpha', {10}, ...
        'x' , {28000}            , ...
        'y' , {0}               ,...
        'dh', {10}               ...
        );
else
    a.maxA = maxA;
    a.t1 = t1;
    a.t2 = t2;
    a.v = v;
    a.H1 = H1;
    a.H2 = H2;
    a.x0 = x0;
    a.y0 = y0;
    a.T1 = T1;
    a.T2 = T2;
    a.t_step = t_step;
    a.lamda = lamda;
    a.x = x;
    a.y = y;
    a. dh = dh;
    a.alpha = alpha;
end


a.dh = a.v*a.t_step;

%% Start modeling

% lamda serves a purpose similar to the characteristic length in RB
% dH = abs(a.H1 - a.H2);
% a.lamda = dH/log(dH);


a.k       = (a.t2-a.t1)/a.t1;
% a.alpha   = 10/a.t2;

eps0    = 8.8542e-12;               % Permitivity of free space
c       = 299792458/1.0003;         % Speed of light in air


% get the horizontal distance from the pulse to sensors
a.x = a.x - a.x0;
a.y = a.y - a.y0;
D=sqrt(a.x^2+a.y^2);

% Counter
i=0;
ln=length(a.T1:a.t_step:a.T2);

if a.H2 > a.H1
    Z = (a.H1:a.dh:a.H2)';
else
    Z = (a.H1:-a.dh:a.H2)';
    a.maxA = -a.maxA;
end


%wbh = waitbar(0,'Please wait...','Name','IBP Modeler');

intgr=0;

for t=a.T1:a.t_step:a.T2
    
    i=i+1;
%     try
%         waitbar(i/ln,wbh,sprintf('Please Wait... %.0f%%',i*100/ln))
%     catch
%         disp('NBP Program Stopped by User!')
%         return
%     end
    
    
    % Current and derivation of current
    [ii,didt]=current(Z,t-abs(Z-a.H1)/a.v-((D^2+Z.^2).^0.5)/c,a);
   
    
   
    % Radiation term
    Y = (-1/2/pi/eps0)*D^2./(c^2*(D^2+Z.^2).^1.5).*didt;
    E_rad(i) = trapz(Z,Y);
    
    % Induction Term
    Y = (1/2/pi/eps0)*((2*Z.^2-D^2)./(c*((Z.^2+D^2).^2))).*ii;
    E_ind(i) = trapz(Z,Y);
    
    % Electro Static Term
    %Y=(1/2/pi/eps0)*(2*Z.^2 - D^2)./((D^2+Z.^2).^2.5).*intI;
    intgr = int_i(Z,t-abs(Z-a.H1)/a.v-(D^2+Z.^2).^0.5/c,a,intgr);
    Y=(1/2/pi/eps0)*(2*Z.^2 - D^2)./((D^2+Z.^2).^2.5).*intgr;
    E_stat(i) = trapz(Z,Y);
    
    %     figure(500)
    %     plot(Y)
    
    
end
% Charge moment
P = trapz(Z,intgr);


t=(a.T1:a.t_step:a.T2);
% figure
% hold all
% tools2fig
% plot(t,E_rad)
% plot(t,E_ind)
% plot(t,E_stat)


%delete(wbh)


function [ii, didt]=current(z,t,a)

ii=zeros(size(z));
didt=zeros(size(z));

for j=1:length(t)
    
    if t(j) <= a.t1
        % MTLEI
                         %ii(j)=a.maxA*exp((z(j)-a.H1)/a.lamda)*exp(-a.alpha^2*(t(j)-a.t1)^2);
                         %didt(j)=exp((z(j)-a.H1)/a.lamda)*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
        
        % TL
        % ii(j)=a.maxA*exp(-a.alpha^2*(t(j)-a.t1)^2);
        % didt(j)=(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
        
        % MTLL
        %                 ii(j)=a.maxA*(1 - (z(j)-a.H1)/a.H1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
        %                 didt(j)=(1 - (z(j)-a.H1)/a.H1)*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
        
        % MTLE
        
        ii(j)  = a.maxA*exp(-abs(z(j)-a.H1)/a.lamda)*exp(-a.alpha^2*(t(j)-a.t1)^2);
        didt(j)= exp(-abs(z(j)-a.H1)/a.lamda)*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
        
        
        
        
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
        ii(j)=a.maxA*exp(-abs(z(j)-a.H1)/a.lamda)*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        didt(j)=exp(-abs(z(j)-a.H1)/a.lamda)*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        
    end
    
end

       
function intgr=int_i(z,t,a,intgr)

        tau=t-a.t_step:a.t_step/1:t;        
        
        inti=zeros(size(z));
        
        for j=1:length(z)

            cr=zeros(size(tau));
            for kk=1:length(tau)
                cr(kk)= current(z(j),tau(kk),a);
            end
           
            inti(j)=trapz(tau,cr);
            
        end

        intgr=inti+intgr;
 

       
    
