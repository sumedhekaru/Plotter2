% Testing whether Kumaraswamy modal is automatically satisfy the charge
% conservation.
function test2

lamda1 = 1.5:0.5:5;
lamda2 = 1.5:0.5:5;

L1 = length(lamda1);
L2 = length(lamda2);



d.maxA = 1000;
d.t1 = 10e-6;
d.t2 = 25e-6;
d.v  = 1.8e8;
d.dt = 1e-6;
d.H1 = 5000;
d.H2 = 6000;
d.AlpInd = 8.7;
d.modal = 'MTLK';
d.dh = d.v*d.dt;


q1s = nan(L1,L2);
q2s = q1s;
qts = q1s;

for i = 1:L1
    d.lamda1 = lamda1(i);
    for j = 1:L2
        d.lamda2 = lamda2(j);
        d.Hm = d.H1 + (d.H2-d.H1)*((d.lamda1 - 1)/(d.lamda1*d.lamda2 - 1))^(1/d.lamda1);
        tic
        [q1s(i,j) q2s(i,j) qts(i,j)] = do_calculations(d);
        toc
    end
end

figure
hold all
plot(-reshape(q1s,1,L1*L2),'ro-','markerfacecolor','r')
plot(reshape(q2s,1,L1*L2),'bo-','markerfacecolor','b')
plot(-reshape(qts,1,L1*L2),'ko-','markerfacecolor','k')
legend('Q_+','-(Q_-)','Q_t')


function [q1 q2 qt] = do_calculations(d)
        
    Z = d.H1:d.dh/10:d.H2;
    T = 0:d.dt/10:d.t2*2;

    L1 = length(T);
    L2 = length(Z);
    
    iii = nan(L1,L2);
    
    i = 0;
    
    wbh = waitbar(0,'Plase wait...','name','\rho_L finder');
    
    
    for t=T
        
        i=i+1;
        
        ii = current(Z,t-(Z-d.H1)/d.v,d);
        
        iii(i,:) = ii;
        try
            waitbar(i/L1,wbh);
        catch
            disp('User stopped processing')
            return
        end
    end
    
    delete(wbh)
    
    Q_passed = zeros(1,L2);
    
      
    for i=1:L2
        Q_passed(i) = trapz(T,iii(:,i));
    end
    
    
    rho = -diff(Q_passed)/(d.dh);
    
    Q_dipo = -diff(Q_passed);
    

   
    %if strcmp(d.modal, 'MTLK')
        
        %Calculate linier charge density from Hm to H2
        [mm ind_m] = min(abs(Z - d.Hm));
        
        q1 = sum(Q_dipo(1:ind_m -1));
        q2 = sum(Q_dipo(ind_m:end));
        qt = q1 + q2;
        %lcd2 = totC2/(d.H2 - d.Hm);
        
        
        %fprintf('Total & Linear charge density of upper part of MTLK\n\t%0.2e\t\t%0.2e\n',totC2,lcd2)
    %end
    
%     tit_str = sprintf ('a = %0.1f   b = %0.1f   +q = %.3e  -q = %0.3e   Q_t = %0.3e', ...
%                         d.lamda1,d.lamda2,totC1,totC2,totC1+totC2);
%     figure
%     plot(Q_dipo,Z(2:end))
%     title(tit_str)
%     

function [ii, didt]=current(z,t,a)


ii      = zeros(size(z));
didt    = zeros(size(z));
a.alpha = a.AlpInd/a.t2;
a.k     = (a.t2-a.t1)/a.t1;

for j=1:length(t)
    
    if t(j) <= a.t1
        
        switch a.modal
            case 'MTLEI'
                ii(j)=a.maxA*exp((z(j)-a.H1)/a.lamda1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
                didt(j)=exp((z(j)-a.H1)/a.lamda1)*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
            case 'TL'
                ii(j)=a.maxA*exp(-a.alpha^2*(t(j)-a.t1)^2);
                didt(j)=(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
            case 'MTLL'
                ii(j)=a.maxA*(1 - (z(j)-a.H1)/(a.H2-a.H1))*exp(-a.alpha^2*(t(j)-a.t1)^2);
                didt(j)=(1 - (z(j)-a.H1)/(a.H2-a.H1))*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
            case 'MTLE'
                ii(j)  = a.maxA*exp(-(z(j)-a.H1)/a.lamda1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
                didt(j)= exp(-(z(j)-a.H1)/a.lamda1)*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
                
            case 'MTLK'
                ii(j)=a.maxA*((z(j)-a.H1)/(a.H2-a.H1)).^(a.lamda1-1).*((1-((z(j)-a.H1)/(a.H2-a.H1)).^a.lamda1).^(a.lamda2-1))*exp(-a.alpha^2*(t(j)-a.t1)^2);
                didt(j)=((z(j)-a.H1)/(a.H2-a.H1)).^(a.lamda1-1).*((1-((z(j)-a.H1)/(a.H2-a.H1)).^a.lamda1).^(a.lamda2-1))*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
                
               % ii(j)=a.maxA*get_power(((z(j)-a.H1)/(a.H2-a.H1)),(a.lamda1-1)).*get_power(get_power(1-((z(j)-a.H1)/(a.H2-a.H1)),a.lamda1),(a.lamda2-1))*exp(-a.alpha^2*(t(j)-a.t1)^2);
               % didt(j)=get_power(((z(j)-a.H1)/(a.H2-a.H1)),(a.lamda1-1)).*get_power(get_power(1-((z(j)-a.H1)/(a.H2-a.H1)),a.lamda1),(a.lamda2-1))*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);

        end
        
    else
        switch a.modal
            case 'MTLEI'
                ii(j)=a.maxA*exp((z(j)-a.H1)/a.lamda1)*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                didt(j)=exp((z(j)-a.H1)/a.lamda1)*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            case 'TL'
                ii(j)=a.maxA*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                didt(j)=(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            case 'MTLL'
                ii(j)=a.maxA*(1 - (z(j)-a.H1)/(a.H2-a.H1))*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                didt(j)=(1 - (z(j)-a.H1)/(a.H2-a.H1))*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            case 'MTLE'
                ii(j)=a.maxA*exp(-(z(j)-a.H1)/a.lamda1)*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                didt(j)=exp(-(z(j)-a.H1)/a.lamda1)*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
            case 'MTLK'
                ii(j)=a.maxA*((z(j)-a.H1)/(a.H2-a.H1)).^(a.lamda1-1).*((1-((z(j)-a.H1)/(a.H2-a.H1)).^a.lamda1).^(a.lamda2-1))*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                didt(j)=((z(j)-a.H1)/(a.H2-a.H1)).^(a.lamda1-1).*((1-((z(j)-a.H1)/(a.H2-a.H1)).^a.lamda1).^(a.lamda2-1))*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        
                %ii(j)=a.maxA*get_power(((z(j)-a.H1)/(a.H2-a.H1)),(a.lamda1-1)).*get_power(get_power(1-((z(j)-a.H1)/(a.H2-a.H1)),a.lamda1),(a.lamda2-1))*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                %didt(j)=get_power(((z(j)-a.H1)/(a.H2-a.H1)),(a.lamda1-1)).*get_power(get_power(1-((z(j)-a.H1)/(a.H2-a.H1)),a.lamda1),(a.lamda2-1))*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        end
    end    
end