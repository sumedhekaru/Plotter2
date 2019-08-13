function MTL_current_distribution

%% Altitude variations of peak currents
z = 5000:.1:6000;
Ip = 1;
lamda1 = 200;
lamda2 = 200;

% MTLL
I = Ip*(1 -(z-z(1))/(z(end) - z(1)));
I = I/max(I);

figure; hold all; box on;
plot(I,z/1000,'LineWidth',2)
ylabel('Altitude (km)')
xlabel('Normalized Peak Current')


% MTLE
I = Ip*exp(-(z-z(1))/lamda1);
I = I/max(I);
plot(I,z/1000,'LineWidth',2)

% MTLEI
I = Ip*exp((z-z(1))/lamda1);
I = I/max(I);
plot(I,z/1000,'LineWidth',2)

% 
% % MTLER
% Hm = 300;
% I = (z-z(1))/Hm^2.*exp(-(z-z(1)).^2/2/Hm^2);
% I = I/max(I);
% plot(z/1000,I,'LineWidth',2)


box on

% % MTLAG
% Hm = 5300;
% H2 = 8000;
% alpha = 0.006;
% ind1 = sum(z < Hm);
% I(1:ind1) = Ip*exp(-(alpha*(z(1:ind1)-Hm)).^2);
% I(ind1+1:end) = Ip*exp(-(alpha*(z(ind1+1:end)-Hm)*Hm/(H2-Hm)).^2);
% I = I/max(I);
% plot(z/1000,I,'LineWidth',2)

% % MTLARice
% nu = 400;
% sig =100;
% H1 = z(1);
% z = z-H1;
% % I = ricepdf(z, nu, sig);
% I = evpdf(z, nu, sig);
% I = I/max(I);
% z = z+ H1;
% plot(z/1000,I,'LineWidth',2)
% [Hm ind] = max(I);
% Hm = z(ind)

% Kumaraswami
%figure(100); hold all
alp = 1.5;
b = 5;
x = (z-z(1))/range(z);
I = x.^(alp-1).*((1-x.^alp).^(b-1));
I = I/max(I);
plot(I,z/1000,'LineWidth',2)
%[Hm ind] = max(I);
%Hm = z(ind);
%mode = z(1)+((alp-1)/(alp*b - 1))^(1/alp)*range(z);


legend({'MTLL'
        'MTLE (\lambda = 200 m)'
        'MTLEI (\lambda = 200 m)'
        %'MTLR (H_m = 5200 m)'
        %'MTLAG'
        %'MTLR'
        'MTLK (\alpha = 1.5, \beta = 5)'})
    
    
%% Plotting 3d currents
a.t_step = 0.1e-6;
a.T1 = 0;
a.T2 = 60e-6;
a.v = 10e7;
Z = z(1):.1:z(end);
a.H1 = z(1);
a.H2 = z(end);
a.D = 0;
a.c = 3e8;
a.t1 = 8e-6;
a.t2 = 40e-6;
a.maxA = 10000;
a.lamda = 200;
a.alpha = 10/a.t2;
a.k       = (a.t2-a.t1)/a.t1;
a.a = 1.5;
a.b = 5;
a.dh = 10;

i = 0;


T = a.T1:a.t_step:a.T2;

iiii = nan(length(T),length(Z));

for t=T
    
    i=i+1;
%     try
%         waitbar(i/ln,wbh,sprintf('Please Wait... %.0f%%',i*100/ln))
%     catch
%         disp('NBP Program Stopped by User!')
%         return
%     end
    
    
    % Current and derivation of current
    %ii = current(Z,t-(Z-a.H1)/a.v-((D^2+Z.^2).^0.5)/c,a);
    
    ii = current(Z,t-(Z-a.H1)/a.v,a);
    
    iiii(i,:) = ii;
    %figure(100)
    %plot(ii)
    
    
end


Q_passed = zeros(1,length(Z));
for i=1:length(Z)
    Q_passed(i) = trapz(T,iiii(:,i));
end

figure
plot(Z(1:end-1),-diff(Q_passed)/(Z(2)-Z(1))*1e3)

iiii = iiii/max(iiii(:));
Z = Z/1000;
T = T*1e6;
figure
H = surf(Z,T,iiii);
set(H, 'linestyle', 'none');
shading interp;
xlim([min(Z) max(Z)])
ylim([min(T) max(T)])
xlabel('Altitude (km)')
ylabel('Time (\mus)')
zlabel('Normalized current')
tools2fig
%colormap gray





function [ii, didt]=current(z,t,a)

ii=zeros(size(z));
didt=zeros(size(z));

for j=1:length(t)
    
    if t(j) <= a.t1
        % MTLEI
          %               ii(j)=a.maxA*exp((z(j)-a.H1)/a.lamda)*exp(-a.alpha^2*(t(j)-a.t1)^2);
          %              didt(j)=exp((z(j)-a.H1)/a.lamda)*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
        
        % TL
        % ii(j)=a.maxA*exp(-a.alpha^2*(t(j)-a.t1)^2);
        % didt(j)=(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
        
        % MTLL
                         ii(j)=a.maxA*(1 - (z(j)-a.H1)/(a.H2-a.H1))*exp(-a.alpha^2*(t(j)-a.t1)^2);
                         didt(j)=(1 - (z(j)-a.H1)/(a.H2-a.H1))*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
        
        % MTLE
        %ii(j)  = a.maxA*exp(-(z(j)-a.H1)/a.lamda)*exp(-a.alpha^2*(t(j)-a.t1)^2);
        %didt(j)= exp(-(z(j)-a.H1)/a.lamda)*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
         %x = (z-z(1))/range(z);
         %x.^(a-1).*((1-x.^a).^(b-1))
        
        % MTLKumaraswami
        %ii(j)=a.maxA*((z(j)-a.H1)/(a.H2-a.H1)).^(a.a-1).*((1-((z(j)-a.H1)/(a.H2-a.H1)).^a.a).^(a.b-1))*exp(-a.alpha^2*(t(j)-a.t1)^2);
        %didt(j)=((z(j)-a.H1)/(a.H2-a.H1)).^(a.a-1).*((1-((z(j)-a.H1)/(a.H2-a.H1)).^a.a).^(a.b-1))*(-2)*a.maxA*a.alpha^2*(t(j)-a.t1)*exp(-a.alpha^2*(t(j)-a.t1)^2);
        
        
        
        
        
    else
        %                 % MTLEI
          %             ii(j)=a.maxA*exp((z(j)-a.H1)/a.lamda)*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
         %                didt(j)=exp((z(j)-a.H1)/a.lamda)*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        %
        %                 % TL
        % ii(j)=a.maxA*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        % didt(j)=(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        %
        %                 % MTLL
                        ii(j)=a.maxA*(1 - (z(j)-a.H1)/(a.H2-a.H1))*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
                         didt(j)=(1 - (z(j)-a.H1)/(a.H2-a.H1))*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        
        % MTLE
        %ii(j)=a.maxA*exp(-(z(j)-a.H1)/a.lamda)*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        %didt(j)=exp(-(z(j)-a.H1)/a.lamda)*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        
         %MTLKumaraswami
         %ii(j)=a.maxA*((z(j)-a.H1)/(a.H2-a.H1)).^(a.a-1).*((1-((z(j)-a.H1)/(a.H2-a.H1)).^a.a).^(a.b-1))*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
         %didt(j)=((z(j)-a.H1)/(a.H2-a.H1)).^(a.a-1).*((1-((z(j)-a.H1)/(a.H2-a.H1)).^a.a).^(a.b-1))*(-2)*a.maxA*(a.alpha/a.k)^2*(t(j)-a.t1).*exp(-(a.alpha/a.k)^2*(t(j)-a.t1)^2);
        
    end
    
end