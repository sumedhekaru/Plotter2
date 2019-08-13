function data = pulse_simulator(a)
% With original + 4 bouncing
%clc
% try
%      a=open('pulse_modeling_data.mat');
% catch
%     Create the parameter structure with default values
%     a=struct( ...
%         'A'     ,{-5000}    , ...
%         't1'    ,{8e-6}     , ...
%         't2'    ,{26e-6}    , ...
%         'v'     ,{2e8}      , ...
%         'H1'    ,{6600}     , ...
%         'H2'    ,{6800}     , ...
%         'x'     ,{5000}     , ...
%         'y'     ,{5000}     , ...
%         'lamda' ,{800}      , ...
%         'step'  ,{1e-6}  , ...
%         'n_b'   ,{1}        , ...
%         't_b'   ,{0.5e-6}   , ...
%         'rho'   ,{-0.5}     , ...
%         'tshift',{0}        , ...
%         'T1'    ,{1e-5}        , ...
%         'T2'    ,{10e-5}        , ...
%         'add_plots_on'  ,{1}     , ...
%         'E_tot_on'      ,{1}     , ...
%         'E_rad_on'      ,{1}     , ...
%         'E_ind_on'      ,{1}     , ...
%         'E_stat_on'     ,{1}     ...
%         );
%
%     save('pulse_modeling_data.mat','-struct', 'a');
% end

tic;
%close all



% a.maxA       = -50000;                    % Maximum Current
% a.t2      = 26e-6;                    % Width of the current
% a.t1      = 8e-6;                     % Current Increasing time
% a.v       = 2.0e8;                    % Speeed
% a.H1      = 6600;                    % Initial Height
% a.H2      = 6800;                    % Final Height
% a.lamda   = 800;                      % lamda serves a purpose similar to the characteristic length in RB
k       = (a.t2-a.t1)/a.t1;
alpha   = 10/a.t2;

eps0    = 8.8542e-12;               % Permitivity of free space
c       = 299792458/1.0003;         % Speed of light in air
% a.step    = .1e-6;                     % duration of a time step

% rho_t   = -0.5;                     % Reflection coefficient at the top
% rho_b   = -0.5;                     % Reflection coefficient at the bottom
% a.coeff_reflect     = -.5;
% a.t_reflect     = .5e-6;                   % Reflection time
% a.N     = 1;

% Distance and the plotting time range
%D=0;     a.T1=0;           a.T2=9e-5;
%D=2800;     a.T1=0;           a.T2=9e-5;
%D=20000;    a.T1=5e-5;     a.T2=15e-5;
%D=50000;    a.T1=1e-4;     a.T2=3e-4;
%D=200000;       a.T1=6.65e-4;     a.T2=7.10e-4;

% a.T1=5e-5;           a.T2=15e-5;

% If there is no reflections lets make reflection time to zero
if a.N == 1
    a.t_reflect = 0;
    a.coeff_reflect = 1;
end

% Sensor numbers to simulate
s_num = [5 ];
%s_num = a.s_num;

% Number of points in the space that simulation should do
ng=length(s_num);

sen_set = open('sensor_setting.mat');
x_temp = a.x;
y_temp = a.y;

% x,y positions of the sensors
a.x = zeros(1,ng);
a.y = zeros(1,ng);

for i = 1:ng;
    a.x(i) = sen_set.x(s_num(i));
    a.y(i) = sen_set.y(s_num(i));
end

% get the distance from the pulse to sensors
a.x = a.x - x_temp;
a.y = a.y - y_temp;




% Find number of raws and columns
if      ng == 1;    raw = 1;    col = 1;
elseif  ng == 2;    raw = 1;    col = 2;
elseif  ng <= 4;    raw = 2;    col = 2;
elseif  ng <= 8;    raw = 2;    col = 4;
elseif  ng == 9;    raw = 3;    col = 3;
elseif  ng == 10;   raw = 2;    col =5;
else
    disp('this can not handle more than 10 plots. Exiting...')
end








for iiii=1:ng
    
    D=sqrt(a.x(iiii)^2+a.y(iiii)^2);
    
    i=0;
    
    
    global intgr
    intgr=0;
    
    
    ln=length(a.T1:a.t_step:a.T2);
    
    Z = (a.H1:(a.H2-a.H1)/100:a.H2)';
    E_ind=zeros(1,ln);
    E_rad=zeros(1,ln);
    E_stat=zeros(1,ln);
    
    wbh = waitbar(0,'Please wait...','Name','NBP Simulator');
    
    for t=a.T1:a.t_step:a.T2
        
        i=i+1;
        try
            waitbar(i/ln,wbh,sprintf('Please Wait... %.0f%%',i*100/ln))
        catch
            disp('NBP Program Stopped by User!')
            return
        end
        
        ii=0;
        didt=0;
        
        % Radiation Term
        for kn=1:a.N
            [ii_old,didt_old]=current(Z,t-(kn-1)*a.t_reflect-(Z-a.H1)/a.v-(D^2+Z.^2).^0.5/c);
            ii=ii+a.coeff_reflect^(kn-1)*ii_old;
            didt=didt+a.coeff_reflect^(kn-1)*didt_old;
        end
        
        % 1/R radiation term
        %Y = (-1/2/pi/eps0)*D^2./(c^2*(D^2+Z.^2).^1.5).*didt;
        
        % 1/R^1.3 radiation term
        Y = (-1/2/pi/eps0)*D^2./(c^2*(D^2+Z.^2).^1.565).*didt;
        
        E_rad(i) = trapz(Z,Y);
        
        % Induction Term
        Y = (1/2/pi/eps0)*((2*Z.^2-D^2)./(c*((Z.^2+D^2).^2))).*ii;
        E_ind(i) = trapz(Z,Y);
        
        
        % Electro Static Term
        Y=(1/2/pi/eps0)*(2*Z.^2 - D^2)./((D^2+Z.^2).^2.5).*int_i(Z,t-(Z-a.H1)/a.v-(D^2+Z.^2).^0.5/c);
        E_stat(i) = trapz(Z,Y);
        
        % Dipole moment
        if i==ln
            Y=int_i(Z,t-(Z-a.H1)/a.v-(D^2+Z.^2).^0.5/c);
            dipole=trapz(Z,Y);
        end
    end
    
    
    %Get data for current vs t
    
    Z=[a.H1,(a.H1+a.H2)/2,a.H2];
    % Z(2) = 7183;
    
    I=zeros(ln,3);
    II=zeros(ln,a.N);
    i=0;
    if a.addi_plots(3)==1 || a.addi_plots(2)==1
        for t=a.T1:a.t_step:a.T2
            
            i=i+1;
            try
                waitbar(i/ln,wbh,sprintf('Please Wait! Collecting data for i vs t graph... %.0f%%',i*100/ln))
            catch
                disp('NBP Program Stopped by User!')
                return
            end
            
            for kki=1:a.N
                %             II(i,kki)=a.coeff_reflect^(kki-1)*current(Z(1),t-a.t_reflect*(kki-1)-(Z(1)-a.H1)/a.v-(D^2+Z(1).^2).^0.5/c);
                %             I(i,:)=I(i,:)+a.coeff_reflect^(kki-1)*current(Z,t-a.t_reflect*(kki-1)-(Z-a.H1)/a.v-(D^2+Z.^2).^0.5/c);
                II(i,kki)=a.coeff_reflect^(kki-1)*current(Z(1),t-a.t_reflect*(kki-1)-(Z(1)-a.H1)/a.v);
                I(i,:)=I(i,:)+a.coeff_reflect^(kki-1)*current(Z,t-a.t_reflect*(kki-1)-(Z-a.H1)/a.v);
                
            end
            
            
        end
    end
    
    
    delete(wbh)
    
    
    t=(a.T1:a.t_step:a.T2)+a.t_shift;
    
    %chose the current figure
    
    if sum([a.plots_on a.addi_plots])        
        if iiii==1
            fg1=gcf;
            subplot(raw,col,iiii)
        else
            figure(fg1)
            subplot(raw,col,iiii)
        end
    end
    
    E_tot = E_rad + E_stat + E_ind;
    
    data.E_tot = E_tot;
    data.t     = t;
    data.E_rad = E_rad;
    data.E_stat = E_stat;
    data.E_ind  = E_ind;
    
    %E_tot = E_rad + E_ind;
    
    if a.plots_on(4)==1
        hold all
        plot(t,E_tot)
        [mmax1 ind] =  max(E_tot);
        tmax1 = t(ind);
        
        [mmin1 ind] =  min(E_tot);
        tmin1 = t(ind);
        
        hold off
    end
    
    
    if a.plots_on(1)==1
        hold all
        plot(t,E_rad)
        [mmax2 ind] =  max(E_rad);
        tmax2 = t(ind);
        
        [mmin2 ind] =  min(E_rad);
        tmin2 = t(ind);
        
        
        hold off
    end
    
    if a.plots_on(3)==1
        hold all
        plot(t,E_ind)
        hold off
    end
    
    if a.plots_on(2)==1
        hold all
        plot(t,E_stat)
        hold off
    end
    
%     box on
%     xlabel('Time (s)')
%     ylabel('dE (V/m)')
    %legend('Rad + Ind','Rad','Ind')
%     title(['When sensor at ' num2str(a.x) ' m'])
%     hold all
%     plot(tmax1,mmax1,'or','markerfacecolor','r')
%     plot(tmin1,mmin1,'or','markerfacecolor','r')
%     plot(tmax2,mmax2,'og','markerfacecolor','g')
%     plot(tmin2,mmin2,'og','markerfacecolor','g')
    
%     str = sprintf('+ peak shift = %.2f us \n- peak shift = %0.2f us', ...
%         (tmax1-tmax2)*1e6,(tmin1-tmin2)*1e6);
%     annotation('textbox',[.6,.3,.3,.3],...
%         'String',str,'FitBoxToText','on')
    


    
    if a.addi_plots(3)==1
        
        if iiii==1
            fg2=figure;
            title('Current at the channel')
            xlabel('Time (s)')
            ylabel('Current (A)')
        else
            figure(fg2)
        end
        
        hold all
        plot(t',I(:,1),t',I(:,2),t',I(:,3))
        legend( regexp(num2str([a.H1,(a.H1+a.H2)/2,a.H2]),'  ','split'))
        hold off
        
    end
    
    if a.addi_plots(2) == 1
        
        leg={};
        
        if iiii==1
            fg3=figure;
            title('Bouncing current')
            xlabel('Time (s)')
            ylabel('Currentt (A)')
        else
            figure(fg3)
        end
        
        for i=1:a.N
            plot(t',II(:,i))
            hold all
            leg{i}=sprintf('%i',i);
        end
        
        leg=[leg 'sum'];
        plot(t,sum(II'))
        if iiii==1
            legend(leg)
        end
        hold off
        
    end
    
    if a.addi_plots(1) ==1
        
        
        %t=0.0000419:0.000000001:0.000043;
        %t=a.T1:a.t_step:a.T2;
        
        tit=sprintf('I=%iA \nD=%im \nT1=%.1fus \nt2=%.1us \nH1=%im \nH2=%im\nlamda = %.1f \nv=%.1gm/s \ndip-mom = %.1fC.m\nn_b=%i'...
            ,a.maxA,D,a.t1*1e6,a.t2*1e6,a.H1,a.H2,a.lamda,a.v,dipole,a.N);
        
        
        if iiii==1
            fg4=figure;
        else
            figure(fg4)
        end
        
        
        subplot(2,3,1)
        hold all
        plot(t,E_rad)
        xlabel('time (s)')
        ylabel('E-Radiation (V/m)')
        title('Radiation')
        hold off
        
        %figure
        subplot(2,3,2)
        hold all
        plot(t,E_ind)
        xlabel('time (s)')
        ylabel('E-Induction (V/m)')
        title('Induction')
        hold off
        
        %figure
        subplot(2,3,3)
        hold all
        plot(t,E_stat)
        xlabel('time (s)')
        ylabel('E-Static (V/m)')
        title('Static')
        hold off
        
        %figure
        subplot(2,3,4)
        hold all
        plot(t,E_rad+E_stat+E_ind)
        xlabel('time (s)')
        ylabel('E-Total (V/m)')
        title('Total')
        hold off
        
        % figure
        subplot(2,3,5)
        plot(t,E_rad)
        hold all
        plot(t,E_ind)
        plot(t,E_stat)
        plot(t,E_rad+E_stat+E_ind)
        title('Componets')
        %legend('Radiation','Induction','Static','Total')
        
        subplot(2,3,6)
        text(0,.5,tit)
        axis off
        box on
        hold off
        
    end
    
    
    if a.addi_plots(4) ==1
        
        tit=sprintf('I=%iA \nD=%im \nt1=%.1fus \nt2=%.1us \nH1=%im \nH2=%im\nlamda = %.1f \nv=%.1gm/s \ndip-mom = %.1fC.m\nn_b=%i'...
            ,a.maxA,D,a.t1*1e6,a.t2*1e6,a.H1,a.H2,a.lamda,a.v,dipole,a.N);
        
        figure
        
        text(0.25,0.5,tit)
        title('Pulse Simulation Parameters')
        axis off
        box on
        
    end
    
end

tit=sprintf('I=%iA \nD=%im \nt1=%.1fus \nt2=%.1us \nH1=%im \nH2=%im\nlamda = %.1f \nv=%.1gm/s \ndip-mom = %.1fC.m\nn_b=%i'...
    ,a.maxA,D,a.t1*1e6,a.t2*1e6,a.H1,a.H2,a.lamda,a.v,dipole,a.N);

%text(0,0.5,tit)


toc

    function [ii didt]=current(z,t)
        
        ii=zeros(size(z));
        didt=zeros(size(z));
        
        
        for j=1:length(t)
            %             if t(j) < 0;
            %                 ii(j)=0;
            %                 didt(j)=0;
            if t(j) <= a.t1
                % MTLEI
                %                 ii(j)=a.maxA*exp((z(j)-a.H1)/a.lamda)*exp(-alpha^2*(t(j)-a.t1)^2);
                %                 didt(j)=exp((z(j)-a.H1)/a.lamda)*(-2)*a.maxA*alpha^2*(t(j)-a.t1)*exp(-alpha^2*(t(j)-a.t1)^2);
                
                %                 % TL
                %                 ii(j)=a.maxA*exp(-alpha^2*(t(j)-a.t1)^2);
                %                 didt(j)=(-2)*a.maxA*alpha^2*(t(j)-a.t1)*exp(-alpha^2*(t(j)-a.t1)^2);
                
                % MTLL
                %                 ii(j)=a.maxA*(1 - (z(j)-a.H1)/a.H1)*exp(-alpha^2*(t(j)-a.t1)^2);
                %                 didt(j)=(1 - (z(j)-a.H1)/a.H1)*(-2)*a.maxA*alpha^2*(t(j)-a.t1)*exp(-alpha^2*(t(j)-a.t1)^2);
                
                % MTLE
                ii(j)=a.maxA*exp(-(z(j)-a.H1)/a.lamda)*exp(-alpha^2*(t(j)-a.t1)^2);
                didt(j)=exp(-(z(j)-a.H1)/a.lamda)*(-2)*a.maxA*alpha^2*(t(j)-a.t1)*exp(-alpha^2*(t(j)-a.t1)^2);
                
                
            else
                %                 % MTLEI
                %                 ii(j)=a.maxA*exp((z(j)-a.H1)/a.lamda)*exp(-(alpha/k)^2*(t(j)-a.t1)^2);
                %                 didt(j)=exp((z(j)-a.H1)/a.lamda)*(-2)*a.maxA*(alpha/k)^2*(t(j)-a.t1).*exp(-(alpha/k)^2*(t(j)-a.t1)^2);
                %
                %                 % TL
                %                 ii(j)=a.maxA*exp(-(alpha/k)^2*(t(j)-a.t1)^2);
                %                 didt(j)=(-2)*a.maxA*(alpha/k)^2*(t(j)-a.t1).*exp(-(alpha/k)^2*(t(j)-a.t1)^2);
                %
                %                 % MTLL
                %                 ii(j)=a.maxA*(1 - (z(j)-a.H1)/a.H1)*exp(-(alpha/k)^2*(t(j)-a.t1)^2);
                %                 didt(j)=(1 - (z(j)-a.H1)/a.H1)*(-2)*a.maxA*(alpha/k)^2*(t(j)-a.t1).*exp(-(alpha/k)^2*(t(j)-a.t1)^2);
                
                % MTLE
                ii(j)=a.maxA*exp(-(z(j)-a.H1)/a.lamda)*exp(-(alpha/k)^2*(t(j)-a.t1)^2);
                didt(j)=exp(-(z(j)-a.H1)/a.lamda)*(-2)*a.maxA*(alpha/k)^2*(t(j)-a.t1).*exp(-(alpha/k)^2*(t(j)-a.t1)^2);
                
            end
        end
    end


    function inti=int_i(z,t)
        
        tau=t-a.t_step:a.t_step/10:t;
        
        inti=zeros(size(z));
        
        for j=1:length(z)
            % t=t-(z-a.H1)/a.v;
            cr=zeros(size(tau));
            for kk=1:length(tau)
                for knn=1:a.N
                    cr(kk)=cr(kk)+a.coeff_reflect^(knn-1)*current(z(j),tau(knn)-a.t_reflect*(k-1));
                end
            end
            
            %length(tau)
            %length(cr)
            
            inti(j)=trapz(tau,cr);
            
        end
        
        inti=inti+intgr;
        
        intgr=inti;
        
    end
end



