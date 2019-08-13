function TL_model_01
% Trying to follow Uman et al. 1975 figure 3

clc

%% User Inputs
a=struct( ...
    'maxA' , {10000}       , ...    % Peak Current
    't1' , {2.5000e-006}   , ...    % Rise time
    't2' , {22.500e-006}   , ...    % Fall time (not total time)
    'v' , {8e7}            , ...    % Speed of the current
    'H1' , {0}             , ...    % Starting height
    'H2' , {4000}          , ...    % Ending height
    'x0' , {0}             , ...    % Horizontal x coordinate of the current
    'y0' , {0}             , ...    % Horizontal y coordinate of the current
    'T1' , {.000e-003}      , ...    % Plotting range begining (aproximately D/c)
    'T2' , {.10000e-003}   , ...    % Plotting range end  (aproximately D/c+t1+t2 or little more)
    't_step' , {1.0000e-007} , ...  % Time step
    'x' , {1000}           , ...    % Sensor x coordinate
    'y' , {0}              , ...    % sensor y coordinate
    'dh', {1.0}            , ...    % height step (use small as possible for smooth E_rad curve
    'tshift_to_zero',{0}   , ...    % make it '1' if you need to shift plot to start from zero (still need apropriate values for T1 and T2)
    'plotsOn', {[1 1 1 1]} , ...    % Turn on plots you need {[Static Induc Radia Total]}
    'plotInCurrentFig', {0}  ...    % make it '1' if you need to plot it on current figure instead of a new figure
    );

%% Start program
eps0    = 8.8542e-12;               % Permitivity of free space
c       = 299792458/1.0003;         % Speed of light in air

% get the distance from the pulse to sensors
a.x = a.x - a.x0;
a.y = a.y - a.y0;
D=sqrt(a.x^2+a.y^2);
a.D = D;

% Counter
i=0;

% Pre create variables to store data
ln=length(a.T1:a.t_step:a.T2);
Z = (a.H2:-a.dh:a.H1)';
E_ind=zeros(1,ln);
E_rad=E_ind;
E_stat=E_ind;

maxEstat = [];

started = 0;

wbh = waitbar(0,'Please wait...','Name','NBP Simulator');

for t=a.T1:a.t_step:a.T2
    
    i=i+1;
    try
        waitbar(i/ln,wbh,sprintf('Please Wait... %.0f%%',i*100/ln))
    catch
        disp('TL model Program Stopped by User!')
        return
    end
    
    % Current, derivation of current, point vise integration of the current  
    [ii,didt,intI]=current(Z,t-(Z-a.H1)/a.v-((D^2+Z.^2).^0.5)/c,a);
    
    % Radiation term
    Y = (-1/2/pi/eps0)*D^2./(c^2*(D^2+Z.^2).^1.5).*didt;    
    E_rad(i) = trapz(Z,Y);
  
    % Induction Term
    Y = (1/2/pi/eps0)*((2*Z.^2-D^2)./(c*((Z.^2+D^2).^2))).*ii;
    E_ind(i) = trapz(Z,Y); 
   
    
    % Electrostatic Term     
     Y = (1/2/pi/eps0)*((2*Z.^2-D^2)./((Z.^2+D^2).^2.5)).*intI;
     E_stat(i) = trapz(Z,Y);
     
     [mm pind] = max(ii);      
     
     if pind ~= 1 && ~started
         started = 1;
     end
     
     if started && pind == 1 && isempty(maxEstat)
         maxEstat = E_stat(i);
     end
     
        
     figure(101)
     subplot(3,1,1)
     plot(E_stat)
     subplot(3,1,2)
     plot(Z,intI)
     ylim([0 0.15])
     subplot(3,1,3)
     plot(Z,ii)
     
     % When the current goes to the top of the channel, intI(2) will be
     % none zero. We will use this property to identify when it goes to top
     % and then keep static field constant after that assuming all the
     % charge accumulate at the top of the channel.
     
%      if intI(2) ~= 0 
%          if isempty(maxEstat)
%              maxEstat = E_stat(i-1);
%              %maxEstat = intI; 
%          end
%     end
     
     if ~isempty(maxEstat)
         E_stat(i) = maxEstat;
         %Y = (1/2/pi/eps0)*((2*Z.^2-D^2)./((Z.^2+D^2).^2.5)).*maxEstat;
         %E_stat(i) = trapz(Z,Y);
     end
     
end

 %Y = (1/2/pi/eps0)*(Q*4000/((4000^2+D^2).^1.5))

delete(wbh)

% Need to plot in the current figure?
if a.plotInCurrentFig
    figure(gcf)
else
    figure
    xlim([a.T1,a.T2])
end

hold all
box on

t=(a.T1:a.t_step:a.T2);

% Need to time shift to start from the origin?
if a.tshift_to_zero
    t = t - D/c;
end

% Plot apropriate plots
lg ={};
if a.plotsOn(1); plot(t,E_stat); lg = [lg 'Stat']; end 
if a.plotsOn(2); plot(t,E_ind) ; lg = [lg 'Induc']; end
if a.plotsOn(3); plot(t,E_rad) ; lg = [lg 'Radi']; end
if a.plotsOn(4); plot(t,(E_rad+E_ind+E_stat)*1) ; lg = [lg 'Total']; end
legend(lg)
%tools2fig


% function val = HH(H)
%   
%     
%     eps0    = 8.8542e-12;               % Permitivity of free space
%     c       = 299792458/1.0003;         % Speed of light in air
%     %size(H)
%     [ii,didt,intI]=current(H,a.t-(H-a.H1)/a.v-((a.D^2+H.^2).^0.5)/c,a);
%     
%     val = (-1/2/pi/eps0)*a.D^2./(c^2*(a.D^2+H.^2).^1.5).*didt;
% end


    

function [ii,didt,intI]=current(z,t,a)

    ii=zeros(size(z));
    didt= ii;
    intI = zeros(size(t));
    
    % Trangular current
    ind1 = sum(t <= 0);
    ind2 = sum(t <= a.t1);
    ind3 = sum(t <= a.t1 + a.t2);
    
    ii(1:ind1)= 0;
    ii(ind1+1:ind2) = a.maxA*t(ind1+1:ind2)/a.t1;
    ii(ind2+1:ind3) = a.maxA - a.maxA*(t(ind2+1:ind3)-a.t1)/a.t2;
    ii(ind3+1:end) = 0;
    

    
    didt(1:ind1)= 0;
    didt(ind1+1:ind2) = a.maxA/a.t1;
    didt(ind2+1:ind3) = -a.maxA/a.t2;
    didt(ind3+1:end) = 0;
    
    if length(t) > 1
        intI = cumtrapz(t,ii);
    end
    
end

% function [ii,didt,intI]=current2(z,t,a)
%     
%     L = length(z);
%     ii=zeros([L 1]);
%     didt= ii;
%     intI = ii;
%     
%     
%     for kk=1:L
%         if t(kk) > 0 && t(kk) <=a.t1
%             ii(kk) = a.maxA*t(kk)/a.t1;
%             didt(kk) = a.maxA/a.t1;
%         elseif t(kk) <= a.t2
%             ii(kk) = a.maxA - a.maxA*(t(i)-a.t1)/a.t2;
%             didt(kk) = -a.maxA/a.t2;
%         end
%     end
%     
%     if length(t) > 1
%         intI = cumtrapz(t,ii);
%     end
%     
% end

end

        
        


