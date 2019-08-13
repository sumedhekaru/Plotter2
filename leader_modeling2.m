function leader_modeling2
% This will only work with a single sensor

clc

%% User inputs
% Sensor number to model
sn = [1 2 6 3];

% Base folder
bf = 'C:/Users/daqop/Desktop/Sumedhe/HSV analasis/20110814/F01/2ndRS/';

% HSVP file name
hsvp_fn = 'F01-14Aug_212917-50k8mm320x240.txt';

% Figures to turn on
histogram       = 0;
I_vs_t          = 0;
q_distribution  = 1;
E_change        = 0;
E_change_indiv  = 0;    % Plot indivisual E-change plots
dE_dt_plots     = 0;    % Plot dE for individual time frames vs. time
dE_dt_plot      = 0;    % single dE_dt plot

% Threshold to remove charge dynamic range to inprove the plot visibility
q_dist_thr = - 0.09;

% save figures?
save_figure = 0;
% Comments to include at the end of the file name
comment    = '';        

%% Start the program

% Load settings file
sen_set = open('sensor_setting.mat');

% Get HSVP data
hsvp_fn = [bf hsvp_fn];
b = hsvp_extract(hsvp_fn);

% E-chage data file
rd.t = [];
rd.v = [];
L0 = length(sn);

for i = 1:L0
    E_fn = sprintf('%s%s_10KHz_leader_data.mat',bf,sen_set.sen_IDs{sn(i)});
    
    % real data
    d = open(E_fn);
    %size(d.t)
    rd.t = [rd.t; d.t];
    rd.v = [rd.v; d.v]; 
    
end

%rd.v(1,:) = smooth(rd.v(1,:));

b =  scr2ldar(b);
k = 1/2/pi/8.85418782e-12;

% remove propagation time
b.t = b.t - sqrt((b.cam_x - b.ref_x).^2+(b.cam_y - b.ref_y).^2) /3.0e8;


[tuni endInd] = unique(b.t);
L = length(tuni);
E = zeros(L0,L-1);

ind = find(b.frameN == b.frameN(1));

% First value
%rd.v = smooth(rd.v);

for i = 1:L0
    % Remove propagation time
    tshift = sqrt((sen_set.x(sn(i)) - b.ref_x)^2 + (sen_set.y(sn(i)) - b.ref_y)^2 ) / 3.0e8;
    rd.t(i,:) = rd.t(i,:) - tshift;
    
    
    vshift = interp1(rd.t(i,:),rd.v(i,:),b.t(ind(1)));
    rd.v(i,:) = rd.v(i,:) - vshift;
end

wbh = waitbar(0,'Working on...','name','Modeling leader');

totQ = 0;
totabsQ = 0;

dQs = zeros(1,L0);
dEms = dQs;

Qs = zeros(1,L-1);
b.qs = b.x;

% Change in E-chane at each step
dEsm2 = nan(L0,L-1);
dEsc2 = dEsm2;


for i = 2:L-1
    
    try
        waitbar(i/L,wbh,sprintf('Working on...(%0.2f%%)',i/L*100));
    catch
        return
    end
    
    ind = find(b.frameN == b.frameN(endInd(i)+1));
    
    for j = 1:L0          
        % measured dE    
        dEms(j) = interp1(rd.t(j,:),rd.v(j,:),tuni(i))-interp1(rd.t(j,:),rd.v(j,:),tuni(i-1));
        %dEms(j) = interp1(rd.t(j,:),rd.v(j,:),tuni(i))-E(j,i-1);
        dEsm2(j,i) = dEms(j);       
    end
    
    % find dQ
    %dQ = find_dQ1(dEms,b,ind,endInd(i),sen_set,sn,-.1,.1);    
    dQ = find_dQ2(dEms,b,ind,endInd(i),sen_set,sn,-.1,.1,0); 
    
    b.qs(ind) = dQ*b.I(ind)/sum(b.I(ind));
    
    % charge 
    %totQ = totQ + dQ*length(b.I(ind));
    totQ = totQ + dQ;    
    totabsQ = totabsQ + abs(dQ);
    Qs(i) = dQ;
    
    for j = 1:L0
        
        x0 = sen_set.x(sn(j));
        y0 = sen_set.y(sn(j));
        
        % calculated dE
        dEc = sum(k*dQ*b.I(ind)/sum(b.I(ind)).*b.z(ind)./((b.x(ind)-x0).^2+(b.y(ind)-y0).^2+b.z(ind).^2).^1.5);    
        %dEc = sum(k*dQ*(zeros(size(b.I(ind)))+1).*b.z(ind)./((b.x(ind)-x0).^2+(b.y(ind)-y0).^2+b.z(ind).^2).^1.5);    
        E(j,i) = E(j,i-1)+dEc;
        %E(j,i) = interp1(rd.t(j,:),rd.v(j,:),tuni(i-1))+dEc;
        dEsc2(j,i) = dEc; 
    end       
end

delete(wbh)

if histogram
    figure
    hist(Qs)    
    save_this(save_figure,bf,'Mod-q_histogram',comment)
end

% Current vs. time plot
if I_vs_t
    figure
    plot(tuni(1:L-1),abs(Qs/20e-6/1000))
    ylabel('Current (kA) Or Charge (\times0.02 C)')
    xlabel('time')
    title('Current vs. time')
    save_this(save_figure,bf,'Mod-I_vs_t',comment)
end


% Space charge distribution plot
if q_distribution

    Nsr = sum(b.qs<q_dist_thr);
    b.qs(find(b.qs < q_dist_thr)) = q_dist_thr;
    b.qs = abs(b.qs);
    figure
    % hold all
    % data = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F18\background.mat');
    % imagesc(data.img,[0, 65535])
    scatter(b.N_x,240-b.N_y,12,b.qs,'filled')
    cbh = colorbar;
    cbhy = get(cbh,'ylabel');
    set(cbhy,'String','dQ','fontsize',10)
    daspect([1 1 1])
    box on
    title(sprintf('Charge distribution along channels (shrinked %i points at %0.4fC )',Nsr,q_dist_thr))
    save_this(save_figure,bf,'Mod-q_distibution',comment)
end


if E_change
    % Field plots
    figure
    hold all
    lg = {};
    shift = 0;
    
    for i=1:L0
        
        plot(rd.t(i,:),rd.v(i,:)+shift,'LineWidth',2)
        plot(tuni(1:L-1),E(i,:)+shift,'--','LineWidth',2)
        shift = shift +  floor(min(E(i,:)/10))*10;
        lg = [ lg sen_set.sen_IDs{sn(i)} sprintf('%s Modeled',sen_set.sen_IDs{sn(i)})];
        
    end
    
    legend(lg,'Location','SouthWest')
    tools2fig
    xlim([min(tuni) max(tuni)])
    title(sprintf('E-change due to leader (%0.2f ms -- %0.2f C)',range(tuni)*1000,totQ))
    box on
    xlabel('t (s)')
    ylabel('V/m')
    
    save_this(save_figure,bf,'Mod-All_E-changes',comment)
end

% Individual field plots
if E_change_indiv
    for i=1:L0
        figure
        hold all
        plot(rd.t(i,:),rd.v(i,:),'LineWidth',2)
        plot(tuni(1:L-1),E(i,:),'--','LineWidth',2)
        legend({sen_set.sen_IDs{sn(i)} ,'Modeled'});
        legend('Location','SouthWest')
        tools2fig
        xlim([min(tuni) max(tuni)])
        title(sprintf('E-change due to leader (%0.2f ms -- %0.2f C)',range(tuni)*1000,totQ))
        box on
        xlabel('t (s)')
        ylabel('V/m') 
        save_this(save_figure,bf,['Mod-' sen_set.sen_IDs{sn(i)} 'E-change'] ,comment)
    end    
end

% dE dt plots
if dE_dt_plots
    for i = 1:L0
        figure
        hold all
        plot(tuni(1:L-10),dEsm2(i,1:L-10))
        xlim([min(tuni) max(tuni(1:L-10))])
        plot(tuni(1:L-1),dEsc2(i,:))
        legend(sprintf('%s Real',sen_set.sen_IDs{sn(i)}),'Model')
        xlabel('Time (s)')
        ylabel('dE (V/m)')
        title('E-Change at each step')
        box on      
        save_this(save_figure,bf,['Mod-' sen_set.sen_IDs{sn(i)} 'dE-dt'] ,comment)
    end
end

% dE_dt_plots together
if dE_dt_plot
    
    % Field plots
    figure
    hold all
    lg = {};
    shift = 0;
    
    for i=1:L0
        
        plot(tuni(1:L-10),dEsm2(i,1:L-10)+shift)
        plot(tuni(1:L-1),dEsc2(i,:)+shift)
        shift = shift +  floor(min(dEsm2(i,1:L-10)/0.5))*0.5;
        lg = [ lg sen_set.sen_IDs{sn(i)} sprintf('%s Modeled',sen_set.sen_IDs{sn(i)})];
        
    end
    
    legend(lg,'Location','SouthWest')
    tools2fig
    xlim([min(tuni(1:L-10)) max(tuni(1:L-10))])
    title(sprintf('E-change at each step(%0.2f ms -- %0.2f C)',range(tuni)*1000,totQ))
    box on
    xlabel('t (s)')
    ylabel('V/m')
    
    save_this(save_figure,bf,['Mod-All_dE-dt'] ,comment)    
    
end

    

str = sprintf('dT                    \t= %0.1f ms\n',range(tuni)*1000);
str = sprintf('%sNumber of frames       \t= %i\n',str,length(tuni));
str = sprintf('%sNumber of HSVPs        \t= %i\n',str,length(b.t));
str = sprintf('%sMax HSVP altitude      \t= %0.2f km\n',str,max(b.z)/1000);
str = sprintf('%sTotal charge tran      \t= %0.2f C\n',str,totQ);
str = sprintf('%sTotal abs charge tran  \t= %0.2f C\n',str,totabsQ);
str = sprintf('%sAvg speed              \t= %0.0f m/s (%0.4fc)\n',str, max(b.z)/range(tuni),max(b.z)/range(tuni)/3.0e8);
fprintf(str)




% figure
% plot3(x,y,z,'r*')
% daspect([1 1 1])
% box on
% xlabel('x (m)')
% ylabel('y (m)')
% zlabel('z (m)')
% 
% figure
% plot(t,z,'rp')



function b = scr2ldar(b)
N_x = b.N_x;
N_y = b.N_y;

f = 8.0e-3;  %Focal length
p = 20.0e-6; %Fixel Size


r = sqrt((((b.ref_x-b.cam_x)).^2+(b.ref_y-b.cam_y ).^2 ).*(1+(p^2.*(b.scr_x-N_x ).^2)./f^2 ));

theta = atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x ))+atan(p*(b.scr_x-N_x)/f);

b.x = b.cam_x + r.*cos(theta);
b.y = b.cam_y + r.*sin(theta);
b.z = b.cam_z + r.*((b.scr_y-N_y).*p)./f;

% Convert frame numbers to time
b.t = b.align_time + (b.frameN - b.align_frame)./b.frame_rate;


    
 function dQ = find_dQ(dEm,b,ind,endInd,x0,y0)
     k = 1/2/pi/8.85418782e-12;
     
dQ = dEm / sum(k*b.I(ind)/sum(b.I(ind)).*b.z(ind)./((b.x(ind)-x0).^2+(b.y(ind)-y0).^2+b.z(ind).^2).^1.5);


function  dQ = find_dQ1(dEms,b,ind,endInd,sen_set,sn,q1,q2)
    
  
    qs = q1:(q2-q1)/5000:q2;
    
    k = 1/2/pi/8.85418782e-12;
    L0 = length(sn);
    dEcs = zeros(1,L0); 
    
    kis = zeros(size(qs));
    
    
    
    for i = 1:length(qs)
        
        for j = 1:L0
            dEcs(j) = qs(i) * sum(k*b.I(ind)/sum(b.I(ind)).*b.z(ind)./((b.x(ind)-sen_set.x(sn(j))).^2+(b.y(ind)-sen_set.y(sn(j))).^2+b.z(ind).^2).^1.5);
            %dEcs(j) = qs(i) * sum(k*(zeros(size(b.I(ind)))+1).*b.z(ind)./((b.x(ind)-sen_set.x(sn(j))).^2+(b.y(ind)-sen_set.y(sn(j))).^2+b.z(ind).^2).^1.5);
        end
        
        kis(i) = sum((dEcs - dEms).^2);
    end
    
    [mm indx] = min(kis);
    
    %figure
    %plot(qs,kis)
    
    dQ = qs(indx);
 
  
 function  [dQ N] = find_dQ2(dEms,b,ind,endInd,sen_set,sn,q1,q2,N)
%%    
    N = N + 1;
    

   
    dq = (q2 - q1)/4;
    qs = q1:dq:q2;
    
    k = 1/2/pi/8.85418782e-12;
    L0 = length(sn);
    dEcs = zeros(1,L0); 
    
    kis = zeros(1,4);
    
    for i = 1:4
        
        for j = 1:L0
            dEcs(j) = qs(i) * sum(k*b.I(ind)/sum(b.I(ind)).*b.z(ind)./((b.x(ind)-sen_set.x(sn(j))).^2+(b.y(ind)-sen_set.y(sn(j))).^2+b.z(ind).^2).^1.5);
            %dEcs(j) = qs(i) * sum(k*(zeros(size(b.I(ind)))+1).*b.z(ind)./((b.x(ind)-sen_set.x(sn(j))).^2+(b.y(ind)-sen_set.y(sn(j))).^2+b.z(ind).^2).^1.5);
        end
        
        %kis(i) = sum(((dEcs - dEms)./dEms).^2);
        kis(i) = sum((dEcs - dEms).^2);
    end
    
    [mm indx] = min(kis);  
    
    if N > 100 || (q2 - q1) < 10000*eps
        dQ = (q1+q2)/2;      
        return
    end
    
    [dQ N] = find_dQ2(dEms,b,ind,endInd,sen_set,sn,qs(indx)-dq,qs(indx)+dq,N);
    

function save_this(save_figure,bf,fn,comment)

    if save_figure
        saveas(gcf,[bf fn '_' comment '.fig'])
        close(gcf)
    end
        
        
 
        

 