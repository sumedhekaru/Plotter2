function leader_modeling5

% This model charge will be deposite at the optimum previous PBFA points.
% Also it will try to match the early part of the leader too.
%clc


%% User inputs
% Sensor number to model
sn = [2 6 3];

% Base folder
bf = 'C:\Users\Sumedhe\Desktop\LeaderModeling2013\F18-4th-RS\';

% HSVP file name
hsvp_fn = 'F18-14Aug_231217-50k8mm320x240.txt';

% PBFA file name
pbfa_fn = 'pbfa.txt';

% Figures to turn on
histogram       = 0;
I_vs_t          = 0;
q_distribution  = 1;
E_change        = 1;
E_change_indiv  = 0;    % Plot indivisual E-change plots
dE_dt_plots     = 0;    % Plot dE for individual time frames vs. time
dE_dt_plot      = 1;    % single dE_dt plot
input_output    = 1;    % Figure with input output data

% Threshold to remove charge dynamic range to inprove the plot visibility
q_dist_thr = - 0.02;

% Last PBFA point
lastPbfa = [83537.2774920	-17431.0	-4545.0	5387.0];

% save figures?
save_figure = 0;
% Comments to include at the end of the file name
comment    = 'blind_upper_leader_modeling';        

%% Start the program
hsvp_fno = hsvp_fn;
% Load settings file
sen_set = open('sensor_setting.mat');

% Get HSVP data
hsvp_fn = [bf hsvp_fn];
b = hsvp_extract(hsvp_fn);

% Obtain pbfa data
b = load_pbfa_data([bf pbfa_fn],b);

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


b =  scr2ldar(b);

% remove propagation time
b.t = b.t - sqrt((b.cam_x - b.ref_x).^2+(b.cam_y - b.ref_y).^2) /3.0e8;

%% Working for first part of the leader
% Last PBFA data point
b.lastPbfa = lastPbfa;

% first HSVP data point
b.firstHSVP = [b.t(1) b.x(1) b.y(1) b.z(1)];

% Number of steps (100us steps)
N0 = round((b.firstHSVP(1) - b.lastPbfa(1))/200e-6);

% Assign charge locations points
tn = b.lastPbfa(1):(b.firstHSVP(1)-b.lastPbfa(1))/N0:b.firstHSVP(1);
xn = b.lastPbfa(2):(b.firstHSVP(2)-b.lastPbfa(2))/N0:b.firstHSVP(2);
yn = b.lastPbfa(3):(b.firstHSVP(3)-b.lastPbfa(3))/N0:b.firstHSVP(3);
zn = b.lastPbfa(4):(b.firstHSVP(4)-b.lastPbfa(4))/N0:b.firstHSVP(4);

b.t = [tn(2:end-1)' ; b.t];
b.x = [xn(2:end-1)' ; b.x];
b.y = [yn(2:end-1)' ; b.y];
b.z = [zn(2:end-1)' ; b.z];
b.I = [tn(2:end-1)' - tn(2:end-1)' + 1; b.I];
b.frameN = [(b.frameN(1)-N0+1:b.frameN(1)-1)'; b.frameN];

N0 = N0 - 2;

[N_x N_y] = ldar2scr(b,xn(2:end-1),yn(2:end-1),zn(2:end-1));
b.N_x = [N_x' ; b.N_x];
b.N_y = [N_y' ; b.N_y];


%% Cotinuing program
[tuni endInd] = unique(b.t);
L = length(tuni);
E = zeros(L0,L-1);

%ind = find(b.frameN == b.frameN(1));
ind = 1;
prevInd = 1;


% First value
%rd.v = smooth(rd.v);

for i = 1:L0
    % Remove propagation time
    tshift = sqrt((sen_set.x(sn(i)) - b.ref_x)^2 + (sen_set.y(sn(i)) - b.ref_y)^2 ) / 3.0e8;
    rd.t(i,:) = rd.t(i,:) - tshift;
    
    
    vshift = interp1(rd.t(i,:),rd.v(i,:),b.t(ind(1)));
    rd.v(i,:) = rd.v(i,:) - vshift;
end


totQ = 0;
totabsQ = 0;

dQs = zeros(1,L0);
dEms = dQs;

Qs = zeros(1,L-1);
b.qs = b.x;

% Change in E-chane at each step
dEsm2 = nan(L0,L-1);
dEsc2 = dEsm2;

% number of PBFA points
L1 = length(b.pbfa(:,1));
ki_vals = zeros(1,L1);
dQ_vals = ki_vals;
dE_vals = zeros(L1,L0);

% charges at each pbfa point
b.pbfa(:,5) = 0;

multiWaitbar( 'Overall Process', 0,'CanCancel','on' );
multiWaitbar( 'Charges for each PBFA',0 ,'Color', [0.2 0.9 0.3],'CanCancel','on');

% speed = nan(1,L-1);

for i = 2:L-1
    
    abort = multiWaitbar( 'Overall Process', i/L);
    
    if abort
        multiWaitbar( 'CloseAll' );
        disp('User stopped the program ''LeaderModeling'' ')
        return        
    end
    
    if i <= N0
        ind = i;
    else    
        ind = find(b.frameN == b.frameN(endInd(i)+1-N0));
    end
    
    for j = 1:L0          
        % measured dE    
        dEms(j) = interp1(rd.t(j,:),rd.v(j,:),tuni(i))-interp1(rd.t(j,:),rd.v(j,:),tuni(i-1));
        %dEms(j) = interp1(rd.t(j,:),rd.v(j,:),tuni(i))-E(j,i-1);
        dEsm2(j,i) = dEms(j);       
    end
    
    % find dQ
    %dQ = find_dQ1(dEms,b,ind,endInd(i),sen_set,sn,-.1,.1); 
    
    for j = 1:L1
        multiWaitbar( 'Charges for each PBFA', j/L1);
        b.xp = b.pbfa(j,2); b.yp = b.pbfa(j,3); b.zp = b.pbfa(j,4);
        [dQ_vals(j) dE_vals(j,:), N ki_vals(j)] = ...
            find_dQ4(dEms,b,ind,endInd(i),sen_set,sn,-.1,.1,0);        
    end
    
    % get the optimum values
    [mm mind] = min(ki_vals);   
    
    dQ   = dQ_vals(mind);
    b.pbfa(mind,5) = b.pbfa(mind,5) + dQ; 
    dEcs = dE_vals(mind,:);
    
    
    b.qs(ind) = dQ*b.I(ind)/sum(b.I(ind));
    
    % charge 
    %totQ = totQ + dQ*length(b.I(ind));
    totQ = totQ + dQ;    
    totabsQ = totabsQ + abs(dQ);
    Qs(i) = dQ;
       
    
    for j = 1:L0        
        E(j,i) = E(j,i-1)+dEcs(j);
        %E(j,i) = interp1(rd.t(j,:),rd.v(j,:),tuni(i-1))+dEc;
        dEsc2(j,i) = dEcs(j); 
    end   
    
%     % trying to find speed
%     % difference matrix
%     sL1 = length(ind);
%     sL2 = length(prevInd);
%     rdif = nan(sL1,sL2);
%     
%     v = [];
%     
%     for i0=1:sL1
%         for j0 = 1:sL2
%             rdif(i0,j0) = sqrt((b.x(ind(i0))-b.x(prevInd(j0)))^2 + ...
%                                (b.y(ind(i0))-b.y(prevInd(j0)))^2 + ...   
%                                (b.z(ind(i0))-b.z(prevInd(j0)))^2);
%         end
%         
%         smm  = min(rdif(i0,:));
%         
%         if smm > 50 && smm < 100
%             v = [v smm/20e-6];
%         end
%     end
%     
%     speed(i) = mean(v);
%     prevInd = ind;
end

% figure
% plot(speed,'r.')


multiWaitbar( 'CloseAll' );

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
    hold all
    % data = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F18\background.mat');
    % imagesc(data.img,[0, 65535])
    
    size(b.qs)
    size(b.N_x)
    
    scatter(b.N_x,240-b.N_y,12,b.qs,'filled')

    cbh = colorbar;
    cbhy = get(cbh,'ylabel');
    set(cbhy,'String','dQ','fontsize',10)
    daspect([1 1 1])
    box on
    title(sprintf('Charge distribution along channels (shrinked %i points at %0.4fC )',Nsr,q_dist_thr))
    
    
    % Plotting PBFA data      
    [N_x N_y] = ldar2scr(b,b.pbfa(:,2),b.pbfa(:,3),b.pbfa(:,4));    
    b.pbfa(:,6) = N_x;
    b.pbfa(:,7) = N_y;
    
    plot(b.pbfa(:,6),240-b.pbfa(:,7),'k.')
    
    indx = find(b.pbfa(:,5) == 0);
    b.pbfa(indx,6) = NaN;
    
   
    scatter(b.pbfa(:,6),240-b.pbfa(:,7),abs(b.pbfa(:,5)/min(b.pbfa(:,5))*400),[1 0 0])
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

if input_output
    figure
    hold all
    box on
    set(gca,'XTickLabel',[],'YTickLabel',[])
    title('Input Output')
    
    inpstr = sprintf('bf = ''%s'';\nhsvp_fn = ''%s'';\n%pbfa_fn = ''%s''\nq_dist_thr = %0.4f;\ncomment = ''%s'';\n',...
        bf,hsvp_fno,pbfa_fn,q_dist_thr,comment);
    
    annotation('textbox', [0.15 0.5 0.72 0.4],'String', inpstr,'Interpreter','none');
    annotation('textbox', [0.15 0.15 0.72 0.3],'String', str,'Interpreter','none');  
    
    save_this(save_figure,bf,['Mod-Input_output'] ,comment)  
    
end



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

    bf = [bf comment '\'];
    
    if ~exist(bf,'dir')
        mkdir(bf)
    end
    
    if save_figure
        fn1 = [bf fn '.fig'];
        
        if exist(fn1, 'file')
            button = questdlg([fn1 ' Exist! Do you want to replace?'],'File exisit','Yes','No','No');
            
            if ~strcmp(button,'Yes')
                return
            end
        end
            
        saveas(gcf,fn1)
        close(gcf)
    end
        
function  [dQ dEcs N] = find_dQ3(dEms,b,ind,endInd,sen_set,sn,q1,q2,N)
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
            dEcs(j) = -qs(i) * sum(k*b.I(ind)/sum(b.I(ind)).*b.z(ind)./((b.x(ind)-sen_set.x(sn(j))).^2+(b.y(ind)-sen_set.y(sn(j))).^2+b.z(ind).^2).^1.5) + ...
                      qs(i) * k* b.zp./((b.xp-sen_set.x(sn(j))).^2+(b.yp-sen_set.y(sn(j))).^2+b.zp.^2).^1.5;                 
        end
        
        %kis(i) = sum(((dEcs - dEms)./dEms).^2);
        kis(i) = sum((dEcs - dEms).^2);
    end
    
    [mm indx] = min(kis);  
    
    if N > 100 || (q2 - q1) < 10000*eps
        dQ = (q1+q2)/2;      
        return
    end
    
    [dQ dEcs N] = find_dQ3(dEms,b,ind,endInd,sen_set,sn,qs(indx)-dq,qs(indx)+dq,N);     

function [N_x N_y] = ldar2scr(b,x,y,z)

    f = 8.0e-3;  %Focal length
    p = 20.0e-6; %Fixel Size

    N_x = b.scr_x-(f*tan(atan((y-b.cam_y)./(x-b.cam_x ))- atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x )) ))/p;

    N_y = b.scr_y - (z - b.cam_z).*f./(p.*sqrt((x-b.cam_x ).^2+(y-b.cam_y ).^2 ))...
        ./sqrt(1+(tan(atan((y-b.cam_y)./(x-b.cam_x ))- atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x )))).^2);
    
     
    
function  [dQ dEcs N ki] = find_dQ4(dEms,b,ind,endInd,sen_set,sn,q1,q2,N)
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
            dEcs(j) = -qs(i) * sum(k*b.I(ind)/sum(b.I(ind)).*b.z(ind)./((b.x(ind)-sen_set.x(sn(j))).^2+(b.y(ind)-sen_set.y(sn(j))).^2+b.z(ind).^2).^1.5) + ...
                      qs(i) * k* b.zp./((b.xp-sen_set.x(sn(j))).^2+(b.yp-sen_set.y(sn(j))).^2+b.zp.^2).^1.5;                 
        end
        
        %kis(i) = sum(((dEcs - dEms)./dEms).^2);
        kis(i) = sum((dEcs - dEms).^2);
    end
    
    [ki indx] = min(kis);  
    
    if N > 100 || (q2 - q1) < 10000*eps
        dQ = (q1+q2)/2;      
        return
    end
    
    [dQ dEcs N ki] = find_dQ4(dEms,b,ind,endInd,sen_set,sn,qs(indx)-dq,qs(indx)+dq,N);     
        
 
        
function b = load_pbfa_data(fn,b)
    
    fid = fopen(fn);
    data=textscan(fid,'%f %f %f %f %f %f %f %s','HeaderLines',2);
    fclose(fid);

    b.pbfa = cell2mat([data(:,2),data(:,3),data(:,4),data(:,5)]);


    
    
    
    
 