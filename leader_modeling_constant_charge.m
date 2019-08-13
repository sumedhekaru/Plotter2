function leader_modeling_constant_charge
% This model charge will be deposite at the optimum previous PBFA points
% 2013 -06-11 This is a veriation of leader modeling 4. Here only constant
% charges will be used intead of intensity biased method.
clc


%% User inputs
% Sensor number to model
sn = [1 2 3 6];

% Base folder
bf = 'C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F02\ki_sqrd_using_22bit\';

% HSVP file name
% hsvp_fn = 'F18-14Aug_231217-50k8mm320x240.txt';
hsvp_fn = 'F02-14Aug_213025-50k8mm320x240.txt';

% PBFA file name
pbfa_fn = 'pbfa.txt';

% Figures to turn on
histogram       = 1;
I_vs_t          = 1;
q_distribution  = 1;
E_change        = 1;
E_change_indiv  = 1;    % Plot indivisual E-change plots
dE_dt_plots     = 1;    % Plot dE for individual time frames vs. time
dE_dt_plot      = 1;    % single dE_dt plot
input_output    = 1;    % Figure with input output data
ki_sqrd_plot    = 1;    % Ki squired with time
save_output     = 1;    % Save model data in to a matlab structure

% Threshold to remove charge dynamic range to inprove the plot visibility
q_dist_thr = - 0.007;

% save figures?
save_figure = 1;
% Comments to include at the end of the file name
comment    = 'model_with_PBFA_10kHz_data';        

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
    
    %fprintf('%0.6f\t%0.6f\t%0.6f\n',min(rd.t(i,:)),max(rd.t(i,:)),b.t(ind(1)))
    vshift = interp1(rd.t(i,:),rd.v(i,:),b.t(ind(1)));
    rd.v(i,:) = rd.v(i,:) - vshift;
end


totQ = 0;
totabsQ = 0;

dQs = zeros(1,L0);
dEms = dQs;

Qs = zeros(1,L-1);
ki_squires = NaN(1,L-1);
b.qs = b.x-b.x;

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

% Save outputs at each step
stepOut(L).dEms = [];
stepOut(L).dEcs = [];
stepOut(L).ind  = [];
stepOut(L).kis  = [];
stepOut(L).dQ   = [];
stepOut(L).Ems  = [];
stepOut(L).Ecs  = [];
stepOut(L).totQ = [];
stepOut(L).pbfaInd = [];



stepOut(1).Ems  = zeros(1,L0);
stepOut(1).Ecs  = zeros(1,L0);

for i = 2:L-1
    
    abort = multiWaitbar( 'Overall Process', i/L);
    
    if abort
        multiWaitbar( 'CloseAll' );
        disp('User stopped the program ''LeaderModeling'' ')
        return        
    end
    
    ind = find(b.frameN == b.frameN(endInd(i)+1));
    
    for j = 1:L0          
        % measured dE    
        dEms(j) = interp1(rd.t(j,:),rd.v(j,:),tuni(i))-interp1(rd.t(j,:),rd.v(j,:),tuni(i-1));
        %dEms(j) = interp1(rd.t(j,:),rd.v(j,:),tuni(i))-E(j,i-1);
        dEsm2(j,i) = dEms(j);       
    end
    
    stepOut(i).ind  = ind; 
    stepOut(i).dEms = dEms;
    stepOut(i).Ems  = stepOut(i-1).Ems + dEms;
    
    % find dQ
    %dQ = find_dQ1(dEms,b,ind,endInd(i),sen_set,sn,-.1,.1); 
    
    for j = 1:L1
        multiWaitbar( 'Charges for each PBFA', j/L1);
        b.xp = b.pbfa(j,2); b.yp = b.pbfa(j,3); b.zp = b.pbfa(j,4);
        [dQ_vals(j) dE_vals(j,:) N ki_vals(j)] = ...
            find_dQ4(dEms,b,ind,endInd(i),sen_set,sn,-.1,.1,0);        
    end
    
    % get the optimum values
    [mm mind] = min(ki_vals);   
    
    ki_squires(i) = mm;
    dQ   = dQ_vals(mind);
    b.pbfa(mind,5) = b.pbfa(mind,5) + dQ; 
    dEcs = dE_vals(mind,:);
    
    stepOut(i).kis      = mm;
    stepOut(i).dQ       = dQ;
    stepOut(i).pbfaInd  = mind;
    stepOut(i).dEcs     = dEcs;   
    
    
    b.qs(ind) = b.qs(ind)+ dQ*b.I(ind)/sum(b.I(ind));
       
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
    
    stepOut(i).totQ = totQ;
    stepOut(i).Ecs  = E(:,i);
end


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

output.tuni = tuni;
output.Qs = Qs;



% Space charge distribution plot
qs = b.qs;
if q_distribution
    
    Nsr = sum(qs<q_dist_thr);
    qs(find(qs < q_dist_thr)) = q_dist_thr;
    qs(find(qs > 0)) = 0;
    qs = abs(qs);
    
    figure
    hold all
    % data = open('C:\Users\daqop\Desktop\Sumedhe\HSV analasis\20110814\F18\background.mat');
    % imagesc(data.img,[0, 65535])
    scatter(b.N_x,240-b.N_y,12,qs,'filled')

    cbh = colorbar;
    cbhy = get(cbh,'ylabel');
    set(cbhy,'String','dQ','fontsize',10)
    daspect([1 1 1])
    box on
    title(sprintf('Charge distribution along channels (shrinked %i points at %0.4fC )',Nsr,q_dist_thr))
    
    
    % Plotting PBFA data      
    b = ldar2scr(b);
    plot(b.pbfa(:,6),240-b.pbfa(:,7),'k.')
    
    indx = find(b.pbfa(:,5) == 0);
    tempy = b.pbfa(indx,6);
    tempy(indx) = NaN;
    
   
    scatter(tempy,240-b.pbfa(:,7),abs(b.pbfa(:,5)/min(b.pbfa(:,5)))*1000,[1 0 0])
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

output.rd = rd;
output.E = E;


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

output.dEsm2 = dEsm2;
output.dEsc2 = dEsc2;
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

output.totQ = totQ;
output.totabsQ = totabsQ;


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

if ki_sqrd_plot
    figure
    plot(tuni(1:L-7),ki_squires(1:L-7));
    xlabel('Time (s)')
    ylabel('\chi^2 (V^2/m^2)')
    title(sprintf('Max = %0.3f   Min = %0.3f  Avg = %0.8f ',...
        max(ki_squires(1:L-7)),min(ki_squires(1:L-7)),nanmean(ki_squires(1:L-7))))
    
    save_this(save_figure,bf,['Mod-' '-ki-sqrd'] ,comment)
end

% Save outputs
output.b = b;
output.ki_squires = ki_squires;
output.stepOut = stepOut;

if save_output
    save([bf 'ModelDataOut' '_' comment '.mat'],'-Struct','output')
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
        
        kis(i) = sum(((dEcs - dEms)./dEms).^2);
        %kis(i) = sum((dEcs - dEms).^2);
    end
    
    [mm indx] = min(kis);  
    
    if N > 100 || (q2 - q1) < 10000*eps
        dQ = (q1+q2)/2;      
        return
    end
    
    [dQ N] = find_dQ2(dEms,b,ind,endInd,sen_set,sn,qs(indx)-dq,qs(indx)+dq,N);
    

function save_this(save_figure,bf,fn,comment)

    if save_figure
        fn1 = [bf fn '_' comment '.fig'];
        
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
        
        kis(i) = sum(((dEcs - dEms)./dEms).^2);
        %kis(i) = sum((dEcs - dEms).^2);
    end
    
    [mm indx] = min(kis);  
    
    if N > 100 || (q2 - q1) < 10000*eps
        dQ = (q1+q2)/2;      
        return
    end
    
    [dQ dEcs N] = find_dQ3(dEms,b,ind,endInd,sen_set,sn,qs(indx)-dq,qs(indx)+dq,N);     

function b = ldar2scr(b)

    x = b.pbfa(:,2);
    y = b.pbfa(:,3);
    z = b.pbfa(:,4);

    f = 8.0e-3;  %Focal length
    p = 20.0e-6; %Fixel Size

    b.pbfa(:,6) = b.scr_x-(f*tan(atan((y-b.cam_y)./(x-b.cam_x ))- atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x )) ))/p;

    b.pbfa(:,7) = b.scr_y - (z - b.cam_z).*f./(p.*sqrt((x-b.cam_x ).^2+(y-b.cam_y ).^2 ))...
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
    
    % Measurement errors
    mErr = [0.32 0.31 0.31 NaN 0.18 0.29 0.53 0.16 0.13 0.17 NaN]/998;
    
    % Degrees of freedom
    %df = 2;
    
    for i = 1:4
        
        for j = 1:L0
            dEcs(j) = -qs(i) * sum(k*b.I(ind)/sum(b.I(ind)).*b.z(ind)./((b.x(ind)-sen_set.x(sn(j))).^2+(b.y(ind)-sen_set.y(sn(j))).^2+b.z(ind).^2).^1.5) + ...
                      qs(i) * k* b.zp./((b.xp-sen_set.x(sn(j))).^2+(b.yp-sen_set.y(sn(j))).^2+b.zp.^2).^1.5;                 
        end
        
        % kis(i) = sum(((dEcs - dEms)./dEms).^2);
        kis(i) = sum(((dEcs - dEms)./mErr(sn)).^2)/2;
        % kis(i) = sum((dEcs - dEms).^2);
        % kis(i) = sum(abs(dEcs - dEms)./abs(dEms)); % Lu et al 2011
    end
    
    [ki indx] = min(kis);  
    
    
    if N > 100 || (q2 - q1) < 10000*eps || sum(isnan(kis))> 0
        dQ = (q1+q2)/2;      
        return
    end
    
    [dQ dEcs N ki] = find_dQ4(dEms,b,ind,endInd,sen_set,sn,qs(indx)-dq,qs(indx)+dq,N);     
        
 
        
function b = load_pbfa_data(fn,b)
    
    fid = fopen(fn);
    data=textscan(fid,'%f %f %f %f %f %f %f %s','HeaderLines',2);
    fclose(fid);

    b.pbfa = cell2mat([data(:,2),data(:,3),data(:,4),data(:,5)]);


    
    
    
    
 