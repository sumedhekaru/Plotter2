function guess_IBP_para

%% User inputs

% Real data folder
rdf = 'C:\Users\sumedhe\Desktop\CG-20120814_77636_333/';

% Initial parameters
maxA=1000;
H1=5267;
H2=6500;
x0=-9534;
y0= 6080;
t_step=1.0000e-006;
dh=10;  

% Ranges of values
t1 = 7.0e-6;
t2 = 26.5e-6;
v = 1.0e7;
lamda = 20;
alpInd =25;

sns = [1, 2, 3, 6, 8, 9, 10];

%% Starting program
% sensor locations
sns_x =  1.0e+04 *[ -0.7524, 0.3372, -0.3149, -0.4266, -1.1658, -0.2701, -2.7640, ...
                 -6.0091, 0.1825, -5.7394, -2.0637];
             
sns_y =    1.0e+04 *[1.6555, 0.4446, -0.6838, 0.9545, -1.7020, 0.2631, 4.9254, ...
                    -3.3983, -5.3008, 1.1923, 0.1569];

sns_IDs = {'K02','K14','K24','WSB','BCC','K17','EDW','STC','FLT','OVD'};


sns_x = sns_x(sns);
sns_y = sns_y(sns);
sns_IDs = sns_IDs(sns);

% Load real data
L0 = length(sns);

for i = 1:L0
    rd(i) = open([rdf sns_IDs{i} '.mat']);
end

% Generate data
cdd = generate_data(maxA,t1,t2,v,H1,H2,alpInd,x0,y0,...
    t_step,lamda,sns_x,sns_y,dh,rd,L0);


% Get fit values
[diff, norm_diff, scale] = compare_data(rd, L0,cdd);

tit_str = sprintf('Manual Para search  \nmaxA   = % 6.1f kA\nH1     = % 6i\nH2     = % 6i\nx0     = % 6i\ny0     = % 6i\nt-step = % 6.1f\ndh     = % 6i\nt1     = % 6.1f\nt2     = % 6.1f\nv      = % 6.1e\nlamda = % 6i\nalpInd = % 6.1f\n\nFit values...\ndiff = % 6.1f\nNdiff = % 6.1f\n',...
                maxA*scale/1000,H1,H2,x0,y0,t_step*1e6,dh,t1*1e6,t2*1e6,v,lamda,alpInd,diff,norm_diff);

% Plot figures
plot_figures(rd,cdd,sns_IDs,tit_str)



       
       
%% Functions       
function data = generate_data(maxA,t1,t2,v,H1,H2,alpInd,x0,y0,...
           t_step,lamda,sns_x,sns_y,dh,rd,L0)
   
   scale = zeros(1,L0);
   
   parfor j=1:L0
       % Changing parameters
       alpha = alpInd/t2;
       
       % Pulse arrival time to the sensor location
       t0 = sqrt((x0-sns_x(j))^2+(y0-sns_y(j))^2+H1^2)/3e8;
       rdt = rd(j).t;
       rdy = rd(j).y;
             
       [t,E_stat,E_ind,E_rad] = IBP_modeler2(maxA,t1,t2,v,H1,H2,x0,y0,...
           t0 - 50e-6, t0 + 100e-6 ,t_step,lamda,sns_x(j),sns_y(j),dh,alpha);
       
       E_tot = E_stat+E_ind+E_rad;
       
       [hpk1, dt1] = max(E_tot);
       [hpk2, dt2] = max(rdy);
       
       [lpk1, dt3] = min(E_tot);
       [lpk2, dt4] = min(rdy);
   
       %scale(j) = mean([hpk2/hpk1,lpk2/lpk1]);
       scale(j) = (hpk2-lpk2)/(hpk1-lpk1);
       
       
       % time shift
       if hpk2 > -lpk2
            t = t + rdt(dt2) - t(dt1);
       else
            t = t + rdt(dt4) - t(dt3);
       end
        
       
       % save data
       data(j).t = t;
       data(j).y = E_tot;      

   end
   
   data(1).scale = mean(scale);
   
function plot_figures(rd,cdd,sns_IDs,tit_str)

L = length(sns_IDs);



if     L == 1; m = 1; n = 2; 
elseif L == 2; m = 1; n = 3;
elseif L <= 5; m = 2; n = 3;
elseif L <= 7; m = 2; n = 4;
elseif L <= 8; m = 3; n = 3;
else m = 3; n = 4;
end

sz = get(0, 'ScreenSize');
fg = figure(1431);
clf
%set(fg,'visible','off'); 
set(fg, 'Position', [0 0 sz(3) sz(4) ] );

for i = 1:L
    subplot(m,n,i); hold all; box on;
    plot(rd(i).t,rd(i).y);
    plot(cdd(i).t,cdd(i).y*cdd(1).scale);
    xlim([min(rd(i).t) max(rd(i).t)])
    title(sns_IDs{i})
    legend('Real','Calc')
end

subplot(m,n,L+1); box on
text(0.1,0.5,tit_str) 
set(gca,'XTick',[])
set(gca,'YTick',[])

%title([rd(1).rd '    ' sec2hhmmss(rd(1).t(1))])
set(fg,'visible','on')

function [diff, norm_diff, scale] = compare_data(rd, L0,data)
  
   diff = 0;
   norm_diff = 0;
   scale = data(1).scale;
   for j = 1:L0
              
       % Get calculated data for real data
       caly = interp1(data(j).t,data(j).y,rd(j).t)*scale;
       
%        figure; hold all;       
%        plot(rd(j).t,caly)
%        plot(rd(j).t,rd(j).y)
       
       % differnce
       diff = diff + sum(abs(caly - rd(j).y));
       
       factor = 1/range(rd(j).y);
       
       norm_diff = norm_diff + sum(abs(caly - rd(j).y))*factor;

   end