function test_NBP_model
% This function is written to test NBP model I submitted for 2016 papar.
% Idea is to pick up parameters that produce the correct charge moment
% values.

%% Test results
% d = open('C:\Users\sumedhe\Desktop\NBP-timing-JGR-2015\Modeling\para_search\20160409-1247NBP_para_search.mat');
% 
% figure
% plot(d.P)
% plot(d.norm_diff)
% ind = 6526;
% 
% H2 = d.H2(ind)
% t1 = d.t1(ind)
% t2 = d.t2(ind)
% AlpInd = d.AlpInd(ind)
% lamda = d.lamda(ind)
% return

%% User inputs
d.version=1;
             d.msg= 'Ready...';
             d.dir = 'C:\Users\sumedhe\Desktop\NBP-timing-JGR-2015\Modeling\20150819\';
         d.sns_str= '2,3,6';
              d.dt= 1;
              d.dh= 50;
              d.T1= 50;
              d.T2= 100;
           d.modal= 'MTLE';
       d.modal_num= 2;
         d.curType= 'Positive';
    d.currType_num= 2;
              d.x0= 4.9533e+03;
              d.y0= -7.6397e+03;
              d.H1= 13090;
              d.Hm= 5267;
              %d.H2= 12090;
              %d.t1= 2;
              %d.t2= 20;
               %d.v= 6;
          %d.lamda1= 235;
          d.lamda2= 2.7000;
          %d.AlpInd= 10;      
            d.maxA= -1000;
               d.P= 0;

% Changing parameters
H2s = 10590:200:11590;
t1s = 2:2:16;
t2s = 10:5:50;
AlpInds = 5:3:25;
lamdas = 50:50:300;
vs = 1:3:13;

% parameter save file
results_file = ['C:\Users\sumedhe\Desktop\para_test\' datestr(now,'yyyymmdd-HHMM') 'NBP_para_search.mat'];

res.H2 = nan;
res.t1 = nan;
res.t2 = nan;
res.AlpInd = nan;
res.lamda = nan;
res.diff = nan;
res.norm_diff = nan;
res.P = nan;

save(results_file,'-Struct','res')

number_of_permutations = length(H2s)*length(t1s)*length(t2s)*length(AlpInds)*length(lamdas)*length(vs)
estimated_time = sec2hhmmss(0.56*number_of_permutations)


%% Starting program
% save previous differences

sns = str2double({'2','3','6'});

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

tmin = inf;
tleng = inf;

for i = 1:L0
    rd(i) = open([d.dir sns_IDs{i} '.mat']);
    [~, tmpind] = min(rd(i).y);
    tmp1 = rd(i).t(tmpind);
    tmp2 = range(rd(i).t);
    
    if tmp1 < tmin;      tmin = tmp1;    end    
    if tmp2 < tleng;     tleng = tmp2;  end
end

% figure
% hold all
% for i = 1:L0
%     plot(rd(i).t,rd(i).y)
% end


% Round tmin to lowest milisecond
%tleng = tleng*1.2;
tmin = round(tmin*1e6)/1e6;
%tleng = round(tleng*2e5)/2e5;


try temp = d.tt2;
catch
    d.tt2 = 20;
    d.dtt2 = .1;
    d.t3 = 1000;
    d.dt3 = 1;
end

H2 = 10590:200:11590;
t1s = 2:2:16;
t2s = 10:5:50;
AlpInd = 5:2:25;
lamdas = 50:50:500;

counter = 0;
% seve after each ## iterations
save_iter = 5;
ind = 0;

diffx = NaN(1,save_iter);
norm_diffx = NaN(1,save_iter);
Px = NaN(1,save_iter);
H2x = NaN(1,save_iter);
t1x = NaN(1,save_iter);
t2x = NaN(1,save_iter);
AlpIndx = NaN(1,save_iter);
lamdax = NaN(1,save_iter);

tic
for ind1 = 1:length(H2);    d.H2 = H2s(ind1);
for ind2 = 1:length(t1s);   d.t1 = t1s(ind2);
for ind3 = 1:length(t2s);     d.t2 = t2s(ind3);
for ind4 = 1:length(AlpInds); d.AlpInd = AlpInds(ind4);
for ind5 = 1:length(lamdas);  d.lamda1 = lamdas(ind5);
for ind6 = 1:length(vs);  d.v = vs(ind6);
    counter = counter + 1;
    ind = ind + 1;
    
    [diff, norm_diff,scale,P,cdd] = gen_data(d,sns_x,sns_y,rd,L0);
    
    diffx(ind) = diff;
    norm_diffx(ind) = norm_diff;
    Px(ind) = P;
    H2x(ind) = H2s(ind1);
    t1x(ind) = t1s(ind2);
    t2x(ind) = t2s(ind3);
    AlpIndx(ind) = AlpInds(ind4);
    lamdax(ind) = lamdas(ind5);
    
    if mod(counter,save_iter) == 0
        t_elapsed = toc;
        % we need to save the results
        %Open the file
        res = open(results_file);
        res.diff = [res.diff diffx];
        res.norm_diff = [res.norm_diff norm_diffx];
        res.P = [res.P Px];
        res.H2 = [res.H2 H2x];
        res.t1 = [res.t1 t1x];
        res.t2 = [res.t2 t2x];
        res.AlpInd = [res.AlpInd AlpIndx];
        res.lamda = [res.lamda lamdax];
        save(results_file,'-Struct','res');
        ind = 0;
        
        % Print some results
        fprintf('\nElapsed time = %0.1f s\n',t_elapsed);
        fprintf('Elalpsed steps = %i (%0.3f %%)\n',counter,counter/number_of_permutations*100);
        rem_time = (t_elapsed/counter)*(number_of_permutations-counter);
        hh = floor(rem_time/3600);
        mm = floor((rem_time-hh*3600)/60);
        fprintf('Remaining time = %2.2i:%2.2i\n',hh,mm);
        
    end

% Ends for 6 for loops
end
end                   
end
end
end
end

                    
%[diff, norm_diff,scale,P,cdd] = gen_data(d,sns_x,sns_y,rd,L0);

%figure(1979)

%set(gcf,'units','normalized','outerposition',[0.05 0.05 0.9 0.9])

clf
hold all
for i = 1:3
    subplot(1,3,i)
    hold all
    plot(cdd(i).t,cdd(i).y*scale)
    plot(rd(i).t,rd(i).y)
    box on
    xlim(tmin+[-25 125]*1e-6)
end

tools2fig

function [diff, norm_diff,scale,P,cdd] = gen_data(d,sns_x,sns_y,rd,L0)

% Generate data
cdd = generate_data(d.maxA,d.t1*1e-6,d.t2*1e-6,d.tt2*1e-6,d.t3*1e-6,d.v*1e7,d.H1,d.Hm,d.H2,d.AlpInd,d.x0,d.y0,...
    d.dt*1e-6,d.lamda1,d.lamda2,sns_x,sns_y,d.dh,rd,L0,d.T1*1e-6,d.T2*1e-6,d.modal);

% Get fit values
[diff, norm_diff, scale] = compare_data(rd, L0,cdd);

P = cdd(1).P*scale;


%toffset = sec2hhmmss(tmin);






function [diff, norm_diff, scale] = compare_data(rd, L0,data)
  
   diffs = nan(1,L0);
   norm_diffs = nan(1,L0);
   scale = data(1).scale;
   
   for j = 1:L0
        
       % Get calculated data for real data
       caly = interp1(data(j).t,data(j).y,rd(j).t)*scale;
       
       % differnce
       diffs(j) =  sqrt(nanmean((caly - rd(j).y).^2));
       
       factor = 1/range(rd(j).y);
       
       norm_diffs(j) = sqrt(nanmean((caly - rd(j).y).^2))*factor;

   end
   
   diff = mean(diffs);
   norm_diff = mean(norm_diffs);
   
function data = generate_data(maxA,t1,t2,tt2,t3,v,H1,Hm,H2,alpInd,x0,y0,...
           t_step,lamda,lamda2,sns_x,sns_y,dh,rd,L0,T1,T2,method)
       
   
   scale = zeros(1,L0);
   
   % if it is MTLK we need to estimate H1 first
   if strcmp(method,'MTLK')       
       x = 0:0.0001:1;
       I = x.^(lamda-1).*((1-x.^lamda).^(lamda2-1));
       I = I/max(I);
       [~, ind] = max(I);       
       H1 = H2 - (Hm-H2)/(x(ind) - 1);
       data.x = x;
       data.I = I;
       data.z = H1 + (H2-H1)*x;
       
       % % Another way to find H1
       % H1 = (Hm - H2*((lamda - 1)/(lamda*lamda2 - 1))^(1/lamda))/...
       %        (1  -  ((lamda - 1)/(lamda*lamda2 - 1))^(1/lamda))
   end
  
   parfor j=1:L0
       t = [];
       E_stat = [];
       E_ind = [];
       E_rad = [];
       P = [];
       
       % Changing parameters
       alpha = alpInd/t2;
       
       % Pulse arrival time to the sensor location
       t0 = sqrt((x0-sns_x(j))^2+(y0-sns_y(j))^2+H1^2)/3e8;
       rdt = rd(j).t;
       rdy = rd(j).y;
       
       switch method
           case 'MTLE'
                  [t,E_stat,E_ind,E_rad,P] = IBP_mod_MTLE(maxA,t1,t2,v,H1,H2,x0,y0,...
                   t0 - T1, t0 + T2 ,t_step,lamda,sns_x(j),sns_y(j),dh,alpha);
           case 'MTLL'
               [t,E_stat,E_ind,E_rad,P] = IBP_mod_MTLL(maxA,t1,t2,v,H1,H2,x0,y0,...
                   t0 - T1, t0 + T2 ,t_step,lamda,sns_x(j),sns_y(j),dh,alpha);

            case 'MTLEI'
               [t,E_stat,E_ind,E_rad,P] = IBP_mod_MTLEI(maxA,t1,t2,v,H1,H2,x0,y0,...
                   t0 - T1, t0 + T2 ,t_step,lamda,sns_x(j),sns_y(j),dh,alpha);
           case 'MTLEID'
               [t,E_stat,E_ind,E_rad] = IBP_mod_MTLEID(maxA,t1,t2,v,H1,Hm,H2,x0,y0,...
                   t0 - T1, t0 + T2 ,t_step,lamda,lamda2,sns_x(j),sns_y(j),dh,alpha); 
           case 'MTLK'
               [t,E_stat,E_ind,E_rad,P] = IBP_mod_MTLK(maxA,t1,t2,v,H1,Hm,H2,x0,y0,...
                   t0 - T1, t0 + T2 ,t_step,lamda,lamda2,sns_x(j),sns_y(j),dh,alpha);
           case 'MTLEL'
               [t,E_stat,E_ind,E_rad,P,Is,didts] = IBP_mod_MTLEL(maxA,t1,t2,tt2,t3,v,H1,H2,x0,y0,...
                   t0 - T1, t0 + T2 ,t_step,lamda,sns_x(j),sns_y(j),dh,alpha);
                data(j).Is = Is;
                data(j).didts = didts;
           otherwise
               disp('This model is not working yet')                  
       end
       
       E_tot = E_stat+E_ind+E_rad;
       
       [hpk1, dt1] = max(E_tot);
       [hpk2, dt2] = max(rdy);
       
       [lpk1, dt3] = min(E_tot);
       [lpk2, dt4] = min(rdy);
   
       %scale(j) = mean([hpk2/hpk1,lpk2/lpk1]);
       scale(j) = (hpk2-lpk2)/(hpk1-lpk1);
       
       
       % time shift
       if hpk2 > -lpk2
            data(j).t = t + rdt(dt2) - t(dt1);
       else
            data(j).t = t + rdt(dt4) - t(dt3);
       end        
       
       % save data
       data(j).y = E_tot;  
       data(j).P = P;

   end
   
   data(1).scale = mean(scale);