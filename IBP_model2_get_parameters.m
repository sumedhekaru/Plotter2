function IBP_model2_get_parameters

%% User inputs

% Initial parameters

maxA=-420;
t1=12.0e-006;
t2=30.0e-006;
v=8.e7;
H1=7600;
H2=8400;
x0=0;
y0=0;
t_step=1.0000e-006;
lamda=140   ;
dh=100;  

sns = [1 2 3 6 7 8 9 10];

% sensor locations
sns_x =  1.0e+04 *[ -0.7524, 0.3372, -0.3149, -0.4266, -1.1658, -0.2701, -2.7640, ...
                 -6.0091, 0.1825, -5.7394, -2.0637];
             
sns_y =  1.0e+04 *[1.6555, 0.4446, -0.6838, 0.9545, -1.7020, 0.2631, 4.9254, ...
                    -3.3983, -5.3008, 1.1923, 0.1569];

N = length(sns);

for i = 1:N 
    
   alpha = 10/t2;
   
   % Pulse arrival time to the sensor location
   t0 = sqrt((x0-sns_x(sns(i)))^2+(y0-sns_y(sns(i)))^2+H1^2)/3e8;
   T1 = t0 - 0.5e-4;
   T2 = t0 + 1.5e-4;
   
   
   [t,E_stat,E_ind,E_rad] = IBP_modeler2(maxA,t1,t2,v,H1,H2,x0,y0,T1,...
       T2,t_step,lamda,sns_x(sns(i)),sns_y(sns(i)),dh,alpha);
   
   figure(100)
   cla
   plot(t,E_stat+E_ind+E_rad)

   % Get features
   feat = get_features(t,E_stat,E_ind,E_rad);
   
   offsets(i) = feat.offset;
   maxs(i)    = feat.tot_max;
   mins(i)    = feat.tot_min;
   maxts(i)   = feat.tot_maxt;
   mints(i)   = feat.tot_mint;
   start_ts(i) = feat.tot_st_t;
   end_ts(i)   = feat.tot_end_t; 
    
end
clc
offsets'
maxs'
mins'
maxts'
mints'
start_ts'
end_ts'

function feat = get_features(t,E_stat,E_ind,E_rad)
    
    % test percentage
    p = 1;

    m1 = max(E_stat); m2 = min(E_stat);
    
    if     abs(m1) > abs(m2);        feat.offset = m1;
    else   feat.offset = m2;
    end
    

    %Total statistics
    E_tot = E_rad + E_stat + E_ind;
    [feat.tot_max, ind1] = max(E_tot);
    [feat.tot_min, ind2] = min(E_tot);
    feat.tot_maxt = t(ind1);
    feat.tot_mint = t(ind2);
    
    % pulse duration 
    st_ind = sum(E_tot(1:ind1) <= p/100*feat.tot_max);
    end_ind = sum(E_tot(ind2:end) <= p/100*feat.tot_min)+ind2;
    
    try feat.tot_st_t = t(st_ind); catch; feat.tot_st_t = NaN; end
    try feat.tot_end_t = t(end_ind); catch; feat.tot_end_t = NaN; end