function IBP_auto_modeler

%% User inputs
% Wave features
tt1 = 6.3000e-006;        % Time up to first peak
tt2 = 21.200e-006;        % Time up to the second peak
tt3 = 4.7200e-005;        % Total pulse duration
pk1 = 1.1629;             % first peak magnitude
pk2 = -0.3371;          % Second peak magnitude

% Time steps
dt = 1e-6;


% Initial parameters
a=struct( ...
    'maxA' , {-420}        , ...
    't1' , {10.6e-006}    , ...
    't2' , {15.0e-006}     , ...
    'v' ,  {9.6e7}             , ...
    'H1' , {7600}           , ...
    'H2' , {8400}          , ...
    'x0' , {0}              , ...
    'y0' , {0}              , ...
    'T1' , {0.1500e-004}    , ...
    'T2' , {2.5e-004}     , ...
    't_step' , {1.0000e-007}, ...
    'lamda', {113}, ...
    'x' , {50000}            , ...
    'y' , {0}               ,...
    'dh', {100}             ...    
    );
a.alpha = 10/a.t2;

% Is the answer acceptable?
tot_accept = 0;
t1_accept = 0;

P = 5; % matching percentage


figure
hold all
n = 0;

while ~tot_accept  
    
    %% Changing t1
    while ~t1_accept
        
    n = n + 1;
   [t,E_stat,E_ind,E_rad] = IBP_modeler(a);
   
   % Scaling
   [E_stat,E_ind,E_rad,a.maxA] = scale_it(E_stat,E_ind,E_rad,a.maxA,pk1,pk2);

   % Get features
   feat = get_features(t,E_stat,E_ind,E_rad);
   
   % Changing a.t1 significantly change the total pulse time
   totT = feat.tot_end_t - feat.tot_st_t;
   if tt3 < (100+P)/100*totT && tt3 > (100-P)/100*totT
       t1_accept = 1;
       fprintf('%i Optimum t1 = %0.1f us\n',n,a.t1*1e6)
   elseif tt3 < totT
       a.t1 = a.t1 + dt; t1_accept = 0;       
       fprintf('%i Increasing t1 to %0.1f us\n',n,a.t1*1e6)
   else
       a.t1 = a.t1 - dt; t1_accept = 0;
       fprintf('%i Decreasing t1 to %0.1f us\n',n,a.t1*1e6)
   end
   
   cla
   plot(t,feat.E_tot) 
    end
   
   %% Changing t2
   n = n+1;
   [t,E_stat,E_ind,E_rad] = IBP_modeler(a);
   
   % Scaling
   [E_stat,E_ind,E_rad,a.maxA] = scale_it(E_stat,E_ind,E_rad,a.maxA,pk1,pk2);

   % Get features
   feat = get_features(t,E_stat,E_ind,E_rad);
   
   % Changing a.t2 significantly change the second peak magnitude
   peak2 = abs(feat.tot_min);
    totT = feat.tot_end_t - feat.tot_st_t;
   if abs(pk2) < (100+P)/100*peak2 && abs(pk2) > (100-P)/100*peak2
        t2_accept = 1;
        fprintf('%i Optimum t2 = %0.1f\n',n,a.t2*1e6)
   elseif pk2 < peak2 && tt3 > totT
       a.t2 = a.t2 + dt; t2_accept = 0;
       fprintf('%i Increasing t2 to %0.1f us\n',n,a.t2*1e6)
   elseif pk2 > peak2 && tt3 < totT
       a.t2 = a.t2 - dt; t2_accept = 0;
       fprintf('%i Decreasing t2 to %0.1f us\n',n,a.t2*1e6)
   else
       fprintf('%i Didn''t change t2\n',n)
       t2_accept =  0;
   end

   %% Increasing Speed makes positive peak more peaker
   n = n+1;
   [t,E_stat,E_ind,E_rad] = IBP_modeler(a);
   
   % Scaling
   [E_stat,E_ind,E_rad,a.maxA] = scale_it(E_stat,E_ind,E_rad,a.maxA,pk1,pk2);

   % Get features
   feat = get_features(t,E_stat,E_ind,E_rad);
   
   % Changing v significantly change the second peak magnitude
   peak1 = abs(feat.tot_max);
  
   if abs(pk1) < (100+P)/100*peak1 && abs(pk1) > (100-P)/100*peak1
        v_accept = 1;
        fprintf('%i Optimum v = %0.0f m/s\n',n,a.v)
   elseif pk1 > peak1
       a.v = a.v + 1e7; v_accept = 0;
       fprintf('%i Increasing v to %0.0f m/s\n',n,a.v)
   else
       a.v = a.v - 1e7; v_accept = 0;
       fprintf('%i Decreasing v to %0.0f m/s\n',n,a.v)
   end
   
   
   
%    figure
%    hold all
   %cla
   plot(t,feat.E_tot)  

   
   if t1_accept && t2_accept && v_accept
       tot_accept = 1;
   end
   
end

function feat = get_features(t,E_stat,E_ind,E_rad)
    
    % test percentage
    p = 1;

    %Radiation term statistics
    [feat.rad_max ind1] = max(E_rad);
    [feat.rad_min ind2] = min(E_rad);
    feat.rad_maxt = t(ind1);
    feat.rad_mint = t(ind2);
    
    
    st_ind = sum(E_rad(1:ind1) <= p/100*feat.rad_max);
    end_ind = sum(E_rad(ind2:end) <= p/100*feat.rad_min)+ind2;
    
    try; feat.rad_st_t = t(st_ind); catch; feat.rad_st_t = NaN; end
    try; feat.rad_end_t = t(end_ind); catch; feat.rad_end_t = NaN; end
    
    % Induction term statistics
    [feat.ind_min ind1] = min(E_ind);
    feat.ind_mint = t(ind1);
    
    st_ind = sum(E_ind(1:ind1) >= p/100*feat.ind_min); 
    end_ind = sum(E_ind(ind1:end) <= p/100*feat.ind_min)+ind1;
    
    try; feat.ind_st_t = t(st_ind); catch; feat.ind_st_t = NaN; end
    try; feat.ind_end_t = t(end_ind); catch; feat.ind_end_t = NaN; end
   
   
    % Static term statistics
    feat.stat_min  = min(E_stat);

    %Total statistics
    E_tot = E_rad + E_stat + E_ind;
    [feat.tot_max ind1] = max(E_tot);
    [feat.tot_min ind2] = min(E_tot);
    feat.tot_maxt = t(ind1);
    feat.tot_mint = t(ind2);
    
    
    st_ind = sum(E_tot(1:ind1) <= p/100*feat.tot_max);
    end_ind = sum(E_tot(ind2:end) <= p/100*feat.tot_min)+ind2;
    
    try; feat.tot_st_t = t(st_ind); catch; feat.tot_st_t = NaN; end
    try; feat.tot_end_t = t(end_ind); catch; feat.tot_end_t = NaN; end
    feat.E_tot = E_tot;
    
function [E_stat,E_ind,E_rad,maxA] = scale_it(E_stat,E_ind,E_rad,maxA,pk1,pk2)
   % Lets scale
   E_tot = E_stat + E_ind + E_rad;
   range1 = range(E_tot);
   range2 = abs(pk1 - pk2);
   maxA = maxA*range2/range1;
   %E_tot = E_tot*range2/range1;
   E_rad = E_rad*range2/range1;
   E_ind = E_ind*range2/range1;
   E_stat = E_stat*range2/range1;
