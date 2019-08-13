function Ip = calculate_Ip(Vp,x,y,z,sn,sen_set,type)

% E-change meters were cross calibrated with CGLSS to find peak currents.
if isempty(Vp)
   Ip = NaN;
   return
end

%% Using Hilbert Transform using channel 3
if type == 1

    factor = [450.5 372.2 400.0 NaN NaN 393.9 646.9 202.0 162.2 206.3];
    
    R = sqrt((x - sen_set.x(sn)).^2+(y - sen_set.y(sn)).^2+(z - sen_set.z(sn)).^2);
    
    Ip = Vp.*R./100000.*factor(sn);
end

Ip = nanmean(Ip);

   


