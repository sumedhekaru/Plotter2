function find_peak_current

% Import plotter data
h=guidata(findall(0,'Tag','plotter2'));

opt.textOut = false;
opt.graphOut = false;
data = find_max_dE(h,opt);

Ep  = data.Ep;
%sn  = data.sn;
R   = data.R/1000;
L   = length(R);
Ip  = nan(1,L);

for i = 1:L
    if R(i) > 20 && R(i) < 70
        % Best fit parmeters from Nadee for differnt distances
        if      R(i) >= 40;    n = 0.71;   m = 0.1915;     c = 0.0700;
        elseif  R(i) >= 30;    n = 0.71;   m = 0.1915;     c = 0.0700;
        elseif  R(i) >= 20;    n = 0.70;   m = 0.1870;     c = 6.8232;
        else                   n = 0.88;   m = 0.0594;     c = -4.8609;
        end
        
        n = 0.71;   m = 0.1915;     c = 0.0700;
        
        Ip(i)= m*(Ep(i)*R(i)^n) + c;
    else
        %do nothing
    end
end

Ip'

mu = nanmean(Ip)

