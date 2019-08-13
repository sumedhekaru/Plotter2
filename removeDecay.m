function Vnn = removeDecay(vn,tau,dt)

% Remove decay according to the tau dt


indx = isnan(vn);
vn(indx) = 0;

Vnns = (vn - vn*exp(-dt/tau));

Vnn = vn+cumsum(Vnns);
vn(indx) = NaN;