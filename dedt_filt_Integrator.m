function [t,y_subRegion,filtered] = dedt_filt_Integrator(t,y_dEdt,fp)


y_dEdt(isnan(y_dEdt)) = 0;

fs = 1/(t(2)-t(1));
i = 1;
while ~isnan(fs) && i < length(fs)
    i = i+1;
    fs = 1/(t(i+1)-t(i));
end

%% filter spec
N=4;  
Astop = 100;
F_stop = 2*fp/fs;

%% filter type
[z,p,k] =cheby2(N,Astop,F_stop,'high');
[s,L] = zp2sos(z,p,k);
H=dfilt.df2sos(s,L);
%% remove end effect of filter: 
M= impzlength(H)*2;
y_pad = padarray(y_dEdt,[0 M],'both');
y_filter=filtfilthd(H,y_pad);

y_subRegion=y_filter(M+1:M+length(y_dEdt));

filtered = 1;

