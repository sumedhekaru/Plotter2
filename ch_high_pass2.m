function [t,y] = ch_high_pass2(t,y,cut_off_f)

% Creted from Sampath's file sent on 2015-06-15 to my email.

% Get the frequency of the data
fs = NaN;
i = 0;
while isnan(fs)
    i = i+1;
    fs = 1/(t(i+1)-t(i));
end


%% filter spec
N =4;  
Astop = 100;
cut_off_f;
F_stop = 2*cut_off_f/fs;

%% filter type
[z,p,k] =cheby2(N,Astop,F_stop,'high');
[s,L] = zp2sos(z,p,k);
H = dfilt.df2sos(s,L);
y = filtfilthd(H,y);
