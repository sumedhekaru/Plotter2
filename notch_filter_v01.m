function yn = notch_filter_v01(t,y,f0,bw,at)

fs = 1/(t(2) - t(1));
i = 2;
while isnan(fs)
    fs = 1/(t(i+1)-t(i));
end

wo = f0/(fs/2);
bw = bw/(fs/2);
[b, a] = iirnotch(wo,bw,at);
yn = filtfilthd(b,a,y);

