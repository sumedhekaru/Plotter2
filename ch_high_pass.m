function [t,y,filtered] = ch_high_pass(t,y,f)

% This function will do high-pass butter worth filtering and if it is
% successful filtered = 1. If not filtered = 0;

filtered = 0;

y_org = y;

y(isnan(y))=0;
Fs = NaN;

% filter out low frequencies
L = length(t);

if L == 0
    % disp('There are no data to filter')
    return
end


i = 1;
while isnan(Fs) && i < L
    Fs = 1/(t(i+1)-t(i));
    i = i+1;
end

[z,p] = butter(5,f/(Fs/2),'high'); % Create a High-pass butterworth filter;

try
    yF = filtfilt(z,p,y);    % filter the data.
catch
     return
end

if isnan(yF)
    filtered = 2;
    % Filtereing doesn't work - too smal frequency
    y = y_org;
else    
    filtered = 1;    
    y = yF;
end

                