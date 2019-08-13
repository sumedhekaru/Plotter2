function [width, ttrail, tlead, halfPeak] = fwhm(x,y)

% function width = fwhm(x,y)
%
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
%
% 2014-10-16 Sumedhe Karunarathne writen for positive ppulses

[maxy, ind] = max(y);
halfPeak = maxy/2;
width = NaN;
ttrail = NaN;
tlead = NaN;


% Walk backword until we find 50% level or begining of data
i = ind - 1;
while i > 0 && y(i) > halfPeak 
    i = i-1;    
end

if i == 0 % We have not go trhough 50% yet, can't find FWHM
    return
else
    ttrail = x(i) + (x(i+1)-x(i))*(halfPeak-y(i)) / (y(i+1)-y(i));    
end


% Walk foward until we find 50% peak or end of data
i = ind + 1;
L = length(y);

while i < L+1 && y(i) > halfPeak
    i = i+1;
end

if i == L+1 % We have not go trough 50% yet. can't find FWMH
    return
else
    tlead = x(i-1) + (x(i)-x(i-1))*(halfPeak-y(i-1)) / (y(i)-y(i-1));    
end

width = tlead - ttrail;


