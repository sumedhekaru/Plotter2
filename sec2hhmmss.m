function [hhmmss hh mm ss dss]=sec2hhmmss(sec_from_midnight)
% convert hhmmss format from second from midnight
%
% Usage Ex. hhmmss = sec2hhmmss(7123)
% will return 
%
% hhmmss =
%
%    01:58:43 

    s=sec_from_midnight;
    
    hh = floor(s/3600);
    mm = floor(mod((s/60), 60));
    ss = floor(mod(s,60));
    dss = num2str(mod(s,60)-ss);
    %ms = mod(s,1)*1e3;

    hhmmss=sprintf('%02d:%02d:%02d%s',hh,mm,ss,dss(2:end));
    dss = str2double(dss);