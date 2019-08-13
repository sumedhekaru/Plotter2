function [tn,En]=fieldMillExtract3(fn,t1,t2,tshift)


% Read data from the file
fid = fopen(fn);
data=textscan(fid,'%f %f:%f:%f %f','HeaderLines',0);
fclose(fid);

%data=cell2mat(data);
%t=data(:,2)*3600+data(:,3)*60+data(:,4);
%E=data(:,5);


t=data{2}*3600+data{3}*60+data{4}+tshift;
E=data{5};




% find data point position as
% tn, En will only have the date between time interval given
lol=size(t,1)- nnz(t>t1)+1;             % Index of the lower matrix element
ul=nnz(t<t2);                           % Index of the upper matrix element  


tn=t(lol:ul);                           % t data for given time interval
En=E(lol:ul);                           % E data for given time intervel

% figure
% plot(tn,En)
% size(tn)
% size(En)
% 
% format long
% tn
% En