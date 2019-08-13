function data2text(fn)
% This function will convert ch1 data to text optional input is the file
% name you want to put. Input for the data coming from plotter2 and it
% should be run before you run this. 
%
%   Author: Sumedhe Karunarathne
% 
% Modification history
%   2014-04-02 Initial program written

%% Get plotter 2 input
h=guidata(findall(0,'Tag','plotter2'));
g = h.g;
sen_set = h.sen_set;


% selected sensor
i = find(g.chgraphs == 1);
i = i(1);

sns = ceil(i/3);


if nargin < 1
    
    % Get the desktop folder
    dsktop = winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop');
    
    [hms hh mm ss] = sec2hhmmss(g.t1);
    
    % saving name
    fn = sprintf('%s/%s_%s%s%s_%2.2i%2.2i%2.2i.txt',...
        dsktop,sen_set.sen_IDs{sns},g.YYYY{:},g.MM{:},g.DD{:},hh,mm,ss);
end

%% load data
% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

factor = g.factor;
tshift=sen_set.t_shift;
vshift=sen_set.vshift;


ext=mod(i,3);
if ext==0
    ext=3;
end

% Finding the stattion ID
sid=sen_set.sen_IDs{ceil(i/3)};

% finding the file name
filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
    bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);

hfilename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
    bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);


% If file is not exist we don't need to do anything
if ~exist(filename,'file') || ~exist(hfilename,'file')
    fprintf('%s does not exist. Returning...',filename)
    fprintf('%s does not exist. Returning...',hfilename)  
    return
end

[t y ch_freq_str] = FA_Extract1(filename,hfilename,g.t1,g.t2,tshift(i),sen_set,i);
y = y*factor(i)+vshift(i) ;

figure
plot(t,y)

% Write data in to a text file
fID = fopen(fn,'w');
for i = 1:length(t)
    fprintf(fID,'%13.7f\t%5.4f\n',t(i),y(i));
end

fclose(fID);
    