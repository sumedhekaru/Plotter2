function show_pulse_intervals
% This function was written to plot pulse group times obtained from 
% 'get_pulse_interval' program.
%
%   Requirements: None
%
%   Usage : show_pulse_intervals
%           Once run, it will show the pulse grop time as vertical lines.
%           Times will be obtained from the file
%           'PulseIntervalData-20110814.csv'.
%
%   History:
%       2014-02-25 Created by Sumedhe Karunarathne
%       2014-03-05 Modified to include PulseIntervalData2 files


% get plotter2 data
h=guidata(findall(0,'Tag','plotter2'));
g = h.g;


%% Type 1 files
% Open file to read data
fn = sprintf('PulseIntervalData-%s%s%s.txt',g.YYYY{:},g.MM{:},g.DD{:});


if ~exist(fn,'file')
    % nothing to do if there is no file
    fprintf('file %s not found\n',fn);
else
    % if exist let's read the data
    fID = fopen(fn,'r');
    data = textscan(fID,'%s %f %f');
    fclose(fID);
    
    % sort data acording to time
    [t,I] = sort(data{2});
    
    % data{1} = data{1}(I);
    % data{2} = data{2}(I);
    % data{3} = data{3}(I);
    
    
    % Figure out times belongs to this plot
    tl = xlim;
    lol = sum(t < tl(1))+1;
    ul  = sum(t < tl(2));
    
    yl = ylim;
    
    t = t(lol:ul);
    
    % format before plot lines so that it really is a one line (easy to delete)
    
    L = length(t);
    Es = [];
    ts = [];
    
    for i = 1:L
        ts = [ts t(i) t(i) NaN];
        Es = [Es yl NaN];
    end
    
    plot(ts,Es,'r--')
end

%% Type 2 files
% Open file to read data
fn = sprintf('PulseIntervalData2-%s%s%s.txt',g.YYYY{:},g.MM{:},g.DD{:});


if ~exist(fn,'file')
    % nothing to do if there is no file
    fprintf('file %s not found\n',fn);
else
    % if exist let's read the data
    fID = fopen(fn,'r');
    data = textscan(fID,'%s %f %f');
    fclose(fID);
    
    % sort data acording to time
    [t,I] = sort(data{2});
    
    % data{1} = data{1}(I);
    % data{2} = data{2}(I);
    % data{3} = data{3}(I);
    
    
    % Figure out times belongs to this plot
    tl = xlim;
    lol = sum(t < tl(1))+1;
    ul  = sum(t < tl(2));
    
    yl = ylim;
    
    t = t(lol:ul);
    
    % format before plot lines so that it really is a one line (easy to delete)
    
    L = length(t);
    Es = [];
    ts = [];
    
    for i = 1:L
        ts = [ts t(i) t(i) NaN];
        Es = [Es yl NaN];
    end
    
    plot(ts,Es,'b--')
end

    

