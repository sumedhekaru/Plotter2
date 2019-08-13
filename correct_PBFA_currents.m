function correct_PBFA_currents

% Old PBFA current calculations has been done using ch3 data. However, most
% of the time peak current came out were smaller than CGLSS estimation.
% Then we switch to using ch1 dat for current calculations. This function
% is used to calculate peak currents of ALREADY calculated PBFA points.
clc

%% User inputs
% Base folder up to PBFA (date will be added later)
a.bf = 'C:\data\2011-08-05 -- 2011-08-16\data\PBFA\';

% Make sure you run the plotter 
h=guidata(findall(0,'Tag','plotter2'));
sen_set = h.sen_set;
g = h.g;

%% Start the program
for i  = 0:47
   fn = sprintf('%s%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        a.bf,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},g.MM{:},g.DD{:},floor(i/2),mod(i,2)*30);
    
   % if the file available load it
   if exist(fn,'file')
        fid = fopen(fn);
        data=textscan(fid,'%f %f %f %f %f %f %f %s %f %f %f %f %f','HeaderLines',2);
        fclose(fid);
        
        L = length(data{:,1});
        
        for j = 1:L
            t0 = data{2}(j);  x0 = data{3}(j);  y0 = data{4}(j);  
            z0 = data{5}(j);  Ip = data{7}(j);
            
            Ipn = calculate_newIp(h,t0,x0,y0,z0)

        end
        
        return

       
   else
       fprintf('File not available: %s\n',fn)
   end
end

function Ipn = calculate_newIp(h,t0,x0,y0,z0)
g = h.g;
settings = h.sen_set;

% turn on all ch1 data
g.chgraphs = [1 0 0 1 0 0 1 0 0 0 0 0 0 0 0 ...
           1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
       
hms = sec2hhmmss(t0);
g.hh = str2double(hms(1:2));
g.mm = floor(str2double(hms(4:5))/5)*5;
g.ss = 0;

% Generate file names
[ch_fn h_fn] = generate_ch_fn(settings,g);

% Time shifts and calibration factors
tshift=settings.t_shift;
factor = g.factor;

figure
hold all
% Load data
for i = 1:3:30
    
    if ~isempty(ch_fn{i})
        [t y] = FA_Extract1(ch_fn{i},h_fn{i},t0-1e-3,t0+1e-3,tshift(i),settings,i);
        
        if ~isempty(t)
            
            [yF yH] = hill_tra(t,y,2000);
            yH = abs(yH);
            
            yPeakH = max(yH);
            yPeakF = min(yF);
            
            plot(t,yH)
            
            %         if abs(yPeakF*factor(i)) > 1
            %             Eps(sn) = yPeakH;
            %         end
            
        end
    end
end

Ipn  = 0;
       
       
function [ch_fn h_fn] = generate_ch_fn(settings,g)
% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

% checking witch graphs are on
ch_on=g.chgraphs;

for i=1:30
    if ch_on(i)
        % finding the file extention number
        ext=mod(i,3);
        if ext==0
            ext=3;
        end
        
        % Finding the stattion ID
        sid=settings.sen_IDs{ceil(i/3)};
        
        % finding the file name
        filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        
        hfilename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        
        
        % If file is not exist don't store the file name
        if ~exist(filename,'file')
            ch_fn{i}='';
            h_fn{i}='';
        else
            ch_fn{i}=filename;
            h_fn{i}=hfilename;
        end
    else
        ch_fn{i}='';
        h_fn{i}='';
    end
end

