function NBP_catogorizer

% This function is written to catogorize NBPs according to their types.
%
% After loading press
%   q - to quit
%   0 - go to next (Can't determine the type) (type 0)
%   1 - Catogory 1 - Clean bipolar pulses
%   2 - Catogory 2 - Bipolar w extra pulses after overshoot peak
%   3 - Catogory 3 - extra pulses befor and afer overshoot
%   4 - Catogory 4 - Bipolar with extra pulse between first peak and zero cross
%   5 - Catogory 5 - extra pulses before leading peak

% Open the NBP file
data = xlsread('C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814-NBP_info2.xlsx');

% Starting index acoording to column 1 
% index 1 will process 1st item
index = 206;

% base foloder for saving plots
bf = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\Types-20141030\';

% data saving file name
dsfn = [bf '20110814-type_info.txt'];

%% Start the program

% Prepare writing file
if ~exist(dsfn,'file')
    fID = fopen(dsfn,'a+');
    fprintf(fID,'Index   TotMeanPower\n');
else
    fID = fopen(dsfn,'a+');
end

% Load data
% get plotter2 data
h=guidata(findall(0,'Tag','plotter2'));
sen_set = h.sen_set;
g = h.g;

% Turn off all plots
g.lpgraphs = zeros(1,60);
g.chgraphs = zeros(1,60);

cc = '';

% Seconds are most probably 0
g.ss = 0;

while ~strcmp(cc,'q')
    
    % to ingnore header we need to add 2 lines to raws
    t = data(index+2,10);
    
    % If t is NaN, there is no PBFA for this NBP and let's go to the next
    % one
    if isnan(t)
        fprintf(fID,'%10.3i\t%10.3e\n',...
            index,NaN);
    else
        % PBFA location
        x0 = data(index+2,11);
        y0 = data(index+2,12);
        z0 = data(index+2,13);
        
        % Time input for plotter
        g.hh = floor(t/3600);
        g.mm = floor(((t - 3600*g.hh)/60)/5)*5;
        g.t1 = t - 0.0006;
        g.t2 = t + 0.0006;
        
        % Sensors to use
        sns = [1 2 3 6 8 9 10];
        
        % Channels to use
        ch = 3;
        
        % Let's turn on the just the channels we need
        g.chgraphs(sns*3 - (3-ch)) = 1;
        
        % Generate figure
        gen_multi_figure(index,g,sen_set,t,x0,y0,z0)
        

        fg = gcf;
        
        
        while ~ismember(cc,{'q','0','1','2','3','4','5'})
            figure(fg)
            cc = get(fg,'CurrentCharacter');
            pause(0.1)
            
            
            % Go to next plot means we do not record the data for this plot
            if ismember(cc,{'1','2','3','4','5','0'})
                
                wave_type = str2double(cc);
                
                % Print info to the file
                fprintf(fID,'%10.3i\t%10.3i\n',...
                    index,wave_type);
                
                % save the figure
                set(fg, 'PaperUnits', 'centimeters');
                set(gcf, 'PaperPosition', [0 0 40 25])
                saveas(fg,sprintf('%s/%1.1i/%4.4i-%1.1i-power_dist.fig',...
                    bf,wave_type,index,wave_type))
                saveas(fg,sprintf('%s/%1.1i/%4.4i-%1.1i-power_dist.png',...
                  bf,wave_type,index,wave_type))
            
            end
            
        end
        
        % reset the character
        if ~strcmp(cc,'q')
            cc = '';
        end
        
        % Delete the current figure
        delete(fg)
    end
    
    index = index + 1;
    
end


fclose(fID);
disp('done')


function gen_multi_figure(ind,g,sen_set,t0,x0,y0,z0)

% Generate file names
a = file_names_generator(g,sen_set);

% Check unavailable files if you want
% fprintf('%s\n',a.absant_fn{:})

% TIme shifts
tshift=sen_set.t_shift;

% Load ch data

fgh = figure('units','normalized','outerposition',[0.05 0.05 0.95 0.95],'visible','off');

R = sqrt((sen_set.x-x0).^2 + (sen_set.y - y0).^2 + (sen_set.z - z0).^2);
D = sqrt((sen_set.x-x0).^2 + (sen_set.y - y0).^2)/1000;

spot = [1 1 1 2 2 2 3 3 3 0 0 0 0 0 0 5 5 5 6 6 6 7 7 7 9 9 9 10 10 10 11 11 11];
wbh = waitbar(0,'Please wait..','name','plotting');

for i = 1:60
    waitbar(i/60)
    if ~strcmp(a.ch_fn{i},'')
        % Arrival time
        sn = ceil(i/3);
        
        % Lets load 75us from both sides of arrival time
        g.t1 = t0 + R(sn)/3.0e8 - 75e-6;
        g.t2 = t0 + R(sn)/3.0e8 + 75e-6;
        
        
        [t, y, ch_freq_str] = FA_Extract1(a.ch_fn{i},a.h_fn{i},g.t1,g.t2,tshift(i),sen_set,i);
        t = (t - g.t1)*1e6;
        if ~isempty(t) && range(y)>0.001
            
            y = y*g.factor(i);
            offset = nanmean(y(1:10));
            
            subplot(3,4,spot(i))
            box on
            plot(t,y-offset)
            xlim([0, g.t2-g.t1]*1e6)
            title([a.ch_legend{i} ch_freq_str '     D = ' sprintf('%0.1f',D(sn)) ' km'])
        end
    end
end


% Let's plot xy 
subplot(3,4,[4 8 12]);
florida_map
xlim([-120, 120])
ylim([-200, 200])
daspect([ 1 1 1])
box on
plot(x0/1000,y0/1000,'ro','markerfacecolor','r')

plot(sen_set.x([1 2 3 6 7 8 9 10 11])/1000,...
     sen_set.y([1 2 3 6 7 8 9 10 11])/1000,...
      'b*','markerfacecolor','b')

h = annotation('textbox',[0.75,0.8,0.1,0.1],'string',sprintf('Index = %3.3i',ind));
set(h,'fontsize',20,'lineStyle','none')

% Display the distance from the zero as the title
title(sprintf('Index = %i       R_0 = %0.1f',ind,sqrt(x0^2+y0^2+z0^2)/1000))

delete(wbh)
set(fgh,'visible','on')












