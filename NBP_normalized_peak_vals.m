function NBP_normalized_peak_vals

% Load power curves

% file name for NBP data
fn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814-NBP_info2.xlsx';
sheet = 1;
xlRange = 'A3:AD307';
ndata = xlsread(fn, sheet, xlRange);

pbfat = ndata(:,10);
pbfax = ndata(:,11);
pbfay = ndata(:,12);
pbfaz = ndata(:,13);


% directory name for power data
%dn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814_plots\MeanPowers\';
dn = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\PowerCurves5MHz-2\';

% Directory to save stuffs
bf2 = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\PowerCurves5Mhz-3';

% DIrectory name to save

indx = ndata(:,1);

% get plotter2 data
h=guidata(findall(0,'Tag','plotter2'));
sen_set = h.sen_set;
g = h.g;


% To put in the title
pulse_kinds = { ...
    '0: Not claer'
    'A: Clean bipolar'
    'B: Pulses after neg overshoot'
    'C: Pulses befora and after neg overshoot'
    'D: Pulses before zero cross'
    'E: Pulses in rising edge'};

fID = fopen('C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\Types-20141030\20110814-type_info.txt','r');
typeD = textscan(fID,'%f %f');
types = typeD{2};
fclose(fID);


% Turn off all the plots
g.lpgraphs = zeros(1,60);
g.chgraphs = zeros(1,60);

% Channel to work on
ch = 3;

for type = 1:4
   
    %inds = [find(types == 6); find(types == 8)] ;
    switch type
        case 1
            inds = find(types == 1);
        case 2
            inds = [find(types == 2); find(types == 3)];
        case 3
            inds = find(types == 4);
        case 4
            inds = find(types == 5);
    end

    
    for i = 1:length(inds)

        % Power file name
        %pfn = sprintf('%sMeanPowersMeanPower-%3.3i.mat',dn,indx(i));
        pfn = sprintf('%s%4.4i-power_dist.fig',dn,indx(inds(i)));

        if exist(pfn,'file')
            fprintf('Working on %s\n',pfn)
            % b = open(pfn);
            fgh = openfig(pfn,'reuse','invisible');
            d = guidata(fgh);
            close(fgh);
            

            % Let's turn on the just the channels we need
            g.chgraphs(d.sns*3 - (3-ch)) = 1;
            
            % Get the location of this NBP
            ind = inds(i);
            t0 = pbfat(ind);
            x0 = pbfax(ind);
            y0 = pbfay(ind);
            z0 = pbfaz(ind);
           
            % Get the positive peak for these plots
            %get_positive_peak(ind,g,sen_set,t0,x0,y0,z0,ch)

            % Get the power curve
            power_calculation(d.sns,ch,g,sen_set,0,x0,y0,z0,t0)
            
            fg = gcf;
             % save the figure
            set(fg, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperPosition', [0 0 40 25])
            saveas(fg,sprintf('%s/%3.3i-power_dist.fig',...
                bf2,ind))
            saveas(fg,sprintf('%s/%3.3i-power_dist.png',...
                bf2,ind))
            
            delete(fg)


         
        else
            fprintf('File not found %s\n',pfn)
        end
    end   
    
end 


function get_positive_peak(ind,g,sen_set,t0,x0,y0,z0,ch)



% Check unavailable files if you want
% fprintf('%s\n',a.absant_fn{:})

% TIme shifts
tshift=sen_set.t_shift;


% Load ch data
fgh = figure('units','normalized','outerposition',[0.05 0.05 0.95 0.95],'visible','off');
hold all

R = sqrt((sen_set.x-x0).^2 + (sen_set.y - y0).^2 + (sen_set.z - z0).^2);
D = sqrt((sen_set.x-x0).^2 + (sen_set.y - y0).^2)/1000;

wbh = waitbar(0,'Please wait..','name','plotting');

vOffset = 0;

[Dn, sortIn] = sort(D,'descend');

wbcounter = 0;

colors = [
         0         0    1.0000
    1.0000         0         0
         0    1.0000         0
         0         0    0.1724
    1.0000    0.1034    0.7241
    1.0000    0.8276         0
         0    0.3448         0
    0.5172    0.5172    1.0000
    0.6207    0.3103    0.2759
         0    1.0000    0.7586
         0    0.5172    0.5862];

% figure data
figD.t0 = t0;
figD.x0 = x0;
figD.y0 = y0;
figD.z0 = z0;
figD.index = ind;
figD.sns = [];
figD.Rs = [];
figD.Ds = [];
figD.normPeaks = [];
figD.Peaks = [];
figD.lg = {};

for ind1 = sortIn
    
    i = ind1*3 - (3-ch);
    
    wbcounter = wbcounter + 1;
    waitbar(wbcounter/60,wbh)
    
    % Generate file names
    g.hh = floor(t0/3600);
    g.mm = floor(((t0 - 3600*g.hh)/60)/5)*5;
    a = file_names_generator(g,sen_set);
    
    if ~strcmp(a.ch_fn{i},'') && D(ceil(i/3)) > 30
        % Arrival time
        sn = ceil(i/3);
        
        % Lets load 75us from both sides of arrival time
        g.t1 = t0 + R(sn)/3.0e8 - 25e-6;
        g.t2 = t0 + R(sn)/3.0e8 + 75e-6;

        [t, y, ch_freq_str] = FA_Extract1(a.ch_fn{i},a.h_fn{i},g.t1,g.t2,tshift(i),sen_set,i);
        t = (t - g.t1)*1e6;
        if ~isempty(t) && range(y)>0.001
            
            y = y*g.factor(i);
            offset = nanmean(y(1:50));
            yn = y-offset;
            
            % real Peak
            realP = max(yn);
            
            % Normalize
            yn = yn*R(sn)/1000/100;
            
            [peaky, peakI] = max(yn); 
            plot(t,yn,'color',colors(ind1,:),'LineWidth',0.5)
            pH = plot(t(peakI),peaky,'ro');
            set(get(get(pH,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            xlim([0, g.t2-g.t1]*1e6)
            str = [a.ch_legend{i}(1:3) ' (' sprintf('%0.1f',D(sn)) ' km)'];
            ylabel('E-change (V/m)')
            xlabel('Time (\mus)')
            %text(60,nanmean(yn(end-50:end)),str,'HorizontalAlignment','left','VerticalAlignment','bottom')
            
            vOffset = vOffset + 10;

            % Add data to figure handle
            figD.sns = [figD.sns sn];
            figD.Rs = [figD.Rs R(sn)/1000];
            figD.Ds = [figD.Ds D(sn)];
            figD.normPeaks = [figD.normPeaks peaky];
            figD.Peaks = [figD.Peaks realP ];
            figD.lg = [figD.lg str];
        end
    end
    
end

legend(figD.lg)

guidata(gcf,figD)
box on

title(sprintf('Index = %3.3i     %0.6f s   (%0.1f, %0.1f, %0.1f) km',ind,t0,x0/1000,y0/1000,z0/1000))


delete(wbh)
set(fgh,'visible','on')
