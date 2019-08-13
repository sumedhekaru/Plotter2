function NBP_power_curve_generator

% This function is written to generate power curves for NBPs we found on
% August 14 2011. The prgram will load time of NBPs and then generate power
% curves.
%
% After loading press
%   q - to quit
%   r - remove sensor numbers, and then enter sensor numbers and press
%   enter (For example, enter 175, to remove sensor numbers 1,7, and 5.
%   n - to go to next plot (ignore the current plot)
%   d - Done! record data in the file


% Open the NBP file
data = xlsread('C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814-NBP_info2.xlsx');

% Starting index acoording to column 1
index = 1;

% base foloder for saving plots
bf = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814\PowerCurves5MHz-2\';

% data saving file name
dsfn = [bf '20110814-power_info.txt'];

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
        x0 = data(index+2,11);
        y0 = data(index+2,12);
        z0 = data(index+2,13);
        g.hh = floor(t/3600);
        g.mm = floor(((t - 3600*g.hh)/60)/5)*5;
        g.t1 = t - 0.0006;
        g.t2 = t + 0.0006;
        
        sns = [1 2 3 6 8 9 10];
        ch = 3;
        power_calculation(sns,ch,g,sen_set,0,x0,y0,z0)
        
        fg = gcf;
        
        
        while ~ismember(cc,{'q','n','d'})
            figure(fg)
            cc = get(fg,'CurrentCharacter');
            pause(0.1)
            
            
            % Go to next plot means we do not record the data for this plot
            if strcmp(cc,'n')
                
                fprintf(fID,'%10.3i\t%10.3e\n',...
                    index,NaN);
                
                
            elseif strcmp(cc,'r')
                answer = inputdlg_modi('Enter the sensor numbers to remove');
                answer = answer{1};
                L = length(answer);
                
                new_sns = zeros(1,L);
                
                for i = 1:length(answer);
                    new_sns(i) = str2double(answer(i));
                end
                
                % Fix sn = 10
                ind0 = find(new_sns == 0);
                
                if ~isempty(ind0)
                    new_sns(ind0) = 10;
                    new_sns(ind0-1) = [];
                end
                
                
                sns = setdiff(sns,new_sns);
                
                delete(fg);
                
                power_calculation(sns,ch,g,sen_set,0,x0,y0,z0)
                fg = gcf;
                cc = get(fg,'CurrentCharacter');
                
            elseif strcmp(cc,'u')
                answer = inputdlg_modi('Enter the sensor numbers to use');
                answer = answer{1};
                L = length(answer);
                
                new_sns = zeros(1,L);
                
                for i = 1:length(answer);
                    new_sns(i) = str2double(answer(i));
                end
                                
                % Fix sn = 10
                ind0 = find(new_sns == 0);
                
                if ~isempty(ind0)
                    new_sns(ind0) = 10;
                    new_sns(ind0-1) = [];
                end
                
                sns = new_sns;
                
                delete(fg);
                
                power_calculation(sns,ch,g,sen_set,0,x0,y0,z0)
                fg = gcf;
                cc = get(fg,'CurrentCharacter');
                
                
                % Done (means record the current data)
            elseif strcmp(cc,'d')
                
                % Get figure data
                figD = guidata(fg);
                
                % Print info to the file
                fprintf(fID,'%10.3i\t%10.3e\n',...
                    index,figD.totalMeanP);
                % save the figure
                saveas(fg,sprintf('%s%4.4i-power_dist.fig',bf,index))
                saveas(fg,sprintf('%s%4.4i-power_dist.png',bf,index))
                
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



















