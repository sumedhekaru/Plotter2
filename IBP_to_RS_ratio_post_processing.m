function IBP_to_RS_ratio_post_processing
clc
%% key preses
% q -- quit the program
% s -- skip the file
% p -- done with the current figure --> move to the selected folder
% z -- zoom (select two points where you want to zoom, any plot)
% r -- reset zoom (to original view)
% c -- modify RS peak
% x -- modify IBP peak
% w -- estimate background noice
% g -- pure CG (means no IBP can be sean in this figure)

% v -- Nice IBP (please press this if you see nice IBPs with no subpulses)
% b -- Pulse bursts (if you a pulse burst please press this)
% n -- Narrow Bipolar
% d -- dart leader 


%% User inputs
% base folder
bf = 'C:\Users\Sumedhe\Desktop\IBP_to_RS_ratio_auto2\';

% Distance limit (in km)
r1 = 0;
r2 = 100;

% Time region interested
t1 = 78075;
t2 = 79500;

% sensors used
sns = [1 2 3 6 7 8 9 10];

% Save file name for data
sFn = 'IBP_to_RA_ratio_selected.txt';

%% Starting programming

% Create required folders
if ~exist([bf 'dartLeader'],'dir'); mkdir([bf 'dartLeader']); end;
if ~exist([bf 'NBP'],'dir'); mkdir([bf 'NBP']); end;
if ~exist([bf 'niceIBP'],'dir'); mkdir([bf 'niceIBP']); end;
if ~exist([bf 'pulseBursts'],'dir'); mkdir([bf 'pulseBursts']); end;
if ~exist([bf 'selected'],'dir'); mkdir([bf 'selected']); end;


Lsns = length(sns);

% read the info file
fID = fopen([bf 'info.txt'],'r');
data = textscan(fID,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',2);
fclose(fID);

% Columns 
% 1-Index       2-RS-peak	3-IBP-np	4-IBP-pp	5-Noice	
% 6-Delta-t     7-tCGLSS	8-xCGLSS	9-yCGLSS	10-Range	
% 11-RSp(Norm)	12-IBP-np(Norm)	13-IBP-pp(Norm)	14-Noice(Norm)
% Sort according to distance from LDAR2 origin

% Sort data according to distance
dist = sqrt(data{:,8}.^2+data{:,9}.^2);
[dist, ix] = sort(dist);


% Choose the only values within r1 and r2
lol = nnz(dist < r1) + 1;
ul  = nnz(dist < r2);


for i=1:14;   dataN{i} = data{i}(ix(lol:ul),:);  end

data = dataN;
clear dataN;

% Sort  according to time
time = data{7};
[time, ix] = sort(time);


% Choose the only values within r1 and r2
lol = nnz(time < t1) + 1;
ul  = nnz(time < t2);


for i=1:14;   dataN{i} = data{i}(ix(lol:ul),:);  end

% fileNames = ls([bf '*.fig']);

fileNames = dataN{1};
L = length(fileNames);

for i = 1:L
    fname = [fileNames{i} '.fig'];
    %fname = fileNames(i,:);

    key = start_working(bf,fname,sFn);
    
    % is user need to quit?
    if strcmp(key,'q');
        disp('User stopped data processing.')
        return;
    end
end
    




%% Funtions
function key = start_working(bf,fname,sFn)

    %fprintf('\n********************************************************\n')
    fprintf('\n\t\t\tworking on the file %s\n',fname)
    fprintf('********************************************************\n')
    openfig([bf fname],'new','invisible')
    fg = gcf;
    
    % handles to subplots
    h1 = subplot(2,3,1:3);
    h2 = subplot(2,3,4);
    h3 = subplot(2,3,5);
    h4 = subplot(2,3,6);   
    
    
    % Current title
    tit = get(get(h1,'Title'),'String');

%     %screen2png(fg,[bf fname(1:end-3) 'png'])
%     key = 'o';
%     delete(fg)
%     sn = str2double(fname(7:8));
%     save_to_file(bf,'info.txt',tit,sn)
%     return
        
    % limits of figures
    xl1 = get(h1,'xlim');
    xl2 = get(h2,'xlim');
    xl3 = get(h3,'xlim');
    xl4 = get(h4,'xlim');
    
    % lines of first subplots
    h = findobj(h1,'Type','line');  
    

    
    key = 'o';

    while ~ismember(key,{'q','s'})
        figure(fg)
        try
            waitforbuttonpress;
        catch
            disp('User has closed the figure.')
            key = 'q';
            return
        end
            
        key = get(fg,'CurrentCharacter');

        switch key
            % User needs to skip this file'
            case 's' ; disp('skipped.'); delete(fg);

            % User has done with this figuredisp('done');
            case 'p' ;
                msgh = msgbox('Please wait...');
                % Save new figure
                saveas(fg,[bf 'selected\' fname])
                screen2png(fg,[bf 'selected\' fname(1:end-3) 'png'])
                
                % delete old ones
                delete([bf fname],[bf fname(1:end-3) 'png'])   
                
                % Sensor used
                sn = str2double(fname(7:8));
                
                % Save title info to the file
                save_to_file([bf 'selected/'],sFn,tit,sn)
                
                fprintf('Done! - Moved to the ''selected'' folder\n');                
                delete(fg);
                key = 's';
                delete(msgh);
            % Horizontal zoom
            case 'z'
                disp('zooming...')
                % get the location where user wants to zoom
                [t1, ~]  = ginput(1);
                [t2, ~]  = ginput(1);
                
                try
                    xlim([t1 t2])
                catch
                    disp('zooming error occured')
                end
                
                % make the key back to default
                set(fg,'CurrentCharacter','o')
                key = 'o';
            
            % reset zoom figure
            case 'r'
                set(h1,'xlim',xl1);
                set(h2,'xlim',xl2);
                set(h3,'xlim',xl3);
                set(h4,'xlim',xl4);
                
                % make the key back to default
                set(fg,'CurrentCharacter','o')
                key = 'o';                
                disp('Resat.')
            
            % modify CG peak
            case 'c'
                disp('Modifying RS peak...')
                disp('Please click on RS peak.')
                [t1, ~]  = ginput(1);
                 
                try
                    % Get all the data
                    t = get(h(6),'Xdata');
                    y = get(h(6),'Ydata');
                    
                    
                    lol = nnz(t<t1-75e-6)+1;
                    ul  = nnz(t<=t1+75e-6) ;
                    
                    % Obtaining RS peak
                    tRS = t(lol:ul);
                    yRS = y(lol:ul);
                    
                    [RS_minY , RS_min_ind] = min(yRS);
                    RS_mint = tRS(RS_min_ind);
                    if RS_min_ind > 10; RS_minYR =  range(yRS(RS_min_ind-10:RS_min_ind));
                    else                RS_minYR =  range(yRS(1:RS_min_ind)); end;
                    
                    % delete old RS peak
                    delete(h(5))
                    
                    % add new RS peak
                    hold all
                    h(5) = plot(RS_mint,RS_minY,'ro');
                    
                    % New title
                    tit(24:30) = sprintf('%07.3f',RS_minYR);
                    title(tit)
                    
                    % Rs plot
                    subplot(h4); cla(h4);
                    plot(tRS,yRS,'b')
                    plot(RS_mint,RS_minY,'ro')
                    plot([RS_mint, RS_mint],[RS_minY,RS_minY + RS_minYR],'r')
                    xl4 = [min(tRS), max(tRS)];
                    xlim(xl4);
                catch
                    disp('Uknown error occured. It is better to skip this file.')
                end
                                                   
                % make the key back to defaulict
                set(fg,'CurrentCharacter','o')
                key = 'o'; 
                
            % modify IBP peak
            case 'x'
                disp('Modifying IBP peak...')
                disp('Please click on IBP peak.')
                [t1, ~]  = ginput(1);
                 
                 try
                    % Get all the data
                    t = get(h(6),'Xdata');
                    y = get(h(6),'Ydata');
                    
                    lol = nnz(t<t1-75e-6)+1;
                    ul  = nnz(t<=t1+75e-6) ;
                    
                    % Obtaining RS peak
                    tIBP = t(lol:ul);
                    yIBP = y(lol:ul);
                    
                    
                    [IBP_minY , IBP_min_ind] = min(yIBP);
                    IBP_mint = tIBP(IBP_min_ind);
                    if IBP_min_ind > 20;  IBP_minYR =  range(yIBP(IBP_min_ind-20:IBP_min_ind));
                    else                  IBP_minYR =  range(yIBP(1:IBP_min_ind));   end
                    
                    % Obtaining IBP peak to peak
                    if IBP_min_ind < 50;      ind1 = 1;
                    else ind1 = IBP_min_ind - 50;  end;
                    
                    LIBP = length(tIBP);
                    if IBP_min_ind + 50 > LIBP;      ind2 = LIBP;
                    else ind2 = IBP_min_ind + 50;  end;
                    
                    
                    if IBP_min_ind + 20 > LIBP;      ind3 = LIBP;
                    else ind3 = IBP_min_ind + 20;  end;
                    
                    tP_IBP = tIBP(IBP_min_ind:ind3);
                    yP_IBP = yIBP(IBP_min_ind:ind3);
                    [IBP_maxY , IBP_max_ind] = max(yP_IBP);
                    IBP_maxt = tP_IBP(IBP_max_ind);
                    
                    % delete old IBP peak
                    delete(h(4))
                    
                    % add new IBP peak
                    hold all
                    h(4) = plot(IBP_mint,IBP_minY,'ro');
                    
                    % New title
                    tit(47:53) = sprintf('%07.3f',IBP_minYR);
                    tit(70:76) = sprintf('%07.3f',range([IBP_minY,IBP_maxY]));
                    title(tit)
                    
                    % IBP plot
                    subplot(h3); cla; hold all
                    plot(tIBP(ind1:ind2),yIBP(ind1:ind2),'b')
                    plot(IBP_mint,IBP_minY,'ro')
                    plot([IBP_mint, IBP_mint],[IBP_minY,IBP_minY + IBP_minYR],'r')
                    plot(IBP_maxt,IBP_maxY,'ro')
                    plot([IBP_maxt, IBP_maxt],[IBP_minY,IBP_maxY],'r')
                    xl3 = [min(tIBP), max(tIBP)];
                    xlim(xl3);
                    
                catch
                    disp('Uknown error occured. It is better to skip this file')
                end
                                                   
                % make the key back to default
                set(fg,'CurrentCharacter','o')
                key = 'o'; 
            
            % Estimate background noice
            case 'w'
                disp('Estimating background noice...')
                disp('Please select two points to select a noice reagion.')
                
                [~, topv]  = ginput(1);
                [~, botv]  = ginput(1);              
                
                
                
                % add new noice to the title
                tit(92:98) = sprintf('%07.3f',topv-botv);
                
                subplot(h1); hold all;
                title(tit)
                
                % Redraw noice levels in the plot
                delete(h(1))
                delete(h(2))
                delete(h(3))
                h(1) = plot(xlim,[topv topv],'r');
                h(2) = plot(xlim,[botv botv],'r');
                h(3) = plot(xlim,[botv botv],'r'); % draq one more line to cover h3
                
                               
                % make the key back to default
                set(fg,'CurrentCharacter','o')
                key = 'o'; 
                
            % If no IBPs found in the current figure, lets make those as nans    
            case 'g'
                
                % Make IBP amplitudes to NaNs
                subplot(h1)
                tit(47:53) = sprintf('%07.3f',NaN);
                tit(70:76) = sprintf('%07.3f',NaN);
                title(tit)
                
                % IBP subplot should be irased as it is not matter any more
                cla(h3)
                
                % IBP peak in the figure should be deleted
                try delete(h(4)); end
                
                % make the key back to default
                set(fg,'CurrentCharacter','o')
                key = 'o';
                
            % a nice IBP has deteced, save data for future use    
            case 'v'
                msgh = msgbox('Please wait...');
                % Save new figure
                saveas(fg,[bf 'niceIBP\' fname])
                screen2png(fg,[bf 'niceIBP\' fname(1:end-3) 'png'])                
                
                % make the key back to default
                set(fg,'CurrentCharacter','o')
                key = 'o';
                delete(msgh)
             
            % a pulse burst was detected
            case 'b'
                msgh = msgbox('Please wait...');
                % Save new figure
                saveas(fg,[bf 'pulseBursts\' fname])
                screen2png(fg,[bf 'pulseBursts\' fname(1:end-3) 'png'])                
                
                % make the key back to default
                set(fg,'CurrentCharacter','o')
                key = 'o';
                delete(msgh)
            case 'n'
                msgh = msgbox('Please wait...');
                % Save new figure
                saveas(fg,[bf 'NBP\' fname])
                screen2png(fg,[bf 'NBP\' fname(1:end-3) 'png'])                
                
                % make the key back to default
                set(fg,'CurrentCharacter','o')
                key = 'o';
                delete(msgh)
                
            case 'd'
                % Save new figure
                msgh = msgbox('Please wait...');
                saveas(fg,[bf 'dartLeader\' fname])
                screen2png(fg,[bf 'dartLeader\' fname(1:end-3) 'png'])                
                
                % make the key back to default
                set(fg,'CurrentCharacter','o')
                key = 'o';
                delete(msgh)
                
            
            % User needs to quit the program
            case 'q'; delete(fg);                
                
        end

    end
    
function save_to_file(bf,sFn,tit,sn)

sns_x =  1.0e+04 *[ -0.7524, 0.3372, -0.3149, -0.4266, -1.1658, -0.2701, -2.7640, ...
                 -6.0091, 0.1825, -5.7394, -2.0637];
             
sns_y =  1.0e+04 *[1.6555, 0.4446, -0.6838, 0.9545, -1.7020, 0.2631, 4.9254, ...
                    -3.3983, -5.3008, 1.1923, 0.1569];

fileN = [bf sFn];

% Write the header if not exist
if ~exist(fileN,'file')
    fID = fopen(fileN,'a+');   
    fprintf(fID,'Index\tRS-peak\tIBP-np\tIBP-pp\tNoice\tDelta-t\ttCGLSS\txCGLSS\tyCGLSS\tRange\tRSp(Norm)\tIBP-np(Norm)\tIBP-pp(Norm)\tNoice(Norm)\n\n');   
else
    fID = fopen(fileN,'a+'); 
end

CGx = str2double(tit(144:150));
CGy = str2double(tit(161:167));
R = sqrt((CGx - sns_x(sn)/1000)^2+(CGy - sns_y(sn)/1000)^2);

RS_p = str2double(tit(24:30))*R/100;
IBP_np = str2double(tit(47:53))*R/100;
IBP_pp = str2double(tit(70:76))*R/100;


fprintf(fID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%7.1f\t%07.3f\t%07.3f\t%07.3f\n', tit(8:15),tit(24:30), tit(47:53), tit(70:76), ...
    tit(92:98),tit(107:112), tit(123:134),tit(144:150), tit(161:167),R,RS_p,IBP_np,IBP_pp);
   
fclose(fID);


