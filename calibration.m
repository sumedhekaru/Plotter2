function handles=calibration(handles)
% THis function will set the calibration factor values
% Note that calibration values are in 'calibration_data.xlsx' file. If you
% want to change calibration values, you may want to go there and modify
% it. Also note that the first raw of the excel file is for manual
% calibration factors.

%try
% clc   
g=handles.g;
a=handles.sen_set;

% In the sensor_setting.mat file contain all the calibration data
%       column1: date yyyymmdd format
%       column 2:16 c
%       Ch1 gain respective to Ch2
%       Ch2 gain respective to KSC field mill
%       Ch3 gain respecive to Ch2
%       Please enter NaN if gain information is not availbale

%data = xlsread('calibration_data.xlsx');
data = csvread('calibration_data.csv',2,0);
try
    year  = g.YYYY{:};
    month = g.MM{:};
    date  = g.DD{:};
catch
    set(handles.cal_info,'String','Calibration data has not loaded yet but will be loaded when you plot')
    return
end
%date value
dvalue=str2double(sprintf('%s%s%s',year,month,date));

% Lets calculate date range average anyway because daily average may not
% availbale everyday.
%lower limit
if a.cal_type == 1
    if strcmp(year,'2010')
        a.cal_start = 20100701;
        a.cal_end   = 20100830;
    elseif strcmp(year,'2011')
        a.cal_start = 20110701;
        a.cal_end   = 20110830;
    else
        % do nothing
    end
end


lol =  sum(data(:,1) < a.cal_start)+1;
ul  =  sum(data(:,1) <= a.cal_end);

if lol < ul
    range_avg = nanmean(data(lol:ul,2:61));
else
    range_avg = data(lol:ul,2:61);
end

ch_str={};
% channel str
for i=1:20
    ch_str = [ ch_str sprintf('%s:1',a.sen_IDs{i})];
    ch_str = [ ch_str sprintf('%s:2',a.sen_IDs{i})];
    ch_str = [ ch_str sprintf('%s:3',a.sen_IDs{i})];
end
    
switch a.cal_type
    % Which type of calibrations are you looking for?
    case 1
        % Use daiy average value
        factor=NaN(1,60);
        
        msg_str1=sprintf('Calibration Type: DA ');
        msg_str2= sprintf('TRA used for (%i:%i):',a.cal_start,a.cal_end);  
        msg_str3 = 'MC used for:';
        
        % finding calibrations
        index = find(data(:,1)==dvalue);
        
        % found the day?
        if isempty(index) == 0
            factor=data(index,2:61);
        end
        
        
        % Now lets replace the values of NaN with Time Range Averge
        for i=1:60
            if isnan(factor(i))
                %if TRA not available lets use manual value
                if isnan(range_avg(i))
                    factor(i)=data(1,i+1);
                    if (handles.g.lpgraphs(i)==1 || handles.g.chgraphs(i)==1)
                        msg_str3 = [msg_str3 ' ' ch_str{i}];
                    end
                else
                    factor(i)=range_avg(i);
                    if (handles.g.lpgraphs(i)==1 || handles.g.chgraphs(i)==1)
                        msg_str2 =[msg_str2 ' ' ch_str{i}];
                    end
                end
                
            end
        end
        
    case 2
        
        msg_str1 = sprintf('Calibration Type: TRA Selected (%i:%i)',a.cal_start,a.cal_end);
        
        factor = range_avg;
        
        msg_str2 = 'MC is used for:';
        msg_str3 ='';
        
        for i=i:60
            if isnan(factor(i)) 
                factor(i)=data(1,i+1);
                if (handles.g.lpgraphs(i)==1 || handles.g.chgraphs(i)==1)
                    msg_str2 = [msg_str2 ' ' ch_str{i}];
                end
            end
        end
        
    case 3
        factor = data(1,2:61); 
        msg_str1 = 'Calibration Type: MC Selected';
        msg_str2 = '';
        msg_str3 = '';
end

if handles.sen_set.plot_calibrated == 0
    msg_str1 = 'Calibration Type: None';
    msg_str2 = '';
    msg_str3 = '';
end

%Lets find the real factors for ch1 and ch2 as they were just gains
%wrt ch2 intially
for i=1:20;
    factor(i*3-2) = factor(i*3-1)/factor(i*3-2);
    factor(i*3)   = factor(i*3-1)/factor(i*3);
end


set(handles.cal_info,'String',sprintf('%s\n%s\n%s',msg_str1,msg_str2,msg_str3))

g.factor = factor;
handles.g=g;
handles.sen_set=a;


