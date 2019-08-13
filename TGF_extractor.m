function TGFd = TGF_extractor(fn)

%fn = 'C:\Users\Sumedhe\Desktop\Fermi TGF near KSC For Marshall\Fermi_TGF_20110801_0257UT\2011-08-01_02-57-17_Duke_VLF.mat';

data = open(fn);
s = data.s;


% Determine the sensor type
SF = s(1).sampling_freq;

if SF == 1000000
    Figure_sensor = 'LF';
elseif SF == 100000
    if strcmpi(s(1).sensor_type,'Quasar Sensor') == 1
        Figure_sensor = 'Quasar';
    else
        Figure_sensor = 'VLF';
    end;
else
    Figure_sensor = 'ULF';
end

% Sampling interval is 0.01 ms
for i = 1:1:length(s)
    
    Event_ID = i;
    
    Figure_station = s(1,Event_ID).station_name;
    
    % Read the calibrated Bphi data
    B_Phi = s(1,Event_ID).sig_dn_Bphi_cal(1,1:end);
    B_r = s(1,Event_ID).sig_dn_Br_cal(1,1:end);
    
    Time_b = s(1,Event_ID).data_start_time;
    %Time_e = s(1,Event_ID).event_time
    
    Time_IT = Time_b(4)*3600 + Time_b(5)*60 + Time_b(6);
    
    Time_UT = Time_IT + (1/SF)*(linspace(1,length(B_Phi),length(B_Phi)) - 1);
    if s(i).event_time(3) <= 2.000 && SF < 10000
        Time_UT_SS = Time_UT - Time_b(4)*3600 - Time_b(5)*60 - 60;
        Time_b(5) = Time_b(5) + 1;
    else
        Time_UT_SS = Time_UT - Time_b(4)*3600 - Time_b(5)*60;
    end;
    
    distance = s(1,Event_ID).lightning_dist;
    eventtime = s(1,Event_ID).event_time;
    
    %         eventtime(3) = 17.409743;
    %         distance = 634;
    
    delaytime = distance/300000;
    Expected_VLF_Time = delaytime + eventtime(3) - sqrt((550)^2 + 000^2)*0/300000;
   
     Time_Event = eventtime(3);
     Time_Arrival = Expected_VLF_Time;
     
    if SF <= 10000
        Time_UT_SS = Time_UT_SS - 0.0014;
    end;
    
   
    absTime = Time_b(4)*3600+Time_b(5)*60;
    TGFd.t = Time_UT_SS + absTime;
    TGFd.B_phi = B_Phi;
    TGFd.B_r   = B_r;
    TGFd.eventTime = Time_Event;
    TGFd.arivlTime = Time_Arrival;
    TGFd.station = Figure_station;
    TGFd.sensor = Figure_sensor;

    
end;