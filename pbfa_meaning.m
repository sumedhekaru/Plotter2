function pbfa_meaning

% After finding peaks, we just find pbfa points using pbfa finder. However,
% we did not have any idea what this position mean. This program is intend
% to answer that problem. 
%   1) Watson and Marshall NBP model will use to generate data at each
%   sensor locations.
%   2) Then highpass filter will apply to those data and then take the
%   hilbert transformation.
%   3) Using the hilbert transformation, the peaks can be determined.
%   4) This peaks information will be use to calculate PBFA position.
%   5) This position will be compared with the original model position.

% Arguments for pulse_simulator.m
%close all
 clc
                 x0 = 0;                    % Original position of the pulse
                 y0 = 0;                    % Original position of the pulse
             a.maxA = 750;                 % Max Current
               a.t1 = (32/6)*1e-6;          % Rise time
               a.t2 = 32e-6;                % Total time
                a.v = 6e7;                  % Pulse speed
               a.H1 = 6400;                 % Initial height
               a.H2 = 7030;                 % Final height
                a.x = 5000;                 % Sensor position
                a.y = 5000;                 % Sensor position
               a.T1 = 0.0e-005;             % Olot time range begining
               a.T2 = 8.5e-004;             % Plot time range end
            a.lamda = 100;                  % Lamda = (H2-H1)/ln560 See Watason & Marshall
           a.t_step = 1.0e-006;             % Time Resolution
                a.N = 1;                    % Bouncing current model (if N >= 2)
        a.t_reflect = 0;                    % Reflection time (When N>=2)
    a.coeff_reflect = -0.500000000000000;   % Reflection coefficient (When N > =2)
          a.t_shift = 0;                    % Manual time shift of curves
         a.is_saved = 0;                    % Related to GUI
          a.save_fn = '';                   % Related to GUI
        a.file_name = 'C:\Users\sumedhe\Desktop\Plotter\pulse_parameter_test.mat';
         a.plots_on = [0 0 0 0];            % Which plots should be on 
       a.addi_plots = [0 0 0 0];

       sen_set = open('sensor_setting.mat');
fg=figure; 
hold all
lg = {};


% Some variables to store things
handles.sen_num = [];
handles.t       = [];
handles.y       = [];
handles.real    = [];
handles.imag    = [];
handles.mag     = [];
handles.pn = 1;

sen_set.x = [1000 2000 5000 10000 20000 40000 80000 100000 150000 200000];
sen_set.sen_IDs = {'1000' '2000' '5000' '10000' '20000' '40000' '80000' '100000' '150000' '200000'};
sen_set.y = sen_set.x-sen_set.x;

for i = 1:10
    a.x = sen_set.x(i)- x0;
    a.y = sen_set.y(i)- y0;
    
    [t, E_tot] = pulse_simulator(a);
    % return
    % filter out low frequencies
    Fs = 1/(t(2)-t(1));
    [z,p] = butter(5,10000/(Fs/2),'high'); % Create a high-pass butterworth filter;
    % [z,p,k] = butter(n,Wn) designs an designs an order n lowpass digital Butterworth filter with normalized
    % cutoff frequency Wn. It returns the zeros and poles in length n column
    % vectors z and p, and the gain in the scalar k
    
    smoothE = filtfilt(z,p,E_tot);    % filter the data.
    %smoothy(1,10)
    % Hilbert transform of voltages
    y_hil = hilbert(smoothE);
    E_real = real(y_hil);
    E_imag = imag(y_hil);
    E_mag = sqrt(imag(y_hil).^2+real(y_hil).^2);
    
    clear peaks_i
    %[ peaks_i(:,1),  peaks_i(:,2)] = max(E_mag);
    %peaks_i = sortrows(peaks_i,-2);
    
    [peakMag peakLoc] = max(E_mag);
    
    %peakLoc = peaks_i(1,1);
    % peakMag = peaks_i(handles.pn,2);
    
    % Store Peak info for future use
    handles.sen_num = [handles.sen_num i];
    handles.t       = [handles.t t(peakLoc)];
    handles.y       = [handles.y E_tot(peakLoc)];
    handles.real    = [handles.real E_real(peakLoc)];
    handles.imag    = [handles.imag E_imag(peakLoc)];
    handles.mag     = [handles.mag E_mag(peakLoc)];
    
    % Calculating peak Current
%     temp1 = abs(E_tot(peakLoc) - E_tot(peakLoc - 1));
%     temp2 = abs(E_tot(peakLoc) - E_tot(peakLoc + 1));
%     
%     if temp1 > temp2
%         temp1*sqrt(a.x^2+a.y^2+((a.H1+a.H2)/2)^2)/100000
%     else
%         temp1*sqrt(a.x^2+a.y^2+((a.H1+a.H2)/2)^2)/100000
%     end
    
    E_mag(peakLoc)*sqrt(a.x^2+a.y^2+((a.H1+a.H2)/2)^2)/100000;
    
    %subplot(2,2,1); hold all; box on;
    figure(fg)
    xlabel('Time (s)'); ylabel('E-Field (V/m)');
    plot(t,E_tot)
    lg = [lg sen_set.sen_IDs{i}];
    
    legend(lg)
    
%     subplot(2,2,2); hold all; box on;
%     xlabel('Time (s)'); ylabel('E-Real Part (V/m)');
%     plot(t,E_real)
%        
%     subplot(2,2,3); hold all; box on;
%     xlabel('Time (s)'); ylabel('E-Imaginary Part (V/m)');
%     plot(t,E_imag)   
%     
%     subplot(2,2,4); hold all; box on;
%     xlabel('Time (s)'); ylabel('E-Magnitude (V/m)');
%     plot(t,E_mag);
    
    
    

end

return

% Plot peaks and Display times
fprintf('\n***********************************\n')
for i = 1:10
    fprintf('\t\t%s:%.9fs\n',sen_set.sen_IDs{i},handles.t(i))
    
    subplot(2,2,1); plot(handles.t(i), handles.y(i),'pr','markerFaceColor','r'); 
    subplot(2,2,2); plot(handles.t(i), handles.real(i),'pr','markerFaceColor','r');
    subplot(2,2,3); plot(handles.t(i), handles.imag(i),'pr','markerFaceColor','r');
    subplot(2,2,4); plot(handles.t(i), handles.mag(i),'pr','markerFaceColor','r');
end
fprintf('\n***********************************\n')

handles.t_accu = a.t_step;
cal_pbfa(handles)
 

function cal_pbfa(handles)

%%%  Input Arguments for pbfa_finder %%%
% arg.inner = [ 1 2 3 6 11];      % Inner sensors
% arg.outer = [ 5 7 8 9 10];      % Outer sensors
% 
% arg.t_in  = getT.t_in;
% arg.t_out = getT.t_out;
% 
% arg.method = 1;

% find inner and outer sensors

% PBFA error calculations
error_cal = 1;
N         = 1000; % Number of itterations

    arg.inner =[];
    arg.outer =[];
    arg.t_in  =[];
    arg.t_out =[];
    arg.method=1;
    
    for i=1:length(handles.t)
        % Let's choose sensors that have time
        if ~isnan(handles.t(i))
            if handles.sen_num(i) < 6
                % considering Ksites, BCC, WSB are inner
                arg.inner = [arg.inner handles.sen_num(i)];
                arg.t_in  = [arg.t_in handles.t(i)];
            else
                arg.outer = [arg.outer handles.sen_num(i)];
                arg.t_out = [arg.t_out handles.t(i)];
            end
        end
    end
    
    if length([arg.t_in arg.t_out]) < 5
        errordlg('Atleast 5 sensors need to calculate a PBFA point','PBFA error')
    else
       [x,y,z,t]=pbfa_finder(arg);
       
       if error_cal
           [dx,dy,dz,dt]=pbfa_error(x,y,z,[arg.outer arg.inner],...
               handles.t_accu,N,1,2);
           
           str = sprintf('x = (%.1f +/- %0.1f)m\ny = (%.1f +/- %0.1f)m \nz = (%.1f +/- %0.1f)m \nt = (%6.9f +/- %0.9f)s\n',...
               x,dx,y,dy,z,dz,t,dt);
          
       else          
           
           str = sprintf('x = %.1fm\ny = %.1fm \nz = %.1fm \nt = %6.9fs\n?',...
               x,y,z,t);
       end
        %button = questdlg(str,'PBFA Finder','Yes','No','No') ;
        fprintf('%s',str)
        fprintf('\n**********************************\n')
        
        handles.pbfa.x = x;
        handles.pbfa.y = y;
        handles.pbfa.z = real(z);
        handles.pbfa.t = t;
              
            
    end




