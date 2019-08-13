function FA_FFT_01(fn,hfn,t1,t2,tshift,settings,sn)

% fn = 'C:/data/2011-08-05 -- 2011-08-16/data//2011/08/14/EDW_20110814_224000.ch3';
% hfn = 'C:/data/2011-08-05 -- 2011-08-16/data//2011/08/14/EDW_20110814_224000.h3';
% t1 =    8.1684e+04;
% t2 =   8.1685e+04;
% tshift = -1;
% settings = open('sensor_setting.mat');
% sn = 21;

% File name for SA data
%fn='Ols_20091207_212101.lp2';
 
% To get optimum results find FFT in full range
% give very small to t1 and very large to t2, then program will
% automatically picks all the points in that file.
 
%t1=0;   %8835;
%t2=500000000; % 8835.5;
 
% Importing data in the file
[tn,vn]=FA_Extract1(fn,hfn ,t1,t2,tshift,settings,sn);



% Number of data points imported
L=length(vn);

if L == 0
    msgbox('No data found in the given time range','FFT failed!')
    return
end
 
% Data frequency 
for i = 1:500
    Fs=1/(tn(i+1)-tn(i));

    if ~isnan(Fs)
        break
    end
end


% Add a high pass filter
[tn, vn] = ch_high_pass(tn,vn,5000);

% Number of discrete Fourier transform points
NFFT = 2^nextpow2(L); % Next power of 2 from length of data
 
% finding the Fourier transform
Y = fft(vn,NFFT)/L; % Y is complex
% frequency
f = Fs/2*linspace(0,1,NFFT/2+1);
 
YY=2*abs(Y(1:NFFT/2+1));
phase=atan(imag(Y(1:NFFT/2+1))./real(Y(1:NFFT/2+1)));
 
% Plot amplitude spectrum.
figure
plot(f,YY) 
title('FA fourier transform')
xlabel('Frequency (Hz)')
ylabel('Normalized apmlitude')


% Find peaks and report those
%try
    hold all
    
    [peakLoc,peakMag]=peakfinder(YY);
    
    freq=f(peakLoc);      
    
    
    % Identifying peaks
    str=sprintf('FFT Peak freq info:\n');
    str2=sprintf('     Freq               Norm Amp           Phase\n');
    str=[str str2];
    for i=1:10
        [m(i),index(i)]=max(YY);
        str2=sprintf('     %6.1f kHz          %3.2e         %1.6f Rad\n',f(index(i))/1000,m(i),phase(index(i)));
        str=[str str2];
        
        YY(1:index(i)+300)=0;
        plot(f(index(i)),m(i),'ro')
        
        if i == 1
            text(f(index(i)),m(i),[num2str(f(index(i))/1000) 'kHz'])
        end
            
    end
    
    hold off
    
    
    datacursormode on;
    
    msgbox(str,'FFT Peaks')
    disp(str)
% catch
%     disp('Finding peaks of FFT graph was not sucusess!!')
% end



 


