function SA_FFT_01(fn,t1,t2,tshift)

% File name for SA data
%fn='Ols_20091207_212101.lp2';
 
% To get optimum results find FFT in full range
% give very small to t1 and very large to t2, then program will
% automatically picks all the points in that file.
 
%t1=0;   %8835;
%t2=500000000; % 8835.5;
 
% Importing data in the file
[tn,vn]=SA_Extract1(fn,t1,t2,tshift);
 
% Number of data points imported
L=length(vn);
 
% Data frequency 
Fs=L/(max(tn)-min(tn)); % ~ 10kHz
 
% Number of discrete Fourier transform points
NFFT = 2^nextpow2(L); % Next power of 2 from length of data
 
% finding the Fourier transform
Y = fft(vn,NFFT)/L; % Y is complex
% frequency
f = Fs/2*linspace(0,1,NFFT/2+1);
 
YY=2*abs(Y(1:NFFT/2+1)); % 
phase=atan(imag(Y(1:NFFT/2+1))./real(Y(1:NFFT/2+1)));
 
% Plot amplitude spectrum.
figure
plot(f,YY) 
title('SA fourier transform')
xlabel('Frequency (Hz)')
ylabel('Normalized apmlitude')
ylim([0 .025]);
 
try
    hold all
    
    [peakLoc,peakMag]=peakfinder(YY(3000:end));
    
    freq=f(peakLoc);
    
    % Identifying peaks
    str=sprintf('FFT Peak freq info:\n');
    str2=sprintf('     Freq               Norm Amp           Phase\n');
    str=[str str2];
    for i=1:10
        [m(i),index(i)]=max(YY);
        str2=sprintf('     %6.1f Hz          %.6f         %1.6f Rad\n',f(index(i)),m(i),phase(index(i)));
        str=[str str2];
        
        YY(1:index(i)+300)=0;
        plot(f(index(i)),m(i),'ro')
    end
    
    hold off
    
    
    datacursormode on;
    
    msgbox(str,'FFT Peaks')
    disp(str)
catch
    errordlg('Finding peaks of FFT graph was not sucusess!!','FFT plot error','modal')
end



 


