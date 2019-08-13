function  [yF, yH] = hill_tra(t,yn,f)
% This function will find the hilbert trabsform of (t,yn) data. f is the
% high pass filter frequency. This will return hilbert transfor yH and 
% filtered data yF 

            yn(isnan(yn))=0;
            yH  =[];
            yF = [];
            Fs = NaN;
            
            % filter out low frequencies
            i = 1;
            L = length(t);
            
            if L == 0
                return
            end
            
            while isnan(Fs) && i < L
                Fs = 1/(t(i+1)-t(i));
                i = i + 1;
            end      
           
            
            [z,p] = butter(5,f/(Fs/2),'high'); % Create a High-pass butterworth filter;
          
            try
                yF = filtfilt(z,p,yn);    % filter the data.
            catch
                yF = [];               
            end
       
            
            if ~isempty(yF)
                % Hilbert transform of voltages
                yH = hilbert(yF);
            end
            
                