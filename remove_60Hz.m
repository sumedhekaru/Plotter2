function vn = remove_60Hz(t,v)


% % Find the sampling frequency
% for i = 1:10
%     Fs = 1/(t(i+1) - t(i));
%     
%     if ~isnan(Fs)
%         break
%     end
% end
%     
% 
% d = designfilt('bandstopiir','FilterOrder',10, ...
%                'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
%                'DesignMethod','butter','SampleRate',Fs);
%            
% vn = filtfilt(d,v);
% %vn = v;

% Changing file names FFI to FF2

% sen_set = open('sensor_setting.mat');
% 
% dir = sprintf('%s/2011/08/09/',sen_set.base_dir);
% 
% files = ls([dir 'FFI*']);
% 
% for i = 1:length(files)
%     file = [dir files(i,:)];
%     toFile = file;
%     toFile(end-22:end-20) = 'FF2'
%     movefile(file,toFile)
% end

%FLT
A = 0.02580;
ph = 227.60;
f = 59.999780;

%STC
A = 0.03280;
ph = 72;
f = 60;

sny = A*sin(2*pi*f*t + ph/180*pi);

vn = v - sny';