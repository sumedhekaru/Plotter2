function data = load_bn_dataV03(fileName,t1,t2)

fileName = 'H:\data\Gage\2015\06\08\17\EE4_20150608_170631.bn4';
% 
% findOut t1 t2 from file name
hh = str2double( fileName(end-9:end-8));
mm = str2double(fileName(end-7:end-6));
ss = str2double(fileName(end-5:end-4));

tRef = hh*3600+mm*60+ss +1;

t1 =  tRef - 0.01;
t2 =  tRef + 0.01;
%t1 = 62650.506;
%t2 = 62650.509;

% skip data
nskip = 4;

fID = fopen(fileName,'r');
% Number of heder variables
%hL = fread(fID,1,'int64');

% Go to the end of the file and revind 8 byte to read the number of header
% lines
fseek(fID,-8,'eof');  
hL = fread(fID,1,'int64');

% Seek to the header and read it
fseek(fID,-8*(hL+1),'eof');
hdr = fread(fID,hL,'int64');

% Decode the header
hd.version          = hdr(1);  
hd.chan.Channel     = hdr(2);
hd.chan.Coupling    = hdr(3);
hd.chan.InputRange  = hdr(4);
hd.chan.Impedance   = hdr(5);
hd.chan.DcOffset    = hdr(6);
hd.chan.Filter      = hdr(7);
hd.chan.sum         = hdr(8);
hd.chan.delta       = [hdr(9), hdr(10)]/1000;
hd.chan.pre         = hdr(11)/1000;
hd.chan.post        = hdr(12)/1000;
hd.st_time          = hdr(13:17);
hd.st_time(6)       = hdr(18)/10^9;
hd.freq             = hdr(19);
hd.PPS_locked       = hdr(20);
hd.acq.SampleRate   = hdr(21);
hd.acq.SampleBits   = hdr(22);
hd.acq.SampleResolution = hdr(23);
hd.acq.SampleSize   = hdr(24);
hd.acq.SegmentCount = hdr(25);
hd.acq.Depth        = hdr(26);
hd.acq.SegmentSize  = hdr(27);
hd.acq.SampleOffset = hdr(28);
hd.acq.TimeStampConfig = hdr(29);
hd.nOfTriggers      = hdr(30);
hd.chan.presN       = hdr(31:31+hd.nOfTriggers-1);
hd.chan.postsN      = hdr(31+hd.nOfTriggers:31+2*hd.nOfTriggers -1);
hd.chan.presT       = hdr(31+2*hd.nOfTriggers:31+3*hd.nOfTriggers -1);
hd.chan.sts         = hdr(31+3*hd.nOfTriggers:end);

% Let's read data (according to user requirements, loading all data is just
% a wast of resources.

figure
hold all;

% Lengths of each trigger
Ls = diff(hd.chan.sts);

hd.freq = hd.acq.SampleRate;

% Let's go trough all the triggers 
for i = 1:hd.nOfTriggers
    
    % Real begining and end time of this
    beg_t = hd.chan.presT(i)/1e9; % devide by 10^9 to convert nano secods to second
    end_t = beg_t + (Ls(i)-1)/hd.freq;
    fprintf('%0.9f\n',beg_t)
    
    % DEBUG - Just load whole triggers to see whether this is working.
    % If you choose to do this, make sure the pre-post triggers are not
    % much long. Long triggers will freez the computer when plotting.
    %fseek(fID,(hd.chan.sts(i)-1)*2,'bof');
    %data.v = (((hd.acq.SampleOffset - double(fread(fID,Ls(i),'int16'))) / hd.acq.SampleResolution) * ...
    %        (hd.chan.InputRange / 2000)) + (hd.chan.DcOffset / 1000);
    
    %data.t = beg_t:1/hd.freq:beg_t + (Ls(i)-1)/hd.freq;        
    %plot(data.t,data.v)
                
    % Find out starting byte for seeking
    if t1 <= beg_t
       nd1 = 0;                     % Number of data points
       stbt = (hd.chan.sts(i)-1)*2; % Number of byte to seek
       stT = beg_t;                 % Starting time
    else
       %Startubg byte for seeking
       nd1 = round((t1 - beg_t)*hd.freq); % starting number of data points
       stbt = (hd.chan.sts(i) + nd1 - 2)*2; % Number of byte to seek
       stT = beg_t + (nd1-1)/hd.freq;
    end
    
    %nd1
    
    % Find out how many points we should read
    if t2 > end_t  
       readCount = Ls(i) - nd1;       % Number of data pointds should be readed
    else
       %disp('working')
       fprintf('%0.9f\n',t2)
       fprintf('%0.9f\n',end_t)
       nd2 = round((end_t - t2)*hd.freq);
       readCount = Ls(i) - nd2 - nd1;         
    end
    
%     Ls(i)
%     nd2
     enT = stT + (readCount-1)/hd.freq;
    
     if nskip > 1
        readCount = readCount/nskip;
        data.t = stT:nskip/hd.freq:stT+nskip*(readCount-1)/hd.freq; % Time vector
     else
         data.t = stT:1/hd.freq:enT;
     end
     
     %nskip
     %t_length = length(data.t)
     %readCount
%     Ls(i)
   
 
    if readCount > 1
        
        %fseek(fID,1,'bof');
        %data.v = fread(fID,readCount,'int16');
        
        % Read data
        fseek(fID,stbt+nskip*126,'bof');
        data.v = (((hd.acq.SampleOffset - double(fread(fID,readCount,'int16',nskip))) / hd.acq.SampleResolution) * ...
            (hd.chan.InputRange / 2000)) + (hd.chan.DcOffset / 1000);
        %dvdt = diff(data.v);
        %length(data.v)
        
        plot(data.t,data.v)
        %xlim([35665.9999999 35666.0000001])
        %xlim([t1 t2])
        %plot(data.t(1:end-1),dvdt)
        %plot(data.t(1:end-2),diff(dvdt))
    end
end




