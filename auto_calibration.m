function auto_calibration
% Auto calibrate of lp2 data

% Load plotter2 data
try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

settings = h.sen_set;
g = h.g;
tshift=settings.t_shift;


% Search n number of CGLSS points within delT time
delT = 1.0;
n = 10; 

% Plotting time range
ptr = 2;

% Sensor number of calibrating
sn = 1;
sid = settings.sen_IDs{sn};

% Guess minimum and maximum gain
maxG = 50000;
minG = 0000;

% Directory to save data
savedir      = 'C:\Users\daqop\Desktop\Sumedhe\Calibration\Auto\';

savedir      = sprintf('%s/%s%s%s/%s/lp2_auto/',savedir,g.YYYY{:},g.MM{:},g.DD{:},sid);

if ~exist(savedir,'dir')
    mkdir(savedir)
end
    
% File name to save data
sv_cal_data = [savedir 'lp2_Calibration_data.txt'];

if exist(sv_cal_data, 'file')
    fID = fopen([sv_cal_data],'a+');
else
    % If not exisit open and write the header
    fID = fopen([sv_cal_data],'a+');
    fprintf(fID,'Plot_name\tGain\tSA_range\tFM_range\n');
end




%% Load CGLSS data
cglss_fn = sprintf('%s/cglss/%s/%s/KSCCGLSS%s%s%s.dat',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.YYYY{:},g.MM{:},g.DD{:});

x0 = 0;
y0 = 0;
z0 = 0;


%CG = CGLSS_extract(cglss_fn,76000,86400,500000,x0,y0,z0);
CG = CGLSS_extract(cglss_fn,68100,86400,500000,x0,y0,z0);
%CG = CGLSS_extract(cglss_fn,0,10800,500000,x0,y0,z0);

%% Find single CGLSS points
t=CG(:,1);                       % time
x=CG(:,2);                         % East
y=CG(:,3);                         % North
z=CG(:,12);                         % Altitude

clear CG

%length(t)
k = 1;
T = [];

if isempty(t)
    fprintf('No CGLSS in the given time range\n')
else
    for i=1:length(t)
        t1=t(i)-delT;
        t2=t(i)+delT;
        %lol=size(CG(:,10),1)- nnz(CG(:,10)>t1)+1;      % Index of the lower matrix element
        lol = sum(t < t1)+1;
        ul = sum(t <= t2);                           % Index of the upper matrix element
        slice=t(lol:ul);                             % CG data for given time interval
        if length(slice) < n + 1
            T(k)=i;
            k=k+1;
        end
    end
end

% Turn on the correct graphs
g.lpgraphs = zeros(1,60);
g.lpgraphs((sn-1)*3+2) = 1;

if isempty(T)
    fprintf('\nNo (%d strock/s) CGLSS were found (Time threshold = %0.3fs)\n',n,delT)
else
    wbh = waitbar(0,['Calibrating ' sid '...'] ,'Name','Auto Calibrator');
    for i = 1:n:k-1
        
        val = i/(k-1);
        wbmsg = sprintf('Calibrationg %s... (%0.2f%%)',sid,val*100);
        waitbar(val,wbh,wbmsg)
        
        ind = T(i);
        tcg = t(ind);
        
        % Propagation time
        dt = sqrt((settings.x(sn)-x(ind))^2 + ... 
                          (settings.y(sn)-y(ind))^2)/3e8;
        
        
        g.t1 = tcg + dt - ptr/2;
        g.t2 = g.t1  + ptr;
        
        % Generate ch file names according to the current CGLSS poing
        hhmmss = sec2hhmmss(tcg);
        g.hh = str2double(hhmmss(1:2));
        g.mm = floor(str2double(hhmmss(4:5))/5)*5;
        
                
        [lp_fn fm_fn] = generate_lp_fn(settings,g);
        
        % Store one data set for plot yy        
        for j=1:60
            % for loop used becuase there may be 15 lp plots
            if strcmp(lp_fn{j},'')==0
                [tn,vn]=SA_Extract1(lp_fn{j},g.t1,g.t2,tshift(j));
                % Summing to 50Hz
                ni = 200;
                
                len = floor(length(tn)/ni);
                
                temp1 = nan(ni,len);
                temp2 = temp1;
                
                for indn = 1:ni
                    tempt= downsample(tn,ni,indn-1);
                    tempy= downsample(vn,ni,indn-1);
                    temp1(indn,:)=tempt(1:len);
                    temp2(indn,:)=tempy(1:len);
                    clear tempt tempy                    
                end
                
                tn = mean(temp1);
                vn = mean(temp2);
                %tn = downsample(tn,200);
                %vn = downsample(vn,200);
            end
        end
        
        
        % Load FM data
        [tf,vf]=fieldMillExtract3(fm_fn,g.t1,g.t2,-1);
        vf = -vf; % Field mill has a different sign convention
        vf = vf-min(vf);
        
        % find gain
        dvf  = range(vf);
        dvn = range(vn);
        gain = dvf/dvn;
        
        if ~isempty(gain)
            % find shift
            maxf = max(vf);
            maxn = max(vn*gain);
            shift = maxn - maxf;
            
            if dvn > 0.02 && dvf > 20 && gain < maxG && gain > minG
                
                % find the time shift
                t0 = min(tn);
                temp = min(tf);
                
                if temp < t0
                    t0 = temp;
                end
                
                % Plot data
                fg = figure;
                hold all
                plot(tn-t0,vn*gain-shift)
                plot(tf-t0,vf)
                
                lg{1} = [ sid ':lp2'];
                lg{2} = 'FieldMill';
                
                str = sprintf('%s-%s-%s\nGain = %0.1f     SA Range = %0.3f      FM Range = %0.1f',...
                    g.YYYY{:},g.MM{:},g.DD{:},gain,dvn,dvf);
                title(str)
                legend(lg)
                xlim([0 ptr])
                xlabel(sprintf('Time after %0.6fs',g.t1))
                ylabel('\Delta E (v/m)')
                box on
                
                secnd = floor(g.t1);
                milsec = round((g.t1-secnd)*1000);
                figfn = sprintf('%s_Auto_%s%s%s_%5.5i_%3.3i.fig',...
                    sid,g.YYYY{:},g.MM{:},g.DD{:},secnd,milsec);
                
                saveas(fg,[savedir figfn])
                delete(fg)
                
                % Write data in to the file
                fprintf(fID,'%s\t%0.1f\t%0.3f\t%0.1f\n',...
                    figfn,gain,dvn,dvf);
                
            end
        end
    end
end
try delete(wbh); end
try fclose(fID); end

function [lp_fn fm_fn] = generate_lp_fn(settings,g)
% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

% checking witch graphs are on
lp_on=g.lpgraphs;
lp_fn = cell(1,60);

for i=1:60
    if lp_on(i)==1
          ext=mod(i,3);
        if ext==0
            ext=3;
        end
        
        % Finding the stattion ID
        sid=settings.sen_IDs{ceil(i/3)};
        
        
        % finding the file name
        filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.lp%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        % Check whether the file is exists
        if exist(filename,'file')~=0
            % Check whether the file is exists
            lp_fn{i}=filename;
        else
            % If file is not exist don't store the file name
            lp_fn{i}='';
        end
    else
        lp_fn{i}='';
    end
end

bfolder=sprintf('%s/FM/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

if g.mm < 30
    mm = 0;
else
    mm = 30;
end

fm_fn = sprintf('%sKSC%s%s%s%s_%2.2i%2.2i.txt',...
       bfolder,sid(2:3),g.YYYY{:},g.MM{:},g.DD{:},g.hh,mm);
