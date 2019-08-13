function TGF_trigger_finder
% This function will use to find number of sensors triggered for TGF list
% sent by Gopeng Lu. This function will find the average LDAR2 X,Y,R
% locations and plot the +/- 100ms worth of data arround TGF and save.
% Figures will be saved as yyyymmdd--hhmmss-sss-XY.fig
%                          yyyymmdd--hhmmss-sss-ch.fig

%% User inputs
% Make sure to select date time from the plotter2 and then plot something
% to test whether it's in correct base folder and correct manual time
% shifts were loaded.

% Base folder to save stuffs
bf = 'C:\Users\Sumedhe\Desktop\TGF\20110814-2\';

% Counter
count = 1;

% Times from Gaopeng
times = {'00:27:17.635'
'00:29:2.577'
'00:30:28.819'
'00:30:28.922'
'00:30:28.952'
'00:32:52.650'
'19:34:38.216'
'19:35:38.151'
'19:35:49.481'
'19:35:58.419'
'19:36:9.041'
'19:36:23.227'
'19:36:23.847'
'19:36:27.089'
'19:36:28.919'
'19:37:48.854'
'19:38:4.560'
'19:38:33.807'
'19:38:35.673'
'19:38:45.699'
'19:39:32.736'
'19:40:12.024'
'19:40:13.072'
'19:40:41.841'
'19:40:44.611'
'19:40:45.357'
'19:40:46.455'
'19:40:51.049'
'19:41:20.521'
'19:41:20.880'
'19:41:22.857'
'19:41:28.097'
'19:41:48.114'
'19:41:58.534'
'19:42:13.917'
'19:42:27.147'
'19:42:42.848'
'19:43:46.631'
'19:44:36.684'
'19:44:37.167'
'19:45:19.101'
'19:46:19.262'
'19:46:24.612'
'19:47:3.293'
'19:47:45.028'
'19:49:5.364'
'19:49:20.260'
'19:49:55.661'
'19:49:57.293'
'19:50:25.674'
'19:50:55.575'
'19:52:3.533'
'19:53:1.238'
'19:54:44.444'
'19:58:20.081'
'20:02:34.843'
'20:04:35.121'
'20:06:39.213'
'20:07:55.235'
'20:09:9.583'
'20:09:48.019'
'20:15:56.588'
'20:18:44.492'
'20:22:4.570'
'20:22:58.689'
'20:28:40.886'
'20:31:39.740'
'20:34:6.248'
'20:34:57.183'
'20:47:23.545'
'20:50:20.693'
'20:55:0.873'
'20:59:6.583'
'21:02:1.769'
'21:17:32.708'
'21:17:41.091'
'21:19:20.783'
'21:20:6.446'
'21:20:49.676'
'21:20:51.055'
'21:20:51.466'
'21:23:9.797'
'21:35:30.266'
'21:36:55.923'
'21:37:12.895'
'21:37:29.820'
'21:39:38.083'
'21:42:1.078'
'21:44:19.656'
'21:45:9.888'
'21:46:2.450'
'21:47:15.923'
'21:47:26.092'
'21:47:29.935'
'21:48:2.919'
'21:48:14.779'
'21:48:20.775'
'21:48:23.192'
'21:48:48.023'
'21:49:45.271'
'21:49:53.801'
'21:50:29.084'
'21:50:32.954'
'21:50:36.296'
'21:50:55.116'
'21:51:1.188'
'21:51:14.403'
'21:51:21.034'
'21:52:18.463'
'21:52:44.168'
'21:53:9.614'
'21:53:29.420'
'21:53:31.916'
'21:53:45.874'
'21:53:54.386'
'21:54:31.169'
'21:56:0.139'
'21:57:25.730'
'21:58:8.536'
'21:58:12.062'
'22:00:42.277'
'22:01:13.908'
'22:01:53.014'
'22:02:11.079'
'22:03:55.339'
'22:03:59.096'
'22:04:43.498'
'22:05:11.454'
'22:05:37.417'
'22:05:54.466'
'22:06:31.226'
'22:06:40.610'
'22:06:44.295'
'22:06:56.796'
'22:07:4.105'
'22:07:31.062'
'22:07:56.238'
'22:08:3.146'
'22:08:13.306'
'22:08:45.356'
'22:08:54.758'
'22:09:13.117'
'22:09:44.696'
'22:10:27.922'
'22:10:47.455'
'22:12:34.798'
'22:12:35.123'
'22:13:14.746'
'22:13:20.399'
'22:13:57.148'
'22:14:42.116'
'22:14:42.661'
'22:14:59.659'
'22:16:9.360'
'22:16:31.843'
'22:16:48.254'
'22:17:32.770'
'22:17:41.053'
'22:17:46.149'
'22:18:43.956'
'22:19:57.705'
'22:20:28.724'
'22:20:43.196'
'22:20:57.903'
'22:24:59.472'
'22:25:9.154'
'22:25:18.725'
'22:25:39.146'
'22:27:10.064'
'22:28:28.840'
'22:28:44.039'
'22:28:44.945'
'22:28:51.739'
'22:29:10.572'
'22:29:16.169'
'22:29:28.121'
'22:29:30.916'
'22:29:44.508'
'22:29:45.459'
'22:29:46.750'
'22:29:52.323'
'22:30:14.877'
'22:30:20.206'
'22:30:20.863'
'22:30:43.242'
'22:30:46.653'
'22:30:58.626'
'22:30:58.635'
'22:31:12.778'
'22:31:18.583'
'22:31:19.727'
'22:31:23.577'
'22:31:37.580'
'22:31:49.681'
'22:32:13.630'
'22:32:42.983'
'22:33:11.559'
'22:33:15.387'
'22:33:33.423'
'22:33:48.163'
'22:34:8.871'
'22:34:12.847'
'22:34:43.216'
'22:34:55.749'
'22:35:0.860'
'22:35:18.672'
'22:35:26.035'
'22:35:36.658'
'22:35:57.976'
'22:36:28.844'
'22:36:32.543'
'22:36:34.423'
'22:36:42.394'
'22:36:47.156'
'22:37:15.706'
'22:37:16.741'
'22:37:28.806'
'22:37:33.320'
'22:38:1.300'
'22:38:50.164'
'22:38:54.072'
'22:39:22.085'
'22:39:38.534'
'22:39:54.922'
'22:40:44.331'
'22:41:4.320'
'22:41:24.357'
'22:41:31.175'
'22:41:51.133'
'22:42:15.933'
'22:42:42.146'
'22:42:50.786'
'22:43:12.009'
'22:43:30.013'
'22:43:47.361'
'22:44:49.103'
'22:44:52.559'
'22:45:10.197'
'22:45:21.807'
'22:45:51.831'
'22:46:10.766'
'22:46:32.349'
'22:46:43.674'
'22:46:44.704'
'22:47:22.536'
'22:47:30.537'
'22:47:40.444'
'22:47:43.190'
'22:47:43.312'
'22:47:50.746'
'22:47:53.444'
'22:48:0.371'
'22:48:8.706'
'22:48:21.360'
'22:49:5.584'
'22:49:21.199'
'22:49:42.877'
'22:49:55.741'
'22:50:15.705'
'22:50:20.304'
'22:50:36.800'
'22:50:47.578'
'22:50:53.772'
'22:50:57.331'
'22:51:36.981'
'22:52:26.125'
'22:52:51.436'
'22:52:56.836'
'22:52:58.649'
'22:53:23.313'
'22:55:15.691'
'22:55:20.826'
'22:55:23.166'
'22:56:52.853'
'22:56:56.364'
'22:57:12.841'
'22:59:43.323'
'23:00:12.661'
'23:01:16.731'
'23:01:43.283'
'23:02:37.974'
'23:02:38.853'
'23:02:40.947'
'23:04:4.345'
'23:04:17.584'
'23:05:13.804'
'23:06:18.092'
'23:06:46.512'
'23:07:10.336'
'23:07:49.399'
'23:09:16.171'
'23:11:19.711'
'23:11:41.692'
'23:11:52.659'
'23:13:2.681'
'23:14:17.982'
'23:14:17.989'
'23:15:39.676'
'23:15:49.850'
'23:16:37.505'
'23:18:1.271'
'23:18:7.354'
'23:18:12.969'
'23:18:27.428'
'23:19:36.156'
'23:20:6.262'
'23:20:59.319'
'23:22:24.618'
'23:22:25.644'
'23:23:22.880'
'23:27:39.959'
'23:30:36.867'
'23:31:25.600'
'23:32:41.597'
'23:32:41.820'
'23:33:48.760'
'23:34:54.104'
'23:35:1.002'
'23:35:8.280'
'23:35:38.797'
'23:35:51.982'
'23:35:57.736'
'23:37:10.530'
'23:37:17.283'
'23:37:57.100'
'23:38:14.262'
'23:38:43.230'
'23:38:56.558'
'23:38:56.567'
'23:39:18.647'
'23:40:5.889'
'23:40:34.206'
'23:40:58.760'
'23:41:42.804'
'23:42:19.246'
'23:42:25.591'
'23:42:25.716'
'23:43:14.672'
'23:43:42.023'
'23:43:42.030'
'23:43:42.097'
'23:44:34.811'
'23:44:35.134'
'23:46:42.997'
'23:46:52.845'
'23:46:53.278'
'23:48:23.284'
'23:48:31.153'
'23:49:16.196'
'23:49:36.913'
'23:50:6.748'
'23:50:37.295'
'23:50:59.678'
'23:54:18.690'
'23:54:46.991'
'23:55:11.465'
'23:55:37.928'
'23:56:14.455'
'23:56:30.876'
'23:56:30.942'
'23:57:13.707'
'23:57:27.489'
'23:58:13.775'
'23:58:17.542'


};


%% Prepare

% Get current plotter2 data
h=guidata(findall(0,'Tag','plotter2'));
g = h.g;

% Turn on all ch3 plots (avoid WSB and add Ch1 of FFI)
g.chgraphs = [0 0 1 0 0 1 0 0 1 0 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 ...
    1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

settings = h.sen_set;
tshift=settings.t_shift;

date = sprintf('%s%s%s',g.YYYY{:},g.MM{:},g.DD{:});

sv_fn = sprintf('%s%ss-TGF_trigger_info.txt',bf,date);

% Create file if not exist and add the header
% Write the header if not exit
if exist(sv_fn,'file')
    fID = fopen(sv_fn,'a+');
else
    fID = fopen(sv_fn,'a+');
    fprintf(fID,'count\tx\ty\tz\tr\ttriggers\n');
end


L = length(times);
wbh = waitbar(0,'Finding triggeres..','name','TGF trigger finder');


%% main loop
for i=1:L
   try 
        val = i/L;
        wbstr = sprintf('Finding triggers... (%0.1f%%)',val);
        waitbar(val,wbh,wbstr)
    catch
        return
   end
   

   t = hhmmss2sec(times{i},0);
   g.t1 = t - 0.05;
   g.t2 = t + 0.05;
   
   %% Loading and plotting ch3 data
   % Generate ch file names according to the current CGLSS poing
    hhmmss = sec2hhmmss(t);
    g.hh = str2double(hhmmss(1:2));
    g.mm = floor(str2double(hhmmss(4:5))/5)*5;
    
    [ch_fn h_fn] = generate_ch_fn(settings,g);
    triggered = 'none';
    
    wbh2 = waitbar(0,'Checking data... ');
    % Check whether sensors have triggered data
    
    peak_dE1 = 0;
    peak_dE_sen = '';
    
    for j=1:60;
        
        waitbar(j/60,wbh2,sprintf('Checking data... (%0.1f%%)',j/6*10))
        if ~isempty(ch_fn{j})
            [tch ych ch_freq_str] = FA_Extract1(ch_fn{j},h_fn{j},g.t1,g.t2,...
                tshift(j),settings,j);
            if ~isempty(tch)
                if strcmp(triggered,'none')
                    triggered = sprintf('%i',ceil(j/3));
                else
                    triggered = sprintf('%s,%i',triggered,ceil(j/3));
                end
                
                peak_dE_temp = range(ych);
                
                if peak_dE_temp > peak_dE1
                    peak_dE1 = peak_dE_temp;
                    peak_dE_sen = ceil(j/3);
                end
                    
                
            end
            
        end
    end

    count = count + 1;
    
     delete(wbh2)
%     plot_all4(g,0)
%     
%     timestr = times{i};
%     inx = find(timestr == ':');
%     hh = str2double(timestr(1:inx(1)-1));
%     mm = str2double(timestr(inx(1)+1:inx(2)-1));
%     ss_tmp = str2double(timestr(inx(2)+1:end));
%     ss  = floor(ss_tmp);
%     ss_dec = round((ss_tmp - ss)*1000);
%     
%     xy_fig_name = sprintf('%s%4.4i-%s--%2.2i%2.2i%2.2i-%3.3i-ch.png',bf,count,date,hh,mm,ss,ss_dec);
%     
%     lgh = findobj(gcf,'Type','axes','Tag','legend');
%     set(lgh,'color','none')
%     
%     %screen_size = get(0, 'ScreenSize');
%     %set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%     
%     saveas(gcf,xy_fig_name)
%     close(gcf)
    
    
    
    
    %% Loading and Plotting LDAR2 data    
    ldar_fn = generate_ldar_fn(settings,g);
    
    [CG,CAL,DLS]=ldarExtract2(ldar_fn,g.t1,g.t2,100000,0,0,0,0);
      
    ts=DLS(:,10);         % time
    xmin = nan;
    ymin = nan;
    zmin = nan;
    rmin = nan;
    
    if ~isempty(ts)
        % find closest ldar point to time "t"
        [minVal,minInd] = min(abs(ts - t));
        xmin = DLS(minInd,6)/1000;
        ymin = DLS(minInd,7)/1000;
        zmin = DLS(minInd,8)/1000;
        rmin = sqrt(xmin^2+ymin^2);
        
    end
%    
%     ldarColorTime1(ldar_fn,'','','',g.t1,g.t2,...
%             str2double(settings.ldar_r),0,0,0,1,[0,0,1,0,0],settings);
   
    %xy_fig_name = sprintf('%s%4.4i-%s--%2.2i%2.2i%2.2i-%3.3i-XY.png',bf,count,date,hh,mm,ss,ss_dec);

    %saveas(gcf,xy_fig_name)
    close(gcf)  
    
    % Write info to the file    
    fprintf(fID,'%i\t%0.3f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%s\t%0.3f\t%i\n',...
        count,t,xmin,ymin,zmin,rmin,triggered,peak_dE1,peak_dE_sen);
    
    
    
end
fclose(fID);
delete(wbh)


function [ch_fn h_fn] = generate_ch_fn(settings,g)
%% Generating ch file names
% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

% checking witch graphs are on
ch_on=g.chgraphs;
ch_fn = cell(1,60);
h_fn = cell(1,10);

for i=1:60
    if ch_on(i)==1
        % finding the file extention number
        ext=mod(i,3);
        if ext==0
            ext=3;
        end
        
        % Finding the stattion ID
        sid=settings.sen_IDs{ceil(i/3)};
        
        % finding the file name
        filename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.ch%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        
        hfilename=sprintf('%s%s_%s%s%s_%2.2i%2.2i%2.2i.h%1.1i', ...
            bfolder,sid,g.YYYY{:},g.MM{:},g.DD{:},g.hh,g.mm,g.ss,ext);
        
        
        % If file is not exist don't store the file name
        if exist(filename,'file')==0 || exist(hfilename,'file')==0
            ch_fn{i}='';
            h_fn{i}='';
        else
            ch_fn{i}=filename;
            h_fn{i}=hfilename;
        end
    else
        ch_fn{i}='';
        h_fn{i}='';
    end
end

function ldar_fn = generate_ldar_fn(sen_set,g)

       
        if g.mm < 30
            ext=0;
        else
            ext=30;
        end
        
        
        dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
            -datenum(str2double(g.YYYY),0,0));
        
        ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
            sen_set.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);
              

      
        
        