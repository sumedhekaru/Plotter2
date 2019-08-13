function varargout = PBFA_semi_auto(varargin)
% PBFA_SEMI_AUTO MATLAB code for PBFA_semi_auto.fig
%      PBFA_SEMI_AUTO, by itself, creates a new PBFA_SEMI_AUTO or raises the existing
%      singleton*.
%
%      H = PBFA_SEMI_AUTO returns the handle to a new PBFA_SEMI_AUTO or the handle to
%      the existing singleton*.
%
%      PBFA_SEMI_AUTO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PBFA_SEMI_AUTO.M with the given input arguments.
%
%      PBFA_SEMI_AUTO('Property','Value',...) creates a new PBFA_SEMI_AUTO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PBFA_semi_auto_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PBFA_semi_auto_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PBFA_semi_auto

% Last Modified by GUIDE v2.5 07-Feb-2015 21:38:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PBFA_semi_auto_OpeningFcn, ...
                   'gui_OutputFcn',  @PBFA_semi_auto_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PBFA_semi_auto is made visible.
function PBFA_semi_auto_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PBFA_semi_auto (see VARARGIN)
clc
% get Plotter2 data
try 
    h=guidata(findall(0,'Tag','plotter2'));
    handles.g = h.g;
    handles.sen_set = h.sen_set;
catch; disp('Run plotter2 first!');
end



% Other data
handles.a.data_save_fn = 'PBFA_semi_auto_data.mat';
handles.a.ch_upsample_factors = [5 1 1]; % Upsampling factor for ch1-3
handles.a.window_size = 350e-6; % Size of a window for cross corr
handles.a.peak_peak_thr = 0.01; % Peak threshold
% Load data
handles = loadData(handles,1);


% colors
handles.a.colors = [                 0                   0   1.000000000000000
   1.000000000000000                   0                   0
                   0   1.000000000000000                   0
                   0                   0   0.172413793103448
   1.000000000000000   0.103448275862069   0.724137931034483
   1.000000000000000   0.827586206896552                   0
                   0   0.344827586206897                   0
   0.517241379310345   0.517241379310345   1.000000000000000
   0.620689655172414   0.310344827586207   0.275862068965517
                   0   1.000000000000000   0.758620689655172
                   0   0.517241379310345   0.586206896551724
                   0                   0   0.482758620689655
   0.586206896551724   0.827586206896552   0.310344827586207
   0.965517241379310   0.620689655172414   0.862068965517241
   0.827586206896552   0.068965517241379   1.000000000000000
   0.482758620689655   0.103448275862069   0.413793103448276
   0.965517241379310   0.068965517241379   0.379310344827586
   1.000000000000000   0.758620689655172   0.517241379310345
   0.137931034482759   0.137931034482759   0.034482758620690
   0.551724137931034   0.655172413793103   0.482758620689655
];

% Do a big data cross correlation (using full length ch1 data)
handles = cross_correlation(handles);

% Devide data in to chunks 
% (second argument specify the chunk length in seconds)
handles = devide_data(handles);

% Find peaks
handles = find_peak_data(handles);

% Choose default command line output for PBFA_semi_auto
handles.output = hObject;



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PBFA_semi_auto wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PBFA_semi_auto_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ccor_show.
function ccor_show_Callback(hObject, eventdata, handles)
d = handles.d;
figure
AX(1) = subplot(2,1,1);
hold all
box on
xlabel('Time (s)'); ylabel('E-change (V/m)'); title('Cross correlation')
lg = {};

for i = 1:20
    indx = i*3-2;
    if ~isempty(d.ch_data(indx).t)
        plot(d.ch_data(indx).t+d.lag(i),d.ch_data(indx).yf,'color',handles.a.colors(i,:))
        lg = [lg d.ch_legend{indx}];
    end
end
clc
legend(lg)

AX(2) = subplot(2,1,2);
hold all
box on
xlabel('Time (s)'); ylabel('E-change (V/m)'); title('Cross correlation')
for i = 1:d.Nw
    for j = 1:20
        indx = j*3 - 2;
        t = d.window(i).ch_data(indx).t;
        if ~isempty(t)
            plot(t + d.lag(j) + d.window(i).lag(j),d.window(i).ch_data(indx).yf ...
                ,'color',handles.a.colors(j,:))            
        end
    end
end

linkaxes(AX,'x')

    


% --- Executes on button press in ccor_use_location.
function ccor_use_location_Callback(hObject, eventdata, handles)

% handles.a.peak_peak_thr = 0.005;
% handles = find_peak_data(handles);
handles.a.max_time_diff = 5e-6;  % Maximum time difference to consider
calculate_PBFA(handles);


% --- Executes on button press in ccorr_manual.
function ccorr_manual_Callback(hObject, eventdata, handles)
handles = find_peak_data(handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
h = msgbox('Clossing... Please wait','PBFA Auto');
d = handles.d;
% Save current data before exit
save(handles.a.data_save_fn,'-Struct','d','-v6')
try
    close(h)
end

fprintf('Closed\n')
delete(hObject);

function handles = loadData(handles,type)

% If type = 1 then program will force to load NEW data according to
% plotter2 inputs. If not, it will attempt to load data set from hard
% drive.
g = handles.g;
if type == 1
    % Turn off unncessary plots
    g.linet = 0;
    g.ldar = 0;
    g.cglss = 0;
    g.pbfa = 0;
    g.nldn = 0;
    g.pbfa_old = 0;
    g.lpgraphs = zeros(1,60);
        
    % We need ch1 data for each ch2/ch3 data selection for cross
    % correlation. Let's turn on ch1 according to ch1 and ch2
    indx = 1:3:60;
    g.chgraphs(indx) = g.chgraphs(indx) | g.chgraphs(indx+1) | g.chgraphs(indx+2);
    
    % Load data
    d = load_data(g,handles.sen_set);
    handles.d = d;
    
    %figure 
    %hold all
    %for i = 1:60
    %    try
    %     plot(d.ch_data(i).y)
    %    end
    %end
    
    
else
    try
        h = msgbox('Loading data. Please wait','PBFA Auto');
        handles.d = load(handles.a.data_save_fn);
        close(h)
        
    catch
        fprintf('Coudn''Load PBFA_semi_auto_data.mat\nLoadin new data instead ...')
        handles = loadData(handles,1);
        close(h)
    end
end


function handles = cross_correlation(handles,force)

if nargin < 2
    force = 0;
end

% If already cross correlated, we don't need to do it again
d = handles.d;

try if d.cCorrelated && ~force ; return ;  end
end

% Let's find the middle distance sensor, Assuming it has the medium
% intensity, and filter data
for i = 1:60    
    try
        d.range(i) =  range(d.ch_data(i).y);
    catch
        d.range(i) = NaN;
    end
end

[dy, I] = sort(d.range(1:3:60));

% ignore FFI, EDW
dy([7,11]) = NaN;

ref_sns = I(round((20 - sum(isnan(dy)))/2));


% reference waveform
refch = ref_sns*3-2;

% Save reference
d.refch = refch;
d.refsns = ref_sns;



% Filter datadata
for i = 1:60
    [yf, yh] = hill_tra(d.ch_data(i).t,d.ch_data(i).y,5000);
    %length(yf)
    %yf = movingAvg(yf,3);
    %length(yf)
    d.ch_data(i).yf = yf;
    d.ch_data(i).yh = yh;
end

% Upscale data
%figure
%hold all

for i = 1:60
    % Upscale ch1 data by factor of n
    n = handles.a.ch_upsample_factors(1);
    if mod(i,3) - 1 == 0 && ~isempty(d.ch_data(i).y)
        d.ch_data(i).yu = interp(d.ch_data(i).yf,n);
        dt = (d.ch_data(i).t(end) - d.ch_data(i).t(1))/(length(d.ch_data(i).yu)-1);
        d.ch_data(i).tu = d.ch_data(i).t(1):dt:d.ch_data(i).t(end);
        %plot(d.ch_data(i).tu,d.ch_data(i).yu)
    end   
end
        

refy = d.ch_data(refch).yf;

% cross corr
for i = 1:20
    % ch1 id according to sns
    chId = i*3 - 2;
    
    if chId ~= refch
        if ~isempty(d.ch_data(chId).y)
            
            % Normalize data to compare
            compy = d.ch_data(chId).yf * d.range(refch)/d.range(chId);
            
            [acor, lag] = xcorr(refy,compy);
           
            [~,I] = max(abs(acor));
            d.lag(i) = lag(I)/1e6;  % Assume 1MHz sampling
            
            
        else
            d.lag(i) = 0;
        end        
    else
        d.lag(i) = 0;
    end
end

handles.d = d;


function handles = devide_data(handles)
tic
d = handles.d;
tw = handles.a.window_size;

T1 = handles.g.t1;
T2 = handles.g.t2;
wbh = waitbar(0,'Calculating correlations','name','Please wait');
% Number of windows 
Nw = ceil((T2 - T1)/tw);
%figure
%hold all
%clc
for i = 1:Nw
    waitbar(i/Nw)
    for j = 1:60
        t = d.ch_data(j).t + d.lag(ceil(j/3));
        if ~isempty(t)
            t1 = T1 + (i-1)*tw;
            t2 = t1 + tw;
            
            lol = sum(t <= t1)+1;
            ul = sum(t < t2);
            
            d.window(i).ch_data(j).t = t(lol:ul) - d.lag(ceil(j/3));
            %d.window(i).ch_data(j).y = d.ch_data(j).y(lol:ul);
            d.window(i).ch_data(j).yf = d.ch_data(j).yf(lol:ul);
            d.window(i).lol = lol;
            d.window(i).ul = ul;
            d.window(i).t1 = t1;
            d.window(i).t2 = t2;
        else
            d.window(i).ch_data(j).t = [];
            %d.window(i).ch_data(j).y = [];
            d.window(i).ch_data(j).yf = []; 
            
        end
    end
    
    % Less do another cross correlation for each window
    refy = d.window(i).ch_data(d.refch).yf;
    
    % cross corr
    for j = 1:20
        % ch1 id according to sns
        chId = j*3 - 2
        
        if chId ~= d.refch
            if ~isempty(d.ch_data(chId).y)
                
                % Normalize data to compare
                compy = d.window(i).ch_data(chId).yf * d.range(d.refch)/d.range(chId);
                
                [acor, lag] = xcorr(refy,compy);
                
                [~,I] = max(abs(acor));
                d.window(i).lag(j) = lag(I)/1e6;  % Assume 1MHz sampling
                
                
            else
                d.window(i).lag(j) = 0;
            end
        else
            d.window(i).lag(j) = 0;
        end
        
        %plot(d.window(i).ch_data(chId).t + d.window(i).lag(j) + d.lag(j),d.window(i).ch_data(chId).yf)
        
    end
    
end

d.Nw = Nw;
handles.d = d;
delete(wbh)


function handles = find_peak_data(handles)

d = handles.d;
figure
hold all
for i = 1:20
    chId = i*3 - 2;
    
    % finding negative peaks
    y = -d.ch_data(chId).yu;
    t = d.ch_data(chId).tu;
    if ~isempty(t)
        [pks,locs] = findpeaks(y,'Threshold',handles.a.peak_peak_thr);
        d.ch_data(chId).pksy = -pks;
        d.ch_data(chId).pkst = t(locs);        
        %plot(t+d.lag(i),-y,'color',handles.a.colors(i,:))
        %plot(t(locs)+d.lag(i),-pks,'o','color',handles.a.colors(i,:))
    end
end


% Find peaks in each window
for i = 1:d.Nw
    
    % Most number of peaks
    mnp = -inf;

    for j = 1:20
        chId = j*3 - 2;
        
        % Get the peaks in this sensor
        t = d.ch_data(chId).pkst + d.lag(j) + d.window(i).lag(j);
       
        if ~isempty(t)
            
            lol = sum(t < d.window(i).t1) + 1;
            ul = sum(t < d.window(i).t2);
            
            d.window(i).ch_data(chId).pksy = d.ch_data(chId).pksy(lol:ul);
            d.window(i).ch_data(chId).pkst = t(lol:ul) - d.lag(j) - d.window(i).lag(j);
                                  
            plot(d.window(i).ch_data(chId).t ,...
                d.window(i).ch_data(chId).yf,'color',handles.a.colors(j,:))
            plot(d.window(i).ch_data(chId).pkst,...
                d.window(i).ch_data(chId).pksy,'o','color',handles.a.colors(j,:))
            
            % Most number of peaks?
            if length(t(lol:ul)) > mnp
                mnp = length(t(lol:ul));
                d.window(i).most_peak_ch = chId;
            end 
        end
    end
end

%legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20')
handles.d = d;

function handles = calculate_PBFA(handles)
d = handles.d;
arg.outer = [];      % Outer sensors
arg.t_out = [];
arg.sen_set = handles.sen_set;
figure
hold all
clc
% Calculate PBFA points in each window
for i = 1:d.Nw
    % Most number of peaks
    sns = ceil(d.window(i).most_peak_ch/3);
    peakTs = d.window(i).ch_data(d.window(i).most_peak_ch).pkst +  ...
        d.lag(sns) + d.window(i).lag(sns);
   
    for ts = peakTs
        sns = ceil(d.window(i).most_peak_ch/3);
        Ts = ts - d.lag(sns) - d.window(i).lag(sns); 
        for j = 1:20
            if j ~= sns(1)
                
                peaksTs2 = d.window(i).ch_data(j*3 - 2).pkst +  ...
                        d.lag(j) + d.window(i).lag(j);
                    
                [dt, ind] = min(abs(peaksTs2 - ts));
                
                if  dt < handles.a.max_time_diff
                    Ts = [Ts, peaksTs2(ind) - d.lag(j) - d.window(i).lag(j)];
                    sns = [sns, j];
                end
            end
        end
        %length(sns)
        if length(sns) > 5
            arg.inner = sns;
            arg.t_in = Ts;
            
            % Method 5 answer
            arg.method = 5;
            [xs,ys,zs,t1]=pbfa_finder(arg);
            ki_sqrd = cal_ki_sqrd(t1,xs,ys,zs,arg,handles.sen_set);
            
            
            if ki_sqrd < 5
                fprintf('%12.6f\t%12.1f\t%12.1f\t%12.1f\t%12.1f\n', ...
                    t1,xs,ys,zs,ki_sqrd)
                plot(t1,zs,'ro')
            end
            
        end
            
        
    end
end

function ki_sqrd = cal_ki_sqrd(t,x,y,z,arg,sen_set)

% Ki -sqrd
tt1 = arg.t_in';

deg_free = length(arg.t_in);

tt2 = t+sqrt((x - sen_set.x(arg.inner)).^2 + ...
    (y - sen_set.y(arg.inner)).^2 + ...
    (z - sen_set.z(arg.inner)).^2)'/3e8;

ki_sqrd = 1/deg_free*sum(((tt1 - tt2)/0.2e-6).^2);



% 
% 
% arg.inner = sns;
% arg.t_in = ts;
% 
% 
% n = length(ts);
% 
% 
% % Method 5 answer
% arg.method = 5;
% [xs,ys,zs,t1]=pbfa_finder(arg);
% %ki_sqrd = cal_ki_sqrd(t1,xs,ys,zs,arg,settings);

% % Method 3 answer
% arg.method = 3;
% [xs2,ys2,zs2,t2]=pbfa_finder(arg);
% ki_sqrd2 = cal_ki_sqrd(t2,xs2,ys2,zs2,arg,settings);
% 
% if ki_sqrd < ki_sqrd2
%     arg.method = 5;
% else
%     arg.method = 3;
%     xs = xs2; ys = ys2; zs = zs2; t1 = t2;
%     ki_sqrd = ki_sqrd2;
% end
% 
% % Determine we are looking for a RS?
% if ~isreal(zs) || zs < 3000
%     arg.method = 6;
%     [xs,ys,zs,t1]=pbfa_finder(arg);
%     pulse_type = 0;
% end
% 
% 
% 
% if n >= 5 && ki_sqrd < 5
%     
%     Ip = Ip_cal2(xs,ys,zs,settings,sns,vHs);
%     
%     snsStr = sprintf('%i,',sort(sns));
%     
%     if pulse_type < 0
%         pulseID = -1000-n-10*arg.method;
%     else
%         pulseID = pulse_type*1000+n+10*arg.method;
%     end
%     
%     fprintf('%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%4.4i\t%0.1f\t%s\tNaN\tNaN\tNaN\tNaN\t%0.1f\n',...
%         index,t1,xs,ys,zs,pulseID,Ip,snsStr(1:end-1),ki_sqrd)
%     
%     fprintf(fID,'%i\t%.7f\t%0.1f\t%0.1f\t%0.1f\t%3.3i\t%0.1f\t%s\tNaN\tNaN\tNaN\tNaN\t%0.1f\n',...
%         index,t1,xs,ys,zs,pulseID,Ip,snsStr(1:end-1),ki_sqrd);
%     
%     %plot(t1,zs/5000,'k*')
%     
%     store_peaks_info(index,[ts',sns'],fID2)
%     
%     index = index+1;
%     
% end



    





