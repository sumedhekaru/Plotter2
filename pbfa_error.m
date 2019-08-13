function [dx,dy,dz,dt]=pbfa_error(x0,y0,z0,sns,t_accu,N,TOAM,EM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function will find the dx,dy,dz,dt according to given inputs. Inputs
%  are 
%     (x0,y0,z0)    - error position 
%     sns           - sensors used (ex sns = [1 2 3 5 12])
%     t_accu         - Time Resolution (Time accuracny)
%     TOAM          - Time Of Arrival Method Used (see pbfa_finder.m)
%     EM            - Error method used
%
% Created: 2011-10-04 Sumedhe Karunarathne
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Input arguments
%     x0      = 5000;
%     y0      = 5000;
%     z0      = 5000;
%     sns     = [1 2 3 4 5];  
%     t_accu  = 0.0000002;
%     N       = 100;
%     TOAM    = 1;
%     EM      = 2;
% %%% End of input arguments 

% Position of the vertual pulse
a1.x0 = x0;
a1.y0 = y0;
a1.z0 = z0;

a1.inner = sns;      % Inner sensors
a1.outer = [];       % Outer sensors (usage of outer sites obasalete at the time of writing)


% find vertual times
getT = pbfaVertualTimeGen(a1);

sen_set = open('sensor_setting.mat');

switch EM
    case 1
    % When time of arrivals have accuracy, we can rount time to match that
    % accuracy. When we do that recalculated position for PBFA is differ
    % from the acctual value. This method will consider the absolute value
    % of that difference as the error.
    
    a2.t_in  = round(getT.t_in/t_accu)*t_accu;
    a2.t_out = [];
    a2.inner = sns;
    a2.outer = [];
    a2.method = TOAM;
            
    [x1,y1,z1,t1]=pbfa_finder(a2);
    
    dx = abs(x1-a1.x0);
    dy = abs(y1-a1.y0);
    dz = abs(z1-a1.z0);
    dt = abs(t1);
    
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
    % Here we do not round the theoritical time of arrival obtained.
    % Instead we add random error to between +/- 1us (if 1us is the timing
    % accuracy). Then we calculate the position and find the error in the
    % position. We can do 1000s of calculation to obtain an error
    % distribution and then use the standard deviation as the error of the
    % each coordinate.
    
  
    n1 = length(getT.t_in);
    n2 = length(getT.t_out);
    
    a2.method = TOAM;
    a2.inner = sns;
    a2.outer = [];
          
    dx_temp = NaN(1,N);
    dy_temp = dx_temp;
    dz_temp = dx_temp;
    dt_temp = dx_temp;
    
    wbh = waitbar(0,'Ready for calculations','name','PBFA error propagation');
    tic
    for k = 1:N
        try
            str = sprintf('Working on %i of %i',k,N);
            waitbar(k/N,wbh,str)
        catch
            return
        end
        
        
        % Add randowm numbers between +/- 1/2 time_accuracy to the time
        [raw1 col1] = size(getT.t_in);
        
        %r1 = -t_accu + 2*t_accu.*rand(1,n1);
        r1 =  normrnd(0,t_accu,raw1,col1);
        %r2 = -t_accu/2 + t_accu.*rand(1,n2)

        
        a2.t_in  = getT.t_in + r1;
        a2.t_out = [];
        [x1,y1,z1,t1]=pbfa_finder(a2);
        
        dx_temp(k) = (x1-x0);
        dy_temp(k) = (y1-y0);
        dz_temp(k) = (z1-z0);
        dt_temp(k) = (t1);
        
    end
            
    dx = std(dx_temp);
    dy = std(dy_temp);
    dz = std(dz_temp);
    dt = std(dt_temp); 
 
    delete(wbh)

end

