function test_NBP_auto_modeler1
tic
%% User inputs

% Real data folder
rdf = 'C:/Users/sumedhe/Desktop/Eack2004-data/';

% Model ID (which will identify the files)
mID = '01';


% Initial parameters
maxA=-100;
H1=7500;
x0=0;
y0=0;
t_step=1.0000e-006;
dh=100;  

% Ranges of values
t1s = 9e-6:4e-6:15e-6;
t2s = 60e-6:4e-6:80e-6;
vs = 6.0e7:2e7:10e7;
lamda1s = 100:500:2100;
lamda2s = 100:500:2100;
H2s    = 8000:500:9500;
Hms    = 7500:500:9500;
alpInds = 1:5:21;

sns = [11,12];

% Computers to use
% cUse = {'sadaq1' 'sadaq2' 'sadaq3' 'sadaq5' 'sadaq6' 'jconnComp' ...
%          'sadaq9' 'sadaq10' 'sadaq11' 'sadaq12' 'sadaqFLT', 'sadaq7'};

    
cUse = {'sadaq7'};


%% Starting program
% sensor locations
sns_x =  1.0e+04 *[ -0.7524, 0.3372, -0.3149, -0.4266, -1.1658, -0.2701, -2.7640, ...
                 -6.0091, 0.1825, -5.7394, -2.0637,2800,200000];
             
sns_y =    1.0e+04 *[1.6555, 0.4446, -0.6838, 0.9545, -1.7020, 0.2631, 4.9254, ...
                    -3.3983, -5.3008, 1.1923, 0.1569,0,0];

sns_IDs = {'K02','K14','K24','WSB','BCC','K17','EDW','STC','FLT','OVD','NEA','FAR'};


sns_x = sns_x(sns);
sns_y = sns_y(sns);
sns_IDs = sns_IDs(sns);

% Load real data
L0 = length(sns);

for i = 1:L0
    rd(i) = open([rdf sns_IDs{i} '.mat']);
end



[t1grid, t2grid, vgrid, lamgrid1,lamgrid2,H2grid,Hmgrid,alphIndgrid] = ...
    ndgrid(t1s,t2s,vs,lamda1s,lamda2s,H2s,Hms,alpInds);

%tdiff = nan(size(t1grid));

N = numel(t1grid);

% Get the range of the values to evaluate
[N1, N2, cID] = get_N_range(N,cUse);

% Exit if the computer not found
if cID < 0; return; end

% Today in DD-MMM-YYYY format
startday = date;

L = N2 - N1 + 1;

% varibles to find differences
diffs = nan(size(t1grid));
scale = diffs;
norm_diffs = diffs;



% Open matlabpool
if matlabpool('size') == 0 
    matlabpool
else
    tstr = datestr(now,'yyyy-mm-dd.HH:MM:SS');
    disp([tstr ' - matlabpool enabled!'])
end

     

% add Number of steps to complte to the progess file
fid = fopen('F:/parfor_progress.txt','w');
fprintf(fid,'%d\n',L);
fclose(fid);

% Zero percemt starting
fprintf(repmat(' ',1,15))

%% Main loop
parfor i = N1:N2
    
   progress(L)
 
   % Generate data for 10 sensors
   data = generate_data(i,maxA,t1grid,t2grid,vgrid,H1,H2grid,Hmgrid,alphIndgrid,x0,y0,...
           t_step,lamgrid1,lamgrid2,sns_x,sns_y,dh,rd,L0); 
   
   
   % Compare modeled data with real data
   if ~isempty(data)
       [diffs(i), norm_diffs(i), scale(i)] = compare_data(rd, L0,data);
   end
 
end


%% Saving the output
output.a.maxA = maxA;
output.a.H1 = H1;
output.a.x0 = x0;
output.a.y0 = y0;
output.a.t_step = t_step;
output.a.dh = dh;

% range parameters
output.b.t1s = t1s;
output.b.t2s = t2s;
output.b.vs = vs;
output.b.lamdas = lamda1s;
output.b.lamdas = lamda2s;
output.b.H2s = H2s;
output.b.Hms = H2s;
output.b.alpInds = alpInds;
output.b.sns = sns;
output.b.sns_x = sns_x;
output.b.sns_y = sns_y;

% Computer specific parameters
output.c.N  = N;
output.c.N1 = N1;
output.c.N2 = N2;
output.c.cID = cID; % Computer ID
output.date = startday;


% Results
output.c.diffs = diffs;
output.c.norm_diffs = norm_diffs;
output.c.scale = scale;



save_fn = [rdf 'ParameterSweeping-' ...
             startday '-' mID '-' cID '.mat'];

save(save_fn,'-Struct','output');

toc


   
   
function progress(N)

    fid = fopen('F:/parfor_progress.txt','a');
    fprintf(fid,'1\n');
    fclose(fid);
    
    fid = fopen('F:/parfor_progress.txt','r');
    is = fscanf(fid, '%d');
    fclose(fid);
    
    val = (length(is)-1)/N*100;
    val = sprintf('%06.2f',val);

    %fprintf([ repmat('\b',1,16) '%06.2f %% done!'],val)
    disp([repmat(char(8),1,15) val ' % done!'])
    
function [N1, N2, cID] = get_N_range(N,cUse)
        
        % Get the computer name
        [~, cname] = system('hostname');
        cname = strtrim(cname);
        
        if strcmp(cname,'sadaq-server')
            tstr = datestr(now,'yyyy-mm-dd.HH:MM:SS');
            disp([tstr ' - This is the server - Will run a test'])
            N1 = 1; N2 = 10; cID = 'Server';
            return
        end
            
        
        % Number of computers
        Nc = length(cUse);
        
        step = round(N/Nc);
        
        % computer number
        cNum = find(ismember(cUse,cname)==1);
       
        
        if isempty(cNum)
            tstr = datestr(now,'yyyy-mm-dd.HH:MM:SS');
            disp([tstr ' - Computer not found!'])
            N1 = 0; N2 = 0; cID = -1;
        else
            
            cID = cUse{cNum};
            N1 = (cNum -1)*step + 1;
            N2 = cNum * step;
            
            if N2 > N
                N2 = N;
            end
            tstr = datestr(now,'yyyy-mm-dd.HH:MM:SS');
            fprintf('%s - This is %s. \n',tstr, cname)
            fprintf('%s - Will run from %d to %d (%d steps)\n',...
                tstr,N1,N2,(N2-N1));
        end
            
    
function [diff, norm_diff, scale] = compare_data(rd, L0,data)
  
   diff = 0;
   norm_diff = 0;
   scale = data(1).scale;
   for j = 1:L0
              
       % Get calculated data for real data
       caly = interp1(data(j).t,data(j).y,rd(j).t)*scale;
       
%        figure; hold all;       
%        plot(rd(j).t,caly)
%        plot(rd(j).t,rd(j).y)
       
       % differnce
       diff = diff + sum(abs(caly - rd(j).y));
       
       factor = 1/range(rd(j).y);
       
       norm_diff = norm_diff + sum(abs(caly - rd(j).y))*factor;

   end

function data = generate_data(i,maxA,t1grid,t2grid,vgrid,H1,H2grid,Hmgrid,alphIndgrid,x0,y0,...
           t_step,lamgrid1,lamgrid2,sns_x,sns_y,dh,rd,L0)
   
   scale = zeros(1,L0);
   
   for j=1:L0
       % Changing parameters
       alpha = alphIndgrid(i)/t2grid(i);
       
       % Pulse arrival time to the sensor location
       t0 = sqrt((x0-sns_x(j))^2+(y0-sns_y(j))^2+H1^2)/3e8;
       rdt = rd(j).t;
       rdy = rd(j).y;
       
       Hm = Hmgrid(i);
       if Hm > H2grid(i)           
           data = [];
           return          
       end
             
       [t,E_stat,E_ind,E_rad] = NBP_modeler(maxA,t1grid(i),t2grid(i),vgrid(i),H1,Hm,H2grid(i),x0,y0,...
           t0 - 100e-6, t0 + 200e-6 ,t_step,lamgrid1(i),lamgrid2(i),sns_x(j),sns_y(j),dh,alpha);
       
       E_tot = E_stat+E_ind+E_rad;
       
       [hpk1, dt1] = max(E_tot);
       [hpk2, dt2] = max(rdy);
       
       lpk1 = min(E_tot);
       lpk2 = min(rdy);
   
       scale(j) = mean([hpk2/hpk1,lpk2/lpk1]);
       
       
       % time shift
       t = t + rdt(dt2) - t(dt1);
       
       % save data
       data(j).t = t;
       data(j).y = E_tot;      

   end

  
   data(1).scale = mean(scale);

    
    
    
