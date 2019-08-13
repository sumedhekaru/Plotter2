function return_stroke_dQ
% The intension of thie program is to find return strokes charge transfer

%% User inputs
% Return stroke time
t = 77636.3501650;

% Ranges to try
H1 = 1000; H2 = 70000; dH = 100;
Q1 =  0.1; Q2 = 100;   dq = 0.1;


%% setting up
% Get current plotter2 data
h=guidata(findall(0,'Tag','plotter2'));
g = h.g;
settings = h.sen_set;
tshift=settings.t_shift;


% Turn of lp plots but turn on all ch3 plots
g.lpgraphs = zeros(1,60);
g.chgraphs = [0 0 1 0 0 1 0 0 1 0 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];


hhmmss = sec2hhmmss(t);
g.hh = str2double(hhmmss(1:2));
g.mm = floor(str2double(hhmmss(4:5))/5)*5;

[ch_fn h_fn] = generate_ch_fn(settings,g);

wbh2 = waitbar(0,'Collecting data... ');
%figure
%hold all
%lg = {};
% Check whether sensors have triggered data

sns = [];
dEs = [];

for j=1:60;
    waitbar(j/60,wbh2,sprintf('Collecting data... (%0.1f%%)',j/6*10))
    if ~isempty(ch_fn{j})
        sn = ceil(j/3);
        g.t1 = t + sqrt(settings.x(sn)^2 + settings.y(sn)^2)/3e8 - 40e-6;
        g.t2 = g.t1 + 80e-6;
        
        [tch ych ch_freq_str] = FA_Extract1(ch_fn{j},h_fn{j},g.t1,g.t2,...
            tshift(j),settings,j);
        
        ych = ych*g.factor(sn);
        
        % find dE
        [m1 ind1] = min(ych);
        [m2 ind2] = max(ych(1:ind1));
        
        dE = m2 - m1;
        
        if ~isempty(dE)
            sns = [sns sn];
            dEs = [dEs dE];
        end
            
        %if ~isempty(tch)
        %    plot(tch,ych)
        %    lg = [lg settings.sen_IDs{sn}];
        %    
        %end
        
    end
end
delete(wbh2)



% Obtain CGLSS x,y
% Closest ldar2 point
ldar_fn = generate_ldar_fn(settings,g);
[CG,CAL,DLS]=ldarExtract2(ldar_fn,t-100e-6,t+100e-6,str2double(settings.ldar_r),0,0,0,0);

ts=DLS(:,10);         % time
x0 = nan;

if ~isempty(ts)
    % find closest ldar point to time "t"
    [minVal,minInd] = min(abs(ts - t));
    x0 = CG(minInd,6);
    y0 = CG(minInd,7);    
end

if isnan(x0)
    disp('No LDAR found, can not continue.')
    return
end
%legend(lg)
%box on


%% Generate data for test
% x0 = 0;
% y0 = 0;
% k0 = 1/2/pi/8.85418782e-12;
% dEs  = k0*4./((settings.x(sns)-x0).^2+(settings.y(sns)-y0).^2+4500^2);



k0 = 1/2/pi/8.85418782e-12;
qs = Q1:dq:Q2;
Hs = H1:dH:H2;

L1 = length(Hs);
L2 = length(qs);
L3 = length(sns);

errors = nan(L1,L2);

for i = 1:L1
    for j=1:L2
        chi_squired = 0;
        for k=1:L3
           zigma =  4; % Sigma for field mills (read page 556 Lightning, rakov and uman)
           E_cal = k0*qs(j)*Hs(i)/((settings.x(sns(k))-x0)^2+(settings.y(sns(k))-y0)^2+Hs(i)^2)^1.5;
         
           chi_squired = chi_squired + ((E_cal - dEs(k))/zigma)^2;
           %waitforbuttonpress 
        end
        errors(i,j) = chi_squired;
    end
end

figure
contourf(errors)


[r,c]=find(errors==min(min(errors)));

Q = qs(c)
H = Hs(r)
chi =  min(min(errors))






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