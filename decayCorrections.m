function decayCorrections

% Slow antenna data has 10s (lp2) and 1s (lp1) electronic decay added to
% their data. The intention of this program to add corrections to lp2 and
% lp3 data.

%% Get plotter2 data
    try; h=guidata(findall(0,'Tag','plotter2'));
    catch; disp('Run plotter2 first!'); return;
    end
    
%% Start program
g = h.g;
settings = h.sen_set;
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

%Is manual Time shift for each sensor activated
if settings.man_tshiftOn==1
    tshift=settings.t_shift;
else
    tshift=zeros(1,60);
end

%Is offset voltages should be included?
if settings.vshiftOn==1
    vshift=settings.vshift;
else
    vshift=zeros(1,60);
end

%Is manual gain should be included?
if settings.gainOn==1
    gain=settings.gain;
else
    gain=zeros(1,60)+1;
end

% Is calibration should be included
if settings.plot_calibrated == 1
    factor = g.factor;
else
    factor = zeros(1,60)+1;
end

tshift1 = zeros(1,20);

% generate lp filenames

lp_on=g.lpgraphs;

lg = {};

%Generating LP and CH Legend names
lp_legend=cell(1,60);
ch_legend=cell(1,60);

% Variable to store tshift info
tshift_str =[];

% Witch sensors has time shift turn on?
tmp = g.lpgraphs + g.chgraphs;


% Do we need to time shift graphs to points?
if ~settings.plot_tshiftOn
    settings.plot_tshiftType = 0;
end

for i=1:20;
    lp_legend{i*3-2}=[settings.sen_IDs{i} ':lp1'];
    lp_legend{i*3-1}=[settings.sen_IDs{i} ':lp2'];
    lp_legend{(i*3)}=[settings.sen_IDs{i} ':lp3'];
    
    ch_legend{i*3-2}=[settings.sen_IDs{i} ':ch1'];
    ch_legend{i*3-1}=[settings.sen_IDs{i} ':ch2'];
    ch_legend{i*3}=[settings.sen_IDs{i} ':ch3'];
    
    if sum(tmp(i*3-2:i*3))> 0 && tshift(i*3)~= 0
        tmps = sprintf(' %s:%0.6fs', settings.sen_IDs{i},tshift(i*3));
        tshift_str = [ tshift_str tmps];
    end
end

% Save absant file extentions
a.absant_fn={};


for i=1:60
    if lp_on(i)==1
        % finding the file extention number
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
            a.lp_fn{i}=filename;
        else
            % If file is not exist don't store the file name
            a.lp_fn{i}='';
            % Absant file
            a.absant_fn=[a.absant_fn filename];
        end
    else
        a.lp_fn{i}='';
    end
end

%% Generating file names for ch (with header files) graphs and check those files are availble

% checking witch graphs are on
ch_on=g.chgraphs;

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
        
        % Check whether the ch file and h file are exists
        if exist(filename,'file')==0
            a.absant_fn=[a.absant_fn filename];
        end
        
        if exist(hfilename,'file')==0
            a.absant_fn=[a.absant_fn hfilename];
        end
        
        % If file is not exist don't store the file name
        if exist(filename,'file')==0 || exist(hfilename,'file')==0
            a.ch_fn{i}='';
            a.h_fn{i}='';
        else
            a.ch_fn{i}=filename;
            a.h_fn{i}=hfilename;
        end
    else
        a.ch_fn{i}='';
        a.h_fn{i}='';
    end
end



% Load and plot lp data


fg=figure;
hold all
box on

% Store one data set for plot yy

for i=1:60
    % for loop used becuase there may be 15 lp plots
    if strcmp(a.lp_fn{i},'')==0
        [tn,vn]=SA_Extract1(a.lp_fn{i},g.t1,g.t2,tshift(i)-tshift1(ceil(i/3)));
        %soundsc(vn,10000,8,[-1 1])
        %[tn,vn,filtered] = lp_high_pass(tn,vn,200);
        vmm = max(vn);
        plot(tn,(vn-vmm)*gain(i)*factor(i)+vshift(i))
        
        %tn = 0:0.00001:5;
        %vn = -5*exp(-tn/1);
        
        vn = removeDecay(tn,vn,lp_legend{i});
        plot(tn,(vn-vmm)*gain(i)*factor(i)+vshift(i))
        
        lg=[lg lp_legend{i}];
        lg=[lg lp_legend{i}];
        empty_plot=false;
    end
end


output.trigrd = [];

for i=1:60
    
   
    if strcmp(a.ch_fn{i},'')==0
        
        [t y ch_freq_str] = FA_Extract1(a.ch_fn{i},a.h_fn{i},g.t1,g.t2,tshift(i)-tshift1(ceil(i/3)),settings,i);
        
        if isempty(t)==0
            
            if isempty(output.trigrd)
                output.trigrd = sprintf('%i',ceil(i/3));
            else
                output.trigrd = sprintf('%s,%i',output.trigrd,ceil(i/3));
            end
              
            
            ch_legend_str = [ch_legend{i} ch_freq_str];
            lg=[lg ch_legend_str];
            lg=[lg ch_legend_str];
            
            mmmm = max(y);
            plot(t,(y-mmmm)*gain(i)*factor(i)+vshift(i))
            vn = removeDecay(t,y,ch_legend{i});
            plot(t,(vn-mmmm)*gain(i)*factor(i)+vshift(i))
            
            %soundsc(y,15000,8,[-1 1])

            empty_plot=false;
        end
    end
end

legend(lg')
   



    