function long_range_lp_fm(t1,t2)

if nargin < 3
    t1 = 80000;
    t2 = 81000;
end


h=guidata(findall(0,'Tag','plotter2'));
settings = h.sen_set;
g = h.g;



%LDAR LINET PBFA marker size
mz = str2double(settings.xt_marker_size);

% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});


%Is manual Time shift for each sensor activated
if settings.man_tshiftOn==1
    tshift=settings.t_shift;
else
    tshift=zeros(1,60);
end


% Location data folder for local access
loc_dir = settings.loc_dir;

vshift=settings.vshift;
gain=settings.gain;
factor = g.factor;

% Generating LP and CH Legend names
lp_legend=cell(1,60);

for i=1:20;
    lp_legend{i*3-2}=[settings.sen_IDs{i} ':lp1'];
    lp_legend{i*3-1}=[settings.sen_IDs{i} ':lp2'];
    lp_legend{(i*3)}=[settings.sen_IDs{i} ':lp3'];
end

% the variable for real legend
lg={};

% checking witch graphs are on
lp_on=g.lpgraphs;

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

figure

for i=1:60
    % for loop used becuase there may be 15 lp plots
    if strcmp(a.lp_fn{i},'')==0
        [tn,vn]=SA_Extract1(a.lp_fn{i},g.t1,g.t2,tshift(i)-tshift1(ceil(i/3)));
        %soundsc(vn,10000,8,[-1 1])
        %[tn,vn,filtered] = lp_high_pass(tn,vn,200);
        
        %vn = remove_60Hz(tn,vn);
        
        if strcmp(settings.remove_time_decay,'on') 
            vn = removeDecay(vn,tDecay(mod(i,3)+1),tn(2)-tn(1));
            lg=[lg sprintf('%s C', lp_legend{i})];
        else
            lg=[lg lp_legend{i}];
        end
        
        plot(tn,vn*gain(i)*factor(i)+vshift(i))
        
        empty_plot=false;
        try
            waitbar(0.25+i*0.005,wbh,'Loading slow antenna data','Name','Plotter Busy')
        catch
            return
        end
    end
end

