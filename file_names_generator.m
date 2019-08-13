function a = file_names_generator(g,settings)


% Base folder
bfolder=sprintf('%s/%s/%s/%s/', ...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:});

% Generating LP and CH Legend names
lp_legend=cell(1,60);
ch_legend=cell(1,60);


for i=1:20;
    lp_legend{i*3-2}=[settings.sen_IDs{i} ':lp1'];
    lp_legend{i*3-1}=[settings.sen_IDs{i} ':lp2'];
    lp_legend{(i*3)}=[settings.sen_IDs{i} ':lp3'];
    
    ch_legend{i*3-2}=[settings.sen_IDs{i} ':ch1'];
    ch_legend{i*3-1}=[settings.sen_IDs{i} ':ch2'];
    ch_legend{i*3}=[settings.sen_IDs{i} ':ch3'];
    
end

%% Generating file names for lp graphs and check those files are availble

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

% outputs
a.lp_legend = lp_legend;
a.ch_legend = ch_legend;

loc_dir = settings.loc_dir;

if g.pbfa ==1 || (settings.plot_tshiftOn && settings.plot_tshiftType == 2)
   
    pbfa_fn=sprintf('%s/PBFA/%s/%s/%s/pbfa_%s%s%s_%2.2i%2.2i.txt',...
        loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},...
        g.MM{:},g.DD{:},g.hh,ext);
    
    if exist(pbfa_fn,'file')~=0
        % Check whether the file is exists
        a.pbfa_fn=pbfa_fn;
    else
        % If file is not exist don't store the file name
        a.pbfa_fn='';
        % Absant file
        if g.pbfa 
            a.absant_fn=[a.absant_fn pbfa_fn];
        end
   end
else
    a.pbfa_fn='';
end


