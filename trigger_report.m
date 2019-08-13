function trigger_report

sid = 'K02';
yyyy = 2011;
MM = 08;
DD = 14;


% Get the date input from the plotter2
try 
    h=guidata(findall(0,'Tag','plotter2'));
    g = h.g;
    sen_set = h.sen_set
catch; disp('Run plotter2 first!'); return
end

L = 0;

for hh = 0:23
    for mm = 0:5:55
        fprintf('%2.2i\t%2.2i\n',hh,mm)
        
        fn = sprintf('%s/%4.4i/%2.2i/%2.2i/%s_%4.4i%2.2i%2.2i_%2.2i%2.2i%2.2i', ...
            sen_set.base_dir,yyyy,MM,DD,sid,yyyy,MM,DD,hh,mm,0);
        
        hfn = [fn '.h1'];
        chfn = [fn '.ch1'];
        
        trigs = [];
        
        if exist(hfn,'file')
            % Load corresponding header file times
            fId = fopen(hfn, 'r');
            trigs  = fread(fId, inf, 'double');
            fclose( fId );
        else
            fprintf('%s not exist!\n',fn)
        end
        
        if ~isempty(trigs)
            for j=1:length(trigs)
                sFa = epp_load_trigfile_time(chfn,trigs(j));
                y=sFa.y_i;
                t=sFa.t_i;
                
                L = length(trigs) + L
                
                
            end            
        end
        
        
    end
end
        


