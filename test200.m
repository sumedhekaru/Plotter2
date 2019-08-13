function test200


try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

g = h.g;
sen_set = h.sen_set;

%tgf_time = 51574.5430;
yyyys = {'2015' '2015' '2015' '2015' '2015' '2015' '2015' '2015'};
mms = {'04' '04' '05' '05' '08' '08' '08' '08'};
dds = {'27' '27' '20' '20' '10' '19' '21' '21'};
tgf_times = [47869.7000 
            47950.8350
            63482.7100
            63482.2090
            77712.6590
            51574.5430
            62736.2530
            62854.0650];

for tshift = -300:300

   

    
    figure
    set(gcf,'units','inches','outerposition',[1 1 14 7],'paperposition',[0 0 17 9])
    sen_set.t_shift(16*3) = tshift;
    sen_set.t_shift(16*3-1) = tshift;
    sen_set.t_shift(16*3-2) = tshift;
    %save('sensor_setting.mat','-Struct','sen_set')
    
    for i = 1:length(tgf_times)
        g.YYYY = {yyyys{i}};
        g.MM = {mms{i}};
        g.DD = {dds{i}};
        [hms hh mm] = sec2hhmmss(tgf_times(i));
        g.hh = hh;
        g.mm = floor(mm/5)*5;
        g.ss = 0;
        g.t1 = tgf_times(i) - 0.1;
        g.t2 = tgf_times(i) + 0.1;
        
        
        subplot(2,4,i)
        hold all
        
        data = load_data(g,sen_set);
        
        % Title
        title(sprintf('%i s --- %0.3f', tshift, tgf_times(i)))
        
        
        box on
      
                
        if isfield(data,'lp_data') && ~isempty(data.lp_data(16*3-2).t)            
            plot(data.lp_data(16*3-2).t - tgf_times(i),data.lp_data(16*3-2).y)
        end
        
        if isfield(data,'ch_data') && ~isempty(data.ch_data(16*3-2).t)            
            plot(data.ch_data(16*3-2).t - tgf_times(i),data.ch_data(16*3-2).y)
        end
        
        plot([0 0]+tgf_times(i)-tgf_times(i),ylim,'k--')
       % xlim([tgf_times(i) - 0.1, tgf_times(i) + 0.1])
        
    end
     
    
    % Save figure
    savefig(sprintf('C:/Users/sumedhe/Desktop/TGF-search/time_shift_search3/%3.3i.fig',tshift));
    print('-dpng', '-r100', sprintf('C:/Users/sumedhe/Desktop/TGF-search/time_shift_search3/%3.3i.png',tshift));
    
    
    close(gcf)
   
end

