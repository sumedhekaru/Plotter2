function IBP_modeler_parameter_test
% User inputs
a=struct( ...
    'maxA' , {-420}        , ...
    't1' , {5.6e-006}    , ...
    't2' , {32.0e-006}     , ...
    'v' ,  {4.6e7}             , ...
    'H1' , {7600}           , ...
    'H2' , {8400}          , ...
    'x0' , {0}              , ...
    'y0' , {0}              , ...
    'T1' , {0.1500e-004}    , ...
    'T2' , {2.5e-004}     , ...
    't_step' , {1.0000e-007}, ...
    'lamda', {113}, ...
    'x' , {50000}            , ...
    'y' , {0}               ,...
    'dh', {10}              , ...
    'alpha' , 10/32e-6 ...
    );


p = 1;                      % Start time end time persentage from peak

para = 10e-6:1e-6:50e-6;
% para = 1000:500:50000; % Changing horizontal distance
% para = 1:50; % Changing alpha
% para = 20:20:1500;          % Changing lamda
% para = 0.4e-6:0.2e-6:30e-6;      % Changing t
% para = 2e6:2e6:3e8;               % Changing v
bf = 'C:\Users\sumedhe\Desktop\IBP_modeling\changing_t2_50000m\';

%% Starting the program

% save('C:\Users\sumedhe\Desktop\test.mat','-Struct','b')

% b = open('C:\Users\sumedhe\Desktop\test.mat');
% E_rad = b.E_rad;
% E_ind = b.E_ind;
% E_stat = b.E_stat;



fn = [bf 'changing_t1.txt'];

if exist(fn,'file')
    fID = fopen(fn,'a+');
else
    fID = fopen(fn,'a+');   
    
    % Write the header
    fprintf(fID,'a = struct( ... \n');
    fprintf(fID,'\t''maxA'' ,{%0.1f}\t, ...\n',a.maxA);
    fprintf(fID,'\t''t1'' ,{%0.7f}\t, ...\n',a.t1);
    fprintf(fID,'\t''t2'' ,{%0.7f}\t, ...\n',a.t2);
    fprintf(fID,'\t''v'' ,{%0.1f}\t, ...\n',a.v);
    fprintf(fID,'\t''H1'' ,{%0.1f}\t, ...\n',a.H1);
    fprintf(fID,'\t''H2'' ,{%0.1f}\t, ...\n',a.H2);
    fprintf(fID,'\t''x0'' ,{%0.1f}\t, ...\n',a.x0);
    fprintf(fID,'\t''y0'' ,{%0.1f}\t, ...\n',a.y0);
    fprintf(fID,'\t''T1'' ,{%0.7f}\t, ...\n',a.T1);
    fprintf(fID,'\t''T2'' ,{%0.7f}\t, ...\n',a.T2);
    fprintf(fID,'\t''t_step'' ,{%0.7f}\t, ...\n',a.t_step);
    fprintf(fID,'\t''lamda'' ,{%0.1f}\t, ...\n',a.lamda);
    fprintf(fID,'\t''x'' ,{%0.1f}\t, ...\n',a.x);
    fprintf(fID,'\t''y'' ,{%0.1f}\t, ...\n',a.y);
    fprintf(fID,'\t''dh'' ,{%0.1f}\t ...\n',a.dh);
    fprintf(fID,');\n\n');
    
    fprintf(fID,'t1\trad_maxt\trad_max\trad_mint\trad_min\trad_st_t\trad_end_t\tind_mint\tind_min\tind_st_t\tind_end_t\tstat_min\n\n');
end


wbh = waitbar(0,'Overall process');

for i = 1:length(para)
    
    waitbar(i/length(para),wbh)
    
   fn2 = sprintf('t1-%2.2i-%2.2i',floor(para(i)*1e6),(para(i)*1e6 - floor(para(i)*1e6))*100);
   %fn2 = sprintf('v-%9.9i',round(para(i)))
   %fn2 = sprintf('lamda-%4.4i',para(i)); 
    
   % a.x = para(i);
    % a.alpha = para(i)/a.t2;
    %a.lamda = para(i);
    %a.t1 = para(i);
    %a.v = para(i);
    a.t2 = para(i);
    
    [t,E_stat,E_ind,E_rad] = IBP_modeler(a);
    b.t = t;
    b.E_rad = E_rad;
    b.E_ind = E_ind;
    b.E_stat = E_stat;
    
    save([bf fn2 '.mat'],'-Struct','b');    
    
    figure
    hold all
    tools2fig
    plot(t,E_rad)
    plot(t,E_ind)
    plot(t,E_stat)
    plot(t,E_rad+E_ind+E_stat)
    legend('Rad','Ind','Stat','Total')
    xlabel('Time (s)')
    ylabel('\DeltaE (V/m)')
    title([sprintf('t1 = %2.2f ',para(i)*1e6) '\mus'])
    %title(sprintf('v = %8.0f m/s',para(i)))
    ylim([-0.02 0.05])
    box on
    
    %Radiation term statistics
    [rad_max ind1] = max(E_rad);
    [rad_min ind2] = min(E_rad);
    rad_maxt = t(ind1);
    rad_mint = t(ind2);
    
    
    st_ind = sum(E_rad(1:ind1) <= p/100*rad_max);
    end_ind = sum(E_rad(ind2:end) <= p/100*rad_min)+ind2;
    
    try; rad_st_t = t(st_ind); catch; rad_st_t = NaN; end
    try; rad_end_t = t(end_ind); catch; rad_end_t = NaN; end
    
    plot(rad_maxt,rad_max,'ro')
    plot(rad_mint,rad_min,'ro')
    plot(rad_st_t,0,'ro')
    plot(rad_end_t,0,'ro')
    
    % Induction term statistics
    [ind_min ind1] = min(E_ind);
    ind_mint = t(ind1);
    
    st_ind = sum(E_ind(1:ind1) >= p/100*ind_min); 
    end_ind = sum(E_ind(ind1:end) <= p/100*ind_min)+ind1;
    
    try; ind_st_t = t(st_ind); catch; ind_st_t = NaN; end
    try; ind_end_t = t(end_ind); catch; ind_end_t = NaN; end
    
    plot(ind_mint,ind_min,'rs')
    plot(ind_st_t,0,'rs')
    plot(ind_end_t,0,'rs')
    
    % Save the figure and close it
    saveas(gcf,[bf fn2 '.png'])
    saveas(gcf,[bf fn2 '.fig'])
    close(gcf)
    
   
    % Static term statistics
    stat_min  = min(E_stat);
    
    
    % Write data in to a text file
    fprintf(fID,'%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\n',...
        para(i),rad_maxt,rad_max,rad_mint,rad_min,rad_st_t,rad_end_t,...
        ind_mint,ind_min,ind_st_t,ind_end_t,...
        stat_min);
    
end
delete(wbh)
fclose(fID);



    

    


