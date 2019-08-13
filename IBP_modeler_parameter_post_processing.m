function IBP_modeler_parameter_post_processing

bf = 'C:\Users\sumedhe\Desktop\IBP_modeling\changing_t2_50000m\';

create_movie = 1;
create_change_t1_plots = 1;

%% Creating the movie
if create_movie 
fns = dir([bf '*.mat']);
[~,idx] = sort([fns.datenum]);
fns = fns(idx);

aviobj = avifile([bf 'example1.avi'],'compression','None','fps',10,'quality',100);


   fig =  figure;
       tools2fig
   
for k=1:length(fns)
    b = open([bf fns(k).name]);
    cla
    hold all

    plot(b.t,b.E_rad,'r')
    plot(b.t,b.E_ind,'g')
    plot(b.t,b.E_stat,'b')
    plot(b.t,b.E_rad+b.E_ind+b.E_stat,'k')
    legend('Rad','Ind','Stat','Total')
    xlabel('Time (s)')
    ylabel('\DeltaE (V/m)')
    %title([sprintf('t1 = %2.2f ',t1s(i)*1e6) '\mus'])
    title(fns(k).name)
    ylim([-0.07 0.1])
    box on
    

    F = getframe(fig);
    aviobj = addframe(aviobj,F);

end

aviobj = close(aviobj);


for i = 1:length(fns)
   fn =  fns(i).name;
end
end

%% Create change t1 plots
if create_change_t1_plots
fID = fopen([bf 'changing_t1.txt'],'r');
data=textscan(fID,'%f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',20);
fclose(fID);

figure
hold all
box on
% plot(data{1}*1e6,(data{6}-data{7})*1e6);
% plot(data{1}*1e6,(data{6}-data{2})*1e6);
plot(data{1},-(data{6}-data{7})*1e6); % Total time
plot(data{1},-(data{6}-data{2})*1e6); % Rise time
plot(data{1},data{3}*2500); % First peak
plot(data{1},-data{5}*2500);  % Second peak

xlabel('t1 (\mus)')
ylabel('time(\mus)')
legend('Pulse Width','Rise time','First peak (x2500)','second peak (x2500)')
end

