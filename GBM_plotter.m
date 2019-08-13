function GBM_plotter


%% USER INPUTS
% GBM CSV filename
fn = 'C:\Users\Sumedhe\Desktop\100807797_Merged.csv'; % Data file name
new_plot = 1;  % Make this 1 if you want a new plot, 0 - plot on current plot
yshift = -58.8;
ygain = 0.001;
tshift = 69435.813579 - 0.001171;


%% Plotting
fID = fopen(fn);
d = textscan(fID,'%f , %f , %f, %f, %f, %f','HeaderLines',9);
fclose(fID);

d = cell2mat(d);



binSize = 50e-6;

binMin = floor(min(d(:,4))/binSize)*binSize;
binMax = ceil(max(d(:,4))/binSize)*binSize;

bins = (binMin:binSize:binMax)+binSize/2;

[nB xB] = hist(d(:,4),bins);

if new_plot
    h = figure;
else
    h = figure(gcf);
    hold all
end

stairs(xB-binSize/2+tshift,nB*ygain+yshift)

%xlim([-400e-6 500e-6])

if new_plot
    xlabel('Time (s)')
    ylabel([sprintf('Counts per %0.1f',binSize*1e6) '\mus'])
    title('All 14 GBM Detectors, Channels 0 to 127')
    box on
    legend('GBM')
else
    
    lg = get(findobj(h,'Type','axes','Tag','legend'),'string');
    
    n1 = find(strcmp(lg, 'LINET'));
    n2 = find(strcmp(lg, 'PBFA'));
    n3 = find(strcmp(lg,'LDAR2'));
    n4 = find(strcmp(lg,'CGLSS'));
    
    n = min([n1 n2 n3 n4]);
    
    lgt1 = lg(1:n-1);
    lgt2 = lg(n:end);
    
    legend([lgt1 'GBM' lgt2])
end



