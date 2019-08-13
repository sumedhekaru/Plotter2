function flash_finder_easy_post_processing

fn = 'C:\Users\sumedhe\Desktop\flash_catogories_test-20140814-75km.txt'; % File name


fID = fopen(fn);
data = textscan(fID,'%f %f %f %f %f %f %f %f %f %s','HeaderLines',2);
fclose(fID);

L = length(data{1});
type = zeros(L,1);

for i = 1:L
    
    if strcmp(data{10}(i),'CG')
        type(i) = 1;
    end
end

data{10} = data{9};
data = cell2mat(data);
data(:,10) = type;

data = sortrows(data,10);

nCG = sum(data(:,10));
nIC = L - nCG

dataIC = data(1:nIC,:);
dataCG = data(nIC+1:end,:);


figure
hist(dataCG(:,8),0.5:20.5)
title(sprintf('Number of RSs per flash (Total = %i CG flashes)',nCG))
ylabel('Number of flashes')
xlabel('Number of RSs')



figure
hist(dataCG(:,9),5:10:400)
title(sprintf('Number of LDAR2 per flash (Total = %i IC flashes)',nIC))
ylabel('Number of flashes')
xlabel('Number of LDAR2 points')



