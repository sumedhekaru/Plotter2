function custom_save(ind)

% saving folder
sbf = 'C:\Users\Sumedhe\Desktop\NBP_project_2013\NBP\20110814_plots\';

% save first figure
saveas(gcf,sprintf('%sMeanPowers/MeanPower-%3.3i',sbf,ind),'fig')
print(gcf,sprintf('%sMeanPowers/MeanPower-%3.3i',sbf,ind),'-dpng')
delete(gcf)

saveas(gcf,sprintf('%sPowers/Powers-%3.3i',sbf,ind),'fig')
print(gcf,sprintf('%sPowers/Powers-%3.3i',sbf,ind),'-dpng')
delete(gcf)

clc

