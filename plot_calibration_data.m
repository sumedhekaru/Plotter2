function plot_calibration_data
% After auto calibrating, after keeping selected files, running of this
% program visualize the calibration data.

% File name with selected calibration data
fn = 'C:\Users\sumedhe\Desktop\Sumedhe\Calibration\Auto\20110801\FFI\ch3_cross\with_K14\selected\selected_ch1_cal_data.txt';

% Read the file
fID = fopen(fn,'r');
data = textscan(fID,'%s %f');
d = data{2};
figure
hist(d)
title(sprintf('(%0.3f\\pm%0.3f)',mean(d),std(d)))
fclose(fID);