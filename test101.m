function test100
%%%%%%%% catogarize NBP according to RADAR
% file name for data
fn = 'C:\Users\Sumedhe\Desktop\20110814-NBP_info2.xlsx';

%% 
sheet = 1;
xlRange = 'A3:AD307';

ndata = xlsread(fn, sheet, xlRange);

type = ndata(:,17);
inds = ndata(:,1);

% base folder for radar data
rbf = 'C:\Users\sumedhe\Desktop\Nadee\Nadee-official-sumpc\NBP- vertical radar scan\2011-08-14\';

for i = 237:237
           
    fn = sprintf('XY-%3.3i',i);
    
    copy_file2([fn '.png'],rbf,type(i))
    copy_file2([fn '.fig'],rbf,type(i))
    
    fn = sprintf('Vertical_scan-%3.3i',i);
    
    copy_file2([fn '.png'],rbf,type(i))
    copy_file2([fn '.fig'],rbf,type(i))
    
end
    
function copy_file(fn,bf,type)

[bf sprintf('%2.2i/',type) fn]

if exist([bf fn],'file')
    copyfile([bf fn],[bf sprintf('%2.2i/',type) fn])
    disp('Success')
else
    disp('nope')
end

function copy_file2(fn,bf,type)

%[bf sprintf('%2.2i/',type) fn]

if exist([bf fn],'file')
    copyfile([bf fn],[bf '/RadarTypes/' sprintf('%2.2i/',type) fn])
    disp('Success')
else
    disp('nope')
end