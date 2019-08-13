function PBFA_convert
% There were 2 files generated for automatic PBFA for begining pulses and
% return strokes. This program will take data in those two files and
% converted in to the dirrectory structure that plotter 2 can understand.
%

% file name for PBFA begining pulses
fn1 = 'C:\Users\sumedhe\Desktop\PBFA_auto_data_test_using_positive_peaks.txt';

% file name for PBFA cg positions
fn2 = 'C:\Users\sumedhe\Desktop\PBFA_CG_data.txt';


% file name for CGLSS
fn3 = 'E:\data\cglss\2011\08\KSCCGLSS20110814.dat';



% Output file folder
out_dir='C:/Users/sumedhe/Desktop/pbfa2/2011/08/14/';

if exist(out_dir,'dir')==0    
    mkdir(out_dir)
end

%% Loading PBFA IBP data
% % Open and read all data in the file
fid = fopen(fn1);
beg_d=textscan(fid,'%f %f %f %f %f %f %f %f','HeaderLines',1);
fclose(fid);

% % lets make last column as number 1 for IBP type
beg_d{8} = zeros(size(beg_d{8}))+1;

beg_d = cell2mat(beg_d);

%% Loading PBFA CG data
fid = fopen(fn2);
cg_d=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f','HeaderLines',1);
fclose(fid);

% % lets make last column as number 2 CG type and data compatible with beg_d
cg_dd{1} = cg_d{1}; % index
cg_dd{2} = cg_d{2}; % t
cg_dd{3} = cg_d{3}; % x
cg_dd{4} = cg_d{4}; % y
cg_dd{5} = zeros(size(cg_d{5})); % z
cg_dd{6} = cg_d{6}; % n
cg_dd{7} = cg_d{5}; % Ip
cg_dd{8} = zeros(size(cg_d{5})) + 2; % type

cg_dd = cell2mat(cg_dd);

% Combine IPB and CG data

data = sortrows([beg_d ; cg_dd],2);

% clear unwanted big data fields
%clear beg_d cg_d cg_dd


%% Loading CGLSS data
cg = CGLSS_extract(fn3,0,86400,1000000,0,0,0);
 
% Total number of colums = 17
% 1 : Occuring time
% 2 : x
% 3 : y
% 12: z
% 11: Current
% 4 : distance to x0,y0,z0 given
% 5,6,7 - Latitude
% 8,9,10 - Logitude
% 13 - detection time
% 17 - Total number of sensors used


% Let's put into files
t1 = -1800;
t2 = 0;

for k=1:48
    
%     cnt1 = 0; % Number of PBFA CG points
%     cnt2 = 0; % Number of PBFA IBP points
     cnt3 = 0; % Number of misidentified IBPs
    cnt4 = 0; % Overall counter
    
    t1 = t1+1800;
    t2 = t2+1800;
    
    lol = sum(data(:,2) < t1) + 1;
    ul  = sum(data(:,2) < t2);
    
    if ul > lol
        
        str=sec2hhmmss(t1);
        fn = [out_dir 'pbfa_20110814_' str(1:2) str(4:5) '.txt'];
        fID = fopen(fn,'w');
        
        for i=lol:ul
            
            % if the point is a CG just write it to the file
            if data(i,8) == 2
                %cnt1 = cnt1 +1;
                cnt4 = cnt4 +1;
                fprintf(fID,'%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
                    cnt4,data(i,2),data(i,3),data(i,4),data(i,5),data(i,6),data(i,7));
            else
                % Lets see this IPB is valid
                dt = abs(cg_dd(:,2) - data(i,2));
                              
                if min(dt) < 10e-6
                    cnt3 = cnt3 + 1;
                else
                    %cnt2 = cnt2+1;
                    cnt4 = cnt4+1;
                    fprintf(fID,'%i\t%0.7f\t%0.1f\t%0.1f\t%0.1f\t%i\t%0.1f\n',...
                        cnt4,data(i,2),data(i,3),data(i,4),data(i,5),data(i,6),data(i,7));
                end
                
            end
        end
        
        fclose(fID);
    end
end

disp(cnt3)







