function [P, ki_sqrd] = point_dipole_1(xyz,sns,dEms,output)
% This function is written to find point 1-D dipole estimation when 
% E change values are known at 1 or more sensors. 
% 
% Modification History:
%   2014-09-04 Created by Sumedhe Karunarathne
%   2015-07-30 Copied from point_dipole program to change to find vertical
%   point dipole moment


%% User Inputs
% Location of the dipole
% a.x0 = 0;  a.y0 = 2000;  a.z0 = 5000;

if nargin < 3
    a.x0 = -9662.0; a.y0 = 14157.0; a.z0 = 6177.0;
    
    % Sensors Used
    %a.sns = [1 2 3 4];
    a.sns = [1 2];
    
    % Measured E-changes
    %a.dEm = [0.9524   35.1329   11.6020    6.5050]; % Values obtained for px = 100; py = 200; pz = 300; for test purpose
    %a.dEm = [0.9   35.1   11.6    6.5]; % Rounded Values obtained for px = 100; py = 200; pz = 300; for test purpose
    a.dEm = [2.7602 -1.1312];
else
    a.x0 = xyz(1); a.y0 = xyz(2); a.z0 = xyz(3);
    a.sns = sns;
    a.dEm = dEms;
end

if nargin <4
    output = 1;
end
    

%% Other inputs
a.k = 1/4/pi/8.854e-12;
a.sen_set = open('sensor_setting.mat');

%% Theoritical estimation of E-change values at given sensors if you know px py and pz
% px = 100; py = 200; pz = 300;
% dEt = calculate_dEt([px py pz],a)

%% Starts the program

% Test value for 3d dipole
x0 = 1;
options=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.001},'TolFun',1e-30,'TolX',1e-30);
f = @(x) myfun(x,a);
%[x temp exitFlag] = fsolve(f,x0,options);
x = fsolve(f,x0,options);

% Re calculate the estimation of dE
dEt = calculate_dEt(x,a);

% Calculated the ki-squire
ki_sqrd = get_ki_sqrd(a,dEt);


%% Display results

if output

    fprintf('\nKi-sqrd = %0.2e\n', ki_sqrd)
    fprintf('P \t\t= [%5.1f %5.1f %5.1f] C.m\n',x)
    fprintf('\nSensor\t Measured \tCalculated \t %% Diff\n')
    fprintf('----------------------------------------\n')
    for i = 1:length(a.sns);
        fprintf('%s\t%13.1f\t%10.1f\t%6.1f%%\n',...
            a.sen_set.sen_IDs{a.sns(i)},a.dEm(i),dEt(i),(a.dEm(i)-dEt(i))/a.dEm(i)*100);
    end
    fprintf('----------------------------------------\n\n')
end

P = x;


function F = myfun(x,a)
% % From Krehbiel_1979 Equation 4
dEt = calculate_dEt(x,a);

% Ki squired
F = (a.dEm - dEt).^2;

function dEt = calculate_dEt(x,a)

% Calculate values of dEs that should be observed for give point dipole
% % From Krehbiel_1979 Equation 4
% When assuming the dipole is vertical, only Pz is available.
R = sqrt((a.sen_set.x(a.sns)-a.x0).^2+(a.sen_set.y(a.sns)-a.y0).^2+(a.sen_set.z(a.sns)-a.z0).^2);
dEt = a.k*(2*x./R.^3 - 6*a.z0./R.^5.*((a.sen_set.z(a.sns)-a.z0)*x));

function ki_sqrd = get_ki_sqrd(a,dEt)

% Measurement errors
mErr = [0.32 0.31 0.31 NaN 0.18 0.29 0.53 0.16 0.13 0.17 NaN];

L0 = length(a.sns);
if L0 > 3; np = L0 - 3;
else       np = 1;
end

ki_sqrd = sum(((a.dEm - dEt)./mErr(a.sns)).^2)/np;