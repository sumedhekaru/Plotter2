function find_positive_RS
clc
% Use CGLSS to find positive RSs

date = '20110814';
t1 = 72000;
t2 = 84600;
rc = 100000;
fn = 'C:\data\2011-08-05 -- 2011-08-16\data\cglss\2011\08\KSCCGLSS20110814.dat';



data=CGLSS_extract(fn,t1,t2,rc,0,0,0);



L = length(data);
k = 0;



for i = 1:L
    if data(i,11) > 0
        k = k + 1;
        fprintf('%4.4i\t%12.6f\t%10.1f\t%10.1f\t%10.1f\t%10.1f\n',...
            k,data(i,1),data(i,2)/1000,data(i,3)/1000,data(i,12)/1000,data(i,11))
    end
end

if k == 0
    disp('No positive RSs found')
end