function find_n(dE_max,dist)

% inorder to find peak current, Nadee is working on plotting 
%    Ip vs dE*r^n
% This n could be any value should be close to 1. This program is intend to
% find the best n value.

dE_max =[ ...
-46.64
-53.92
-25.92
-67.58
-9.96
-10.03
-13.14];

dist = [ ...
16499.2
12606.9
18310.3
10208.4
58962.3
56202.7
45836.9];

% n range
n_range = [1 1.3];

% incriments
inc = 0.01;

n_values = n_range(1):inc:n_range(2);
leng = length(n_values);
stds = nan(1,leng);

for i = 1:5
    temp = dE_max .* dist.^n_values(i)
    stds(i) = range(temp)/mean(temp);
end

[m ind] = min(stds)

    

