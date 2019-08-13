function test_farfor

tic
parfor i=1:1000
    for j=1:200
        x = mean(mean(magic(10)));
    end
end
toc