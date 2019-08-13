function testjavacode
hashcode=1
newhashcode=0;
for i=1:10
%     newhashcode=newhashcode + hashcode + ((-1)^(i-1) *((i+1)/2)^2)
%     disp(newhashcode);
var1 = (-1)^(i-1);
var2 = ((i+1)/2)^2;
fprintf('var1 %f and var2 %f\n',var1,var2);
end
    