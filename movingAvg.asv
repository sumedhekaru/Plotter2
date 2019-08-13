function y_out = movingAvg(y,n)

if length(y) <= n
    y_out = y;
    return
end
L = length(y);

% Semi length
sl = n/2-0.5;
yn = zeros(n,L);
y1 = [repmat(y(1),1,sl), y, repmat(y(end),1,sl)];
for i = 1:n
  yn(i,:) = y1(i:L-sl+i+1);  
end

y_out = mean(yn);
    
% figure
% plot(x,y,'ro-',x,y_out,'ks-')