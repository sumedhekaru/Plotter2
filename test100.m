ke_ae = 8494.7;
r = 150;
angs1 = 0:0.1:12;
angs = angs1 * pi/180;
L = length(angs);
dys = nan(1,L);

for i = 1:L
     
  
        % altitudes of the lowest layer
        %h1 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan(vD(k).a*pi/180);
        %h2 = sqrt((vx - x0/1000).^2 + (vy - y0/1000).^2)*tan((vD(k).a+vD(k+1).a)/2*pi/180);
        
        h1 = sqrt(r^2 + ke_ae^2 + 2*r*ke_ae*sin(angs(i)))-ke_ae;
        h2 = r*tan(angs(i));
        dys(i) = h1 - h2;
end

%figure
hold all
plot(angs1,dys)