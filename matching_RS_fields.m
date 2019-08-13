function matching_RS_fields

%% User inoputs
a.sns = [ 1 2 3 6];                       % Sensor numbers
a.dEs = [-67.9 -73.7 -101.4 -154.9];      % measured E-change

% Rs position
x0 = -12303;
y0 = 1324.2;
z0 = 0;

% Max +/- horizontal distance from the RS location
maxHd = 10000;

% Max and Min search altitudes
maxH = 15000;
minH = 10000;

% x z grid resolution 
dx = 1000;
dz = 500;


%% Starting

% searching grid
xs = x0-maxHd:dx:x0+maxHd;
ys = y0-maxHd:dx:y0+maxHd;
zs = minH:dz:maxH;

% Lengths
L1 = length(xs);
L2 = length(ys);
L3 = length(zs);
N0 = L1*L2*L3;

% storing results
kis = zeros(L1,L2,L3);
rslt(L1,L2,L3).dQ = [];
rslt(L1,L2,L3).dEs = [];

wbh = waitbar(0,'Calculating...','name','Calculating RS fields');
val = 0;

%a.sns = [1 2 6];
%a.dEs = [-11.304 -19.101 -50.207];   
a.printOut = 0;
a.p2 = [x0 y0 z0];  % Lower position [x y z] 

for i=1:L1
    for j=1:L2
        for k=1:L3
            try
                val = val+1;
                msg = sprintf('Calculating... (%0.3f%%)',val/N0*100);
                waitbar(val/N0,wbh,msg)
            catch
                return
            end
                        
            a.p1 = [xs(i)   	ys(j)	zs(k)];  % Higher position [x y z]                       
            argout = dQ_finder(a);
            
            kis(i,j,k) = argout.ki;     
            rslt(i,j,k).dQ = argout.dQ;
            rslt(i,j,k).dEc = argout.dEc;
            
        end
    end
end

delete(wbh);

[mki mInd] = min(kis(:));
[i j k] = ind2sub(size(kis),mInd);

x1 = xs(i);
y1 = ys(j);
z1 = zs(k);

dQ = rslt(i,j,k).dQ;
a.p1 = [x1 y1 z1];
dEc = rslt(i,j,k).dEc;
sen_IDs = {'K02', 'K14', 'K24', 'WSB', 'BCC', 'K17' ,'EDW','STC','FLT','OVD'};

fprintf('\nCalculated charge  \t\t= %0.4f C\n',dQ);
fprintf('Reduced ki-sqrd \t\t= %0.4f \n',mki);
H = a.p1(3)/1000;
fprintf('Charge Moment \t\t\t= %0.04f C km (Verticle Height = %0.2f km)\n',...
    2*dQ*H,H)
D = sqrt(sum((a.p1 - a.p2).^2))/1000;
fprintf('Charge Moment \t\t\t= %0.04f C km (Distance = %0.2f km)\n',...
    2*dQ*D,D)
fprintf('Index \t\t\t\t\t= [%i %i %i] (Search Matrix = [%i %i %i])\n',i,j,k,L1,L2,L3)
fprintf('Upper chrge loc. \t\t= [%0.1f %0.1f %0.1f]\n',x1,y1,z1)

fprintf('_____________________________________________________\n')
fprintf('Sensor\t\tMeasured\t\tCalculated\t\t%%diff\n')
fprintf('=====================================================\n')
for i = 1:length(a.sns)
    fprintf('%s\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
        sen_IDs{a.sns(i)},a.dEs(i),dEc(i),abs(a.dEs(i)-dEc(i))/abs(a.dEs(i)+dEc(i))*200);
end

