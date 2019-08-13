function leader_modeling1
% This will only work with a single sensor


% Load settings file
sen_set = open('sensor_setting.mat');

%% User inputs
% Sensor number to model
sn = [ 1 2 3 6];

% HSVP file name
hsvp_fn = 'C:\Users\sumedhe\Desktop\Leader Modeling\F18-14Aug_231217-50k8mm320x240.txt';



% E-chage data file
rd.t = [];
rd.v = [];
L0 = length(sn);
for i = 1:L0
    E_fn = sprintf('C:/Users/sumedhe/Desktop/Leader Modeling/%s_10KHz_leader_data.mat',sen_set.sen_IDs{sn(i)});
    % real data
    d = open(E_fn);
    rd.t = [rd.t; d.t];
    rd.v = [rd.v; d.v];
end

%rd.v(1,:) = smooth(rd.v(1,:));


%% Start the program

% Get HSVP data
b = hsvp_extract(hsvp_fn);
b.align_time  = 83537.3031492;
b.align_frame = -53497;
b.frame_rate  = 50000;

b =  scr2ldar(b);
k = 1/2/pi/8.85418782e-12;





[tuni endInd] = unique(b.t);
L = length(tuni);
E = zeros(L0,L-1);

ind = find(b.frameN == b.frameN(1));

% First value
%rd.v = smooth(rd.v);
for i = 1:L0
    vshift = interp1(rd.t(i,:),rd.v(i,:),b.t(ind(1)));
    rd.v(i,:) = rd.v(i,:) - vshift;
end

wbh = waitbar(0,'Working on...','name','Modeling leader');

totQ = 0;

dQs = zeros(1,L0);

for i = 2:L-1
    
    try
        waitbar(i/L,wbh,sprintf('Working on...(%0.2f%%)',i/L*100));
    catch
        return
    end
    
    ind = find(b.frameN == b.frameN(endInd(i)+1));
    
    for j = 1:L0
        
        x0 = sen_set.x(sn(j));
        y0 = sen_set.y(sn(j));
        
        % measured dE    
        dEm = interp1(rd.t(j,:),rd.v(j,:),tuni(i))-interp1(rd.t(j,:),rd.v(j,:),tuni(i-1));
        %dEm = interp1(rd.t(j,:),rd.v(j,:),tuni(i))-E(j,i-1);
        
        % find dQ
        dQs(j) = find_dQ(dEm,b,ind,endInd(i),x0,y0);
        
    end
    
    % charge 
    dQ = mean(dQs);
    totQ = totQ + dQ;
    
    for j = 1:L0
        
        x0 = sen_set.x(sn(j));
        y0 = sen_set.y(sn(j));
        
        % calculated dE
        dEc = sum(k*dQ*b.I(ind)/sum(b.I(ind)).*b.z(ind)./((b.x(ind)-x0).^2+(b.y(ind)-y0).^2+b.z(ind).^2).^1.5);    
        E(j,i) = E(j,i-1)+dEc;
    end
    
    
          
end


delete(wbh)


for i=1:L0
    figure
    hold all
    plot(rd.t(i,:),rd.v(i,:),'LineWidth',2)
    plot(tuni(1:L-1),E(i,:),'--','LineWidth',2)
    xlim([min(tuni) max(tuni)])
    
    xlabel('t (s)')
    ylabel('V/m')
    legend(sen_set.sen_IDs{sn(i)},sprintf('Modeled (%0.2f C)',totQ))
    box on
end





% figure
% plot3(x,y,z,'r*')
% daspect([1 1 1])
% box on
% xlabel('x (m)')
% ylabel('y (m)')
% zlabel('z (m)')
% 
% figure
% plot(t,z,'rp')



function b = scr2ldar(b)
N_x = b.N_x;
N_y = b.N_y;

f = 8.0e-3;  %Focal length
p = 20.0e-6; %Fixel Size


r = sqrt((((b.ref_x-b.cam_x)).^2+(b.ref_y-b.cam_y ).^2 ).*(1+(p^2.*(b.scr_x-N_x ).^2)./f^2 ));

theta = atan((b.ref_y-b.cam_y)./(b.ref_x-b.cam_x ))+atan(p*(b.scr_x-N_x)/f);

b.x = b.cam_x + r.*cos(theta);
b.y = b.cam_y + r.*sin(theta);
b.z = b.cam_z + r.*((b.scr_y-N_y).*p)./f;

% Convert frame numbers to time
b.t = b.align_time + (b.frameN - b.align_frame)./b.frame_rate;




% function dQ = find_dQ(dEm,b,ind,endInd,x0,y0)
% 
%     k = 1/2/pi/8.85418782e-12;
% 
%     % Level one
%     dQs = 0.01:0.01:5;
%     size(dQs);
%     L1 = length(dQs);
%     errors = zeros(1,L1);
%     
%     for ii = 1:L1;
%         dEcal = sum(-k*dQs(ii)*b.I(ind)/sum(b.I(ind)).*b.z(ind)./((b.x(ind)-x0).^2+(b.y(ind)-y0).^2+b.z(ind).^2).^1.5);
%         errors(ii) = abs(dEm - dEcal); 
%     end
%     
%     [mmm mInd] = min(errors);
%     
%     % Level two
%     if mInd == 1        
%         dQs = dQs(1)/100:dQs(1)/100:dQs(2);
%     else
%         dQs = dQs(mInd-1):(dQs(mInd+1)-dQs(mInd-1))/100:dQs(mInd+1);
%     end
% 
%     size(dQs);
%     L1 = length(dQs);
%     errors = zeros(1,L1);
%     
%     for ii = 1:L1;
%         dEcal = sum(-k*dQs(ii)*b.I(ind)/sum(b.I(ind)).*b.z(ind)./((b.x(ind)-x0).^2+(b.y(ind)-y0).^2+b.z(ind).^2).^1.5);
%         errors(ii) = abs(dEm - dEcal); 
%     end
%     
%     [mmm mInd] = min(errors);
%     
%     % Level three   
%     if mInd == 1        
%         dQs = dQs(1)/100:dQs(1)/100:dQs(2);
%     else
%         dQs = dQs(mInd-1):(dQs(mInd+1)-dQs(mInd-1))/100:dQs(mInd+1);
%     end
%     
%     size(dQs);
%     L1 = length(dQs);
%     errors = zeros(1,L1);
%     
%     for ii = 1:L1;
%         dEcal = sum(-k*dQs(ii)*b.I(ind)/sum(b.I(ind)).*b.z(ind)./((b.x(ind)-x0).^2+(b.y(ind)-y0).^2+b.z(ind).^2).^1.5);
%         errors(ii) = abs(dEm - dEcal); 
%     end
%     
%     [mmm mInd] = min(errors);
%     
%     dQ = dQs(mInd);
    
 function dQ = find_dQ(dEm,b,ind,endInd,x0,y0)
     k = 1/2/pi/8.85418782e-12;
     
dQ = dEm / sum(k*b.I(ind)/sum(b.I(ind)).*b.z(ind)./((b.x(ind)-x0).^2+(b.y(ind)-y0).^2+b.z(ind).^2).^1.5);

