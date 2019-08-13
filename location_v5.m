function location_v5()

%This program is under develop to find location of lightning

%Sensor Positions K02 K14 K24 WSB BCC respectively
x=[-7524.2     3371.6    -3148.5    -4266.3     -11657.9];
y=[16554.5     4445.6    -6837.6    9544.9      -17019.6];
z=[0            0        0          0       0];

%Speed of light in air
v=299702547;

%Take t from the graph

% K02 time
button = questdlg('Input the time for K02 sensor',...
    'Find Location','OK','No K02 data','Cancel','OK');

if strcmp(button,'Cancel') || strcmp(button,'')
    return
    
elseif strcmp(button,'No K02 data')
    t(1)=NaN;
else
    [t(1),k]=ginput(1);
    %disp(t(1))
end


% K14 time
button = questdlg('Input the time for K14 sensor',...
    'Find Location','OK','No K14 data','Cancel','OK');

if strcmp(button,'Cancel') || strcmp(button,'')
    return
    
elseif strcmp(button,'No K14 data')
    t(2)=NaN;
else
    [t(2),k]=ginput(1);
    %disp(t(2))
end


% K24 time
button = questdlg('Input the time for K24 sensor',...
    'Find Location','OK','No K24 data','Cancel','OK');

if strcmp(button,'Cancel') || strcmp(button,'')
    return
    
elseif strcmp(button,'No K24 data')
    t(3)=NaN;
else
    [t(3),k]=ginput(1);
    %disp(t(3))
end

% WSB time
button = questdlg('Input the time for WSB sensor',...
    'Find Location','OK','No WSB data','Cancel','OK');

if strcmp(button,'Cancel') || strcmp(button,'')
    return
    
elseif strcmp(button,'No WSB data')
    t(4)=NaN;
else
    [t(4),k]=ginput(1);
    %disp(t(4))
end


% K02 time
button = questdlg('Input the time for BCC sensor',...
    'Find Location','OK','No BCC data','Cancel','OK');

if strcmp(button,'Cancel') || strcmp(button,'')
    return
    
elseif strcmp(button,'No BCC data')
    t(5)=NaN;
else
    [t(5),k]=ginput(1);
    %disp(t(5))
end


% off set time
t0=min(t);
t=t-t0;



for i=1:5
    
    switch i
        case 1
            t1=t(1); t2=t(2); t3=t(3); t4=t(4);
            x1=x(1); x2=x(2); x3=x(3); x4=x(4);
            y1=y(1); y2=y(2); y3=y(3); y4=y(4);
            z1=z(1); z2=z(2); z3=z(3); z4=z(4);
            
        case 2
            t1=t(1); t2=t(3); t3=t(4); t4=t(5);
            x1=x(1); x2=x(3); x3=x(4); x4=x(5);
            y1=y(1); y2=y(3); y3=y(4); y4=y(5);
            z1=z(1); z2=z(3); z3=z(4); z4=z(5);
            
        case 3
            t1=t(1); t2=t(2); t3=t(4); t4=t(5);
            x1=x(1); x2=x(2); x3=x(4); x4=x(5);
            y1=y(1); y2=y(2); y3=y(4); y4=y(5);
            z1=z(1); z2=z(2); z3=z(4); z4=z(5);
        
        case 4
            t1=t(1); t2=t(2); t3=t(3); t4=t(5);
            x1=x(1); x2=x(2); x3=x(3); x4=x(5);
            y1=y(1); y2=y(2); y3=y(3); y4=y(5);
            z1=z(1); z2=z(2); z3=z(3); z4=z(5);
            
        case 5
            t1=t(2); t2=t(3); t3=t(4); t4=t(5);
            x1=x(2); x2=x(3); x3=x(4); x4=x(5);
            y1=y(2); y2=y(3); y3=y(4); y4=y(5);
            z1=z(2); z2=z(3); z3=z(4); z4=z(5);
            
        otherwise
            disp('Error Occured')
            return
    end
        
    
    K=[ x1-x2    y1-y2   t1- t2
        x2-x3    y2-y3   t2- t3
        x3-x4    y3-y4   t3- t4
        ];


    gg=[ x1^2+y1^2+z1^2-x2^2-y2^2-z2^2-v^2*(t1^2-t2^2)
        x2^2+y2^2+z2^2-x3^2-y3^2-z3^2-v^2*(t2^2-t3^2)
        x3^2+y3^2+z3^2-x4^2-y4^2-z4^2-v^2*(t3^2-t4^2)
        ]/2;
    
    %Find x,y and T ---- X(1) is x, X(2) is y and X(3) is related to T
    X=K\gg;
    xxx(i)=X(1);
    yyy(i)=X(2);    
    ttt(i)=-X(3)/v^2;    
    zzz(i)=abs(real(sqrt(v^2*(t(1)-ttt(i))^2-(x(1)-X(1))^2-(y(1)-X(2))^2)));

end




%% Generating file name for LDAR data

global g

settings=open('sensor_setting.mat');

if g.mm < 30
    ext=0;
else
    ext=30;
end

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
    settings.base_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,g.hh,ext);



% Loading Ldar data
[CG,CAL,DLS]=ldarExtract(ldar_fn,g.t1,g.t2,str2double(settings.ldar_r),...
    0,0,0,0);

t1=[DLS(:,9);CG(:,9)];
x1=[DLS(:,6);CG(:,6)];
y1=[DLS(:,7);CG(:,7)];
z1=[DLS(:,8);CG(:,8)];

%fprintf('Times\n')
Times=t+t0
    

fprintf('\nt\t\tx\t\ty\t\tz\n')

for i=1:5
    fprintf('%.6f\t%.1f\t\t%.1f\t\t%.1f\n',ttt(i)+t0,xxx(i),yyy(i),zzz(i))
end

fprintf('%.6f\t%.1f\t\t%.1f\t\t%.1f\n\n',mean(ttt)+t0,mean(xxx),mean(yyy),mean(zzz))

for i=1:length(x1)
    fprintf('%.6f\t%.1f\t\t%.1f\t\t%.1f\n',t1(i),x1(i),y1(i),z1(i))
end

figure
plot3(xxx(1),yyy(1),zzz(1),'*');

hold all
plot3(xxx(2),yyy(2),zzz(2),'*');
plot3(xxx(3),yyy(3),zzz(3),'*');
plot3(xxx(4),yyy(4),zzz(4),'*');
plot3(xxx(5),yyy(5),zzz(5),'*');
plot3(mean(xxx),mean(yyy),mean(zzz),'*');

if isempty(x1)==0
    plot3(x1,y1,z1,'p')
    legend('1','2','3','4','5','M','LDAR')
else
    legend('1','2','3','4','5','M')
end
    
hold off


box on
grid on