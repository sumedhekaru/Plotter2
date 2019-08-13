function find_flash_before_after(txyz,R,dt,idt,use_xyz)

% This function will find the closest flash to a given txyz point.
% input:
%       txyz = [time x y z] format
%       R    = Radius window to consider (defulat 10 km)
%       idt  = ignore very close pulses in time in +/-idt window (defulat
%       20us)

% Get plotter2 data
try h=guidata(findall(0,'Tag','plotter2'));
catch; disp('Run plotter2 first!'); return;
end

g = h.g;
sen_set = h.sen_set;


% If user clicked ch data that will potentially freez the computer. Let's
% user ask about that
if sum(g.chgraphs) > 0
    ButtonName = questdlg('Do you really need to plot ch data?', ...
        'Potential Computer Freez', ...
        'Yes', 'No','No');
    
    if ~strcmp(ButtonName,'Yes')
        return
    end
end



% Location data folder for local access
loc_dir = sen_set.loc_dir;

if nargin < 1
    % Get txyz from the user
    txyz = [76454.08803	2269.9	-93685.8	12114.8	1056];
end

if nargin < 2;  R = 10000;  end
if nargin < 3;  dt = 10;    end
if nargin < 4; idt = [150e-6 150e-6]; end



t0 = txyz(1);
x0 = txyz(2);
y0 = txyz(3);
z0 = txyz(4);

% Generate the file name for LDAR2
hh = floor(t0/3600);
mm = floor((t0 - hh*3600)/60);
ss = t0 - hh*3600 - mm*60;


if mm < 30
    ext=0;
else
    ext=30;
end

dnum= (datenum(str2double([g.YYYY,g.MM,g.DD])) ...
    -datenum(str2double(g.YYYY),0,0));

ldar_fn=sprintf('%s/ldar2/%s/%s/%s/ldar2_%s%3.3i%2.2i%2.2i.txt',...
    loc_dir,g.YYYY{:},g.MM{:},g.DD{:},g.YYYY{:},dnum,hh,ext);

%Check whether the file is exists
if ~exist(ldar_fn,'file')~=0
    msgbox(sprintf('File not found.\n%s',ldar_fn),'Flash Finder')
    return
end

if sen_set.ldar_tshiftOn
    xx0 = sen_set.x(sen_set.ldar_tshift_sn);
    yx0 = sen_set.y(sen_set.ldar_tshift_sn);
    zx0 = sen_set.z(sen_set.ldar_tshift_sn);
else
    xx0 = 0 ; yx0 = 0; zx0 = 0;
end



% load +/-dt data
[CG,CAL,DLS]=ldarExtract2(ldar_fn,t0-dt,t0+dt,1000000,xx0,yx0,zx0,0);
% 1 - occuring time
% 6 - x
% 7 - y
% 8 - z
% 9 - Occuring time
% 10 - Detection time
% 11 - horizontal Distance to each pulse



% Find LDAR2 points within R
if use_xyz
    D1s = sqrt((DLS(:,6)-x0).^2+(DLS(:,7)-y0).^2++(DLS(:,8)-z0).^2);
else
    D1s = sqrt((DLS(:,6)-x0).^2+(DLS(:,7)-y0).^2);
end

inds = find(D1s <= R);
tLDAR = DLS(inds,9);
xLDAR = DLS(inds,6);
yLDAR = DLS(inds,7);
zLDAR = DLS(inds,8);
tLDARa = DLS(inds,10);
rLDAR = D1s(inds);

% Find CG points within R
if use_xyz
    D2s = sqrt((CG(:,6)-x0).^2+(CG(:,7)-y0).^2+(CG(:,8)-z0).^2);
else
    D2s = sqrt((CG(:,6)-x0).^2+(CG(:,7)-y0).^2);
end
    
inds = find(D2s <= R);
tCG = CG(inds,9);
xCG = CG(inds,6);
yCG = CG(inds,7);
zCG = CG(inds,8);
tCGa = CG(inds,10);
rCG = D2s(inds);


% Plot 
g.t1 = t0 - dt;
g.t2 = t0 + dt;
g.hh = hh;
g.mm = floor(mm/5)*5;
g.ss = 0;
plot_all5(g,0);
hold all

% Plot NBP location
AH = findall(gcf,'type','axes');

% t0a = t0 + sqrt((x0 - sen_set.x(sen_set.ldar_tshift_sn))^2+...
%     (y0 - sen_set.y(sen_set.ldar_tshift_sn))^2+...
%     (z0 - sen_set.z(sen_set.ldar_tshift_sn))^2)/299792458.0;

t0a = t0;

plot(AH(2),t0a,z0,'r*','markerfacecolor','r','markersize',8)

% Plot all near locations
plot(AH(2),tLDARa,zLDAR,'mo','markerfacecolor','m','markersize',2)
plot(AH(2),tCGa,zCG,'ms','markerfacecolor','m','markersize',5)

%% Nearest pulses -- Time
disp('Nearest pulses in time:')
% Last pulse happened before NBP
indx1 = find(tLDARa<(t0a-idt(1)));
indx2 = find(tCGa<(t0a-idt(1)));

if ~isempty(indx1)
    bf_t1 = tLDARa(indx1(end));
    bf_x1 = xLDAR(indx1(end));
    bf_y1 = yLDAR(indx1(end));
    bf_z1 = zLDAR(indx1(end));
    bf_r1 = rLDAR(indx1(end));
    bf_t1a = tLDARa(indx1(end));
else
    bf_t1 = -inf;
end

if ~isempty(indx2)
    bf_t2 = tCGa(indx2(end));
    bf_x2 = xCG(indx2(end));
    bf_y2 = yCG(indx2(end));
    bf_z2 = zCG(indx2(end));
    bf_r2 = rCG(indx2(end));
    bf_t2a = tCGa(indx2(end));
else
    bf_t2 = -inf;
end


if bf_t1 == -inf && bf_t2 == -inf
    fprintf('\tThere was no flash found within %0.1f km and before %0.3f s\n',R/1000,dt)
elseif bf_t1 > bf_t2
    fprintf('\tIC\t%13.6f\t%8.1f\t%8.1f\t%8.1f\t%6.3f\t%6.1f\n',...
        bf_t1,bf_x1,bf_y1,bf_z1,bf_t1-t0a,bf_r1)
    plot(AH(2),bf_t1a,bf_z1,'r+','markerfacecolor','r','markersize',8)
else
    fprintf('\tCG\t%13.6f\t%8.1f\t%8.1f\t%8.1f\t%6.3f\t%6.1f\n',...
        bf_t2,bf_x2,bf_y2,bf_z2,bf_t2-t0a,bf_r2)
    plot(AH(2),bf_t2a,bf_z2,'r+','markerfacecolor','r','markersize',8)
end

% Pulse happened after NBP

% Last pulse happened before NBP
indx1 = find(tLDARa>(t0a+idt(2)));
indx2 = find(tCGa>(t0a+idt(2)));

if ~isempty(indx1)
    bf_t1 = tLDARa(indx1(1));
    bf_x1 = xLDAR(indx1(1));
    bf_y1 = yLDAR(indx1(1));
    bf_z1 = zLDAR(indx1(1));
    bf_r1 = rLDAR(indx1(1));
    bf_t1a = tLDARa(indx1(1));
else
    bf_t1 = inf;
end

if ~isempty(indx2)
    bf_t2 = tCGa(indx2(1));
    bf_x2 = xCG(indx2(1));
    bf_y2 = yCG(indx2(1));
    bf_z2 = zCG(indx2(1));
    bf_r2 = rCG(indx2(1));
    bf_t2a = tCGa(indx2(1));
else
    bf_t2 = inf;
end


if bf_t1 == inf && bf_t2 == inf
    fprintf('\tThere was no flash found within %0.1f km and after %0.3f s\n',R/1000,dt)
elseif bf_t1 < bf_t2
    fprintf('\tIC\t%13.6f\t%8.1f\t%8.1f\t%8.1f\t%6.3f\t%6.1f\n',...
        bf_t1,bf_x1,bf_y1,bf_z1,bf_t1-t0a,bf_r1)
    plot(AH(2),bf_t1a,bf_z1,'r+','markerfacecolor','r','markersize',8)
else
    fprintf('\tCG\t%13.6f\t%8.1f\t%8.1f\t%8.1f\t%6.3f\t%6.1f\n',...
        bf_t2,bf_x2,bf_y2,bf_z2,bf_t2-t0a,bf_r2)
    plot(AH(2),bf_t2a,bf_z2,'r+','markerfacecolor','r','markersize',8)
end





