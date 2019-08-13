function fracal_v01

%% User inputs
a.Nv = 100;        % Number of verticla grids
a.Nh = 50;        % Number of horizontal grids
a.Vu = 1;         % Upper boundary voltage
a.Vl  = 0;        % Lower boundary voltage
a.delta = 0.001;  % potential calculation escaping value
a.MaxIter = 10000;       % Maximum itterating counts (when calculating potentials)


%% Program starts
V = zeros(a.Nv,a.Nh);

% Find potential everywhere
V = find_voltages(V,a);

% Initial point
xInds = [4 5];round(a.Nh/2);
yInds = [4 7];;

P = findProbs(V,a,xInds,yInds);


function P = findProbs(V,a,xInds,yInds)


L = length(xInds);
dE = zeros(L,8);

for i = 1:L
    % 1 top
    try; dE(i,1) = abs(V(yInds(i),xInds(i))- V(yInds(i)-1,xInds(i))); end 
    % 2 bottom
    try; dE(i,2) = abs(V(yInds(i),xInds(i))- V(yInds(i)+1,xInds(i))); end
    % 3 left 
    V(yInds(i),xInds(i))
    V(yInds(i),xInds(i)-1)
    try; dE(i,2) = abs(V(yInds(i),xInds(i))- V(yInds(i),xInds(i)-1)); end
    % 4 right
    try; dE(i,2) = abs(V(yInds(i),xInds(i))- V(yInds(i),xInds(i)+1)); end
    % 5 top left
    try; dE(i,2) = abs(V(yInds(i),xInds(i))- V(yInds(i)-1,xInds(i)-1))/1.4142; end
    % 6 top right
    try; dE(i,2) = abs(V(yInds(i),xInds(i))- V(yInds(i)-1,xInds(i)+1))/1.4142; end
    % 7 bottom left
    %abs(V(yInds(i),xInds(i))- V(yInds(i)+1,xInds(i)-1)/1.4142)
    try;dE(i,2) = abs(V(yInds(i),xInds(i))- V(yInds(i)+1,xInds(i)-1))/1.4142; end
    % 8 bottom right 
    try; dE(i,2) = abs(V(yInds(i),xInds(i))- V(yInds(i)+1,xInds(i)+1))/1.4142; end    
end

P = dE


%% Find volgates everywhere
function V = find_voltages(V,a)

    delta = inf;
    omega = 1.6;
    
    % Fix the top boundary
    V(1,:) = a.Vu; 
        
    % Fix the bottom boundary;
    V(end,:) = a.Vl;
    Vn = V;
    N =  0;
    
    while delta > a.delta && N < a.MaxIter
        N = N + 1;
       
        % Find deltas for middle chunk
        I = 2:a.Nh-1;   J = 2:a.Nv-1; 
        delta1 = V(J,I+1)+Vn(J,I-1)+V(J+1,I)+Vn(J-1,I)-4*V(J,I);
        % Find deltas for left raw
        I = 1; 
        delta2 = V(J,I+1)+Vn(J,a.Nh)+V(J+1,I)+Vn(J-1,I)-4*V(J,I);
        % Find deltas for right raw
        I = a.Nh; 
        delta3 = V(J,1)+Vn(J,a.Nh)+V(J+1,I)+Vn(J-1,I)-4*V(J,I);
        
        del = [delta2 delta1 delta3];
        
        Vn(J,:) = V(J,:) + omega/4*del;
        V = Vn;
               
        delta = max(max(del));
    end
