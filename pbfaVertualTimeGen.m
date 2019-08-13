function ret = pbfaVertualTimeGen(arg)
% clc
% Speed of light in air
v = 299792458/1.0003;

% % Position of the vertual pulse
% arg.x0 = 1000;
% arg.y0 = 1000;
% arg.z0 = 5000;
% 
% arg.inner = [ 1 2 3 6 11];      % Inner sensors
% arg.outer = [ 5 7 8 9 10];      % Outer sensors

sen_set = arg.sen_set;

% Calculating inner time
ret.t_in = sqrt((sen_set.x(arg.inner)-arg.x0).^2 + ...
                (sen_set.y(arg.inner)-arg.y0).^2 +...
                (sen_set.z(arg.inner)-arg.z0).^2)/v;
            
% Calculating inner time
ret.t_out = sqrt((sen_set.x(arg.outer)-arg.x0).^2 + ...
                (sen_set.y(arg.outer)-arg.y0).^2 +...
                (sen_set.z(arg.outer)-arg.z0).^2)/v;
            
            
    
    
