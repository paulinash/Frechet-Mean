function p = default_params()
%DEFAULT_PARAMS Parameters for Morris–Lecar + slow variable (as in your script).
%
% Keep parameters in one place so experiments can override selectively.

p = struct();

p.Vk  = -84;
p.Vl  = -60;
p.Vca = 120;

p.V1 = -1.2; p.V2 = 18;
p.V3 = 12;   p.V4 = 17.4;

p.gk   = 8;
p.gl   = 2;
p.c    = 20;
p.gca  = 4;
p.gkca = 0.25;

p.y0 = 10;
p.phi = 0.23;
p.mu = 0.2;
p.eps = 0.005;

% Applied current I is set per simulation.
p.I = NaN;

end
