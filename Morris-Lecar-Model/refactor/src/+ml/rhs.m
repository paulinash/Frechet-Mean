function dz = rhs(~, z, p)
%RHS Morris–Lecar + slow variable.
%
% z = [V; w; y]

V = z(1);
w = z(2);
y = z(3);

m_inf = 0.5*(1 + tanh((V - p.V1)/p.V2));
w_inf = 0.5*(1 + tanh((V - p.V3)/p.V4));
tau_w = cosh((V - p.V3)/(2*p.V4));

Ica = p.gca*m_inf*(V - p.Vca);
Kfac = p.gk*w + p.gkca*(y/(y + p.y0));
IK = Kfac*(V - p.Vk);
Ileak = p.gl*(V - p.Vl);

dV = (1/p.c) * (p.I - Ica - IK - Ileak);
dw = p.phi * tau_w * (w_inf - w);
dy = p.eps * (-p.mu*Ica - y);

dz = [dV; dw; dy];
end
