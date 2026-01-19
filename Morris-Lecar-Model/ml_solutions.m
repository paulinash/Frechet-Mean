function ml_solutions
% morris_lecar_book_form_3D_with_speed_control.m
% Morris–Lecar model (book form):
% - includes 2D plots
% - shows full 3D line
% - animated tracking point with real-valued speed control

clear; close all; clc;
set(groot, ...
    'DefaultTextInterpreter','latex', ...
    'DefaultAxesTickLabelInterpreter','latex', ...
    'DefaultLegendInterpreter','latex', ...
    'DefaultAxesFontName','Times', ...
    'DefaultTextFontName','Times', ...
    'DefaultAxesFontSize', 12, ...
    'DefaultTextFontSize', 12, ...
    'DefaultLegendFontSize', 11, ...
    'DefaultLineLineWidth', 1.2, ...
    'DefaultFigureColor','w');
%% --- Parameters (book) ---
Vk  = -84;    Vl = -60;    Vca = 120;
V1  = -1.2;   V2 = 18;
V3  = 12;     V4 = 17.4;
gk  = 8;      gl = 2;
c   = 20;
gca = 4;      gkca = 0.25;
y0  = 10;
phi = 0.23;  
mu = 0.2;
I = 45.5; % [44.8,45.5]
epsilon = 0.005;

% Pack parameters
p.Vk = Vk; p.Vl = Vl; p.Vca = Vca;
p.V1 = V1; p.V2 = V2; p.V3 = V3; p.V4 = V4;
p.gk = gk; p.gl = gl; p.c = c;
p.gca = gca; p.gkca = gkca; p.y0 = y0;
p.phi = phi; p.mu = mu; p.I = I; p.eps = epsilon;

%% --- Time span and ICs ---
tspan = [0 4000];
x0 = [-20; 0.0; 6];

opts = odeset('RelTol',1e-6,'AbsTol',1e-9,'MaxStep',1);

% Integrate
[t, X] = ode15s(@(t,x) ml_book_rhs(t,x,p), tspan, x0, opts);

% Extract
V = X(:,1);
w = X(:,2);
y = X(:,3);

%% --- 2D plots ---
figure('Units','normalized','Position',[0.05 0.05 0.9 0.8]);

subplot(2,2,1)
plot(t, V, 'k');
xlabel('t'); ylabel('V = x_1'); title('Membrane potential V(t)'); box on;

subplot(2,2,2)
plot(t, w);
xlabel('t'); ylabel('w '); title('Recovery variable w'); box on;

subplot(2,2,3)
plot(t, y);
xlabel('t'); ylabel('y (Ca)'); title('Slow variable y(t)'); box on;

subplot(2,2,4)
plot(V, w);
xlabel('V'); ylabel('w = x_2'); title('Phase: V vs w'); box on;



%% --- 3D line + moving point (real-valued speed control) ---
figure('Name','3D Trajectory with Tracker','NumberTitle','off',...
       'Units','normalized','Position',[0.1 0.1 0.65 0.65]);

% Plot full 3D trajectory
plot3(V, w, y, 'b-', 'LineWidth', 1.3);
hold on;

% Create tracking point
tracker = plot3(V(1), w(1), y(1), 'ro', ...
    'MarkerSize', 8, 'MarkerFaceColor', 'r');

xlabel('V = x_1');
ylabel('w = x_2');
zlabel('y (Ca)');
title('3D trajectory (line) with moving point');
grid on;
axis tight;
view([-30, 20]);
box on;

%% --- Animation speed control (real-valued) ---
speedFactor = 1;   % controls the duration of pause (larger value: slower animation)

% Normalize time differences
dt = diff(t);
dt = [dt; dt(end)];  % same length as t

for k = 1:length(t)
    % Move the point
    set(tracker, 'XData', V(k), 'YData', w(k), 'ZData', y(k));
    drawnow;

    % Pause proportional to actual simulation time and speed factor
    pause(speedFactor * dt(k) * 0.01);
end

hold off;

%% --- RHS function ---
function dxdt = ml_book_rhs(~, x, p)
    V = x(1);
    w = x(2);
    y = x(3);

    m_inf = 0.5*(1 + tanh((V - p.V1)/p.V2));
    w_inf = 0.5*(1 + tanh((V - p.V3)/p.V4));
    tau_w = cosh((V - p.V3)/(2*p.V4));

    Ica = p.gca * m_inf * (V - p.Vca);

    Kfactor = p.gk * w + p.gkca * (y/(y + p.y0));
    IK_total = Kfactor * (V - p.Vk);

    Ileak = p.gl * (V - p.Vl);

    dx1 = (1/p.c) * (p.I - Ica - IK_total - Ileak);
    dx2 = p.phi * tau_w * (w_inf - w);
    dy  = p.eps * (-p.mu * Ica - y);

    dxdt = [dx1; dx2; dy];
end
end