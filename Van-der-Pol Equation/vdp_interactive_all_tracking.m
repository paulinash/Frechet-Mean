function vdp_interactive_all_full_params()
% Interactive Van der Pol with all sliders (a,b,c,d,e,f,ε)
% Trajectory and moving point adjust to all parameter changes.

clear; close all; clc;

%% --- Settings
Tfinal = 400; dt = 0.01; M = 500; useArcLength = false;
x0 = [0.5;0]; 

%% --- Figure
fig = figure('Name','Van der Pol Interactive (All Params)','NumberTitle','off',...
    'Units','normalized','Position',[0.05 0.05 0.9 0.9]);
ax = axes('Parent',fig); hold(ax,'on'); grid on;
xlabel(ax,'x'); ylabel(ax,'y'); axis(ax,[-3 3 -4 4]);
ax.Position = [0.05 0.05 0.65 0.9];

[xv, yv] = meshgrid(linspace(-3,3,20), linspace(-4,4,20));

%% --- Initial parameters
a=0; b=1; c=1; d=1; e=0; f=0; epsilon=0.1;

%% --- Plots
h_traj = plot(ax,NaN,NaN,'r','LineWidth',2);
h_dot  = plot(ax,NaN,NaN,'ro','MarkerFaceColor','y','MarkerSize',10);
h_xnull = plot(ax,NaN,NaN,'k--','LineWidth',1.5);
h_yline = xline(ax,a,'b--','LineWidth',1.5);

u = b*yv - c*xv.^3/3 + d*xv + e;
v = epsilon*(a - xv) + f;
m = sqrt(u.^2+v.^2); u=u./m; v=v./m;
h_quiver = quiver(ax,xv,yv,u,v,0.5,'Color',[0.7 0.7 0.7]);

%% --- Sliders
ypos = struct('a',0.38,'b',0.32,'c',0.26,'d',0.20,'e',0.14,'f',0.08,'eps',0.02);
param_range = [-1,4]; epsilon_range = [0.01,2];

[slider_a,label_a] = make_slider(ypos.a,param_range,a,'a');
[slider_b,label_b] = make_slider(ypos.b,param_range,b,'b');
[slider_c,label_c] = make_slider(ypos.c,param_range,c,'c');
[slider_d,label_d] = make_slider(ypos.d,param_range,d,'d');
[slider_e,label_e] = make_slider(ypos.e,param_range,e,'e');
[slider_f,label_f] = make_slider(ypos.f,param_range,f,'f');
slider_eps = uicontrol('Parent',fig,'Style','slider','Units','normalized',...
    'Position',[0.72 ypos.eps 0.20 0.04],'Min',epsilon_range(1),...
    'Max',epsilon_range(2),'Value',epsilon);
label_eps = uicontrol('Parent',fig,'Style','text','Units','normalized',...
    'Position',[0.94 ypos.eps 0.05 0.04],'String',sprintf('ε=%.2f',epsilon));

%% --- Trajectory storage
trajData = []; trajIndex = 1;

%% --- Timer
animTimer = timer('ExecutionMode','fixedRate','Period',0.03,'TimerFcn',@animateDot);

%% --- Update function
function updatePlot(~,~)
    % Read sliders
    a = slider_a.Value; b = slider_b.Value;
    c = slider_c.Value; d = slider_d.Value;
    e = slider_e.Value; f = slider_f.Value;
    epsilon = slider_eps.Value;
    
    % Update labels
    label_a.String = sprintf('a=%.2f',a);
    label_b.String = sprintf('b=%.2f',b);
    label_c.String = sprintf('c=%.2f',c);
    label_d.String = sprintf('d=%.2f',d);
    label_e.String = sprintf('e=%.2f',e);
    label_f.String = sprintf('f=%.2f',f);
    label_eps.String = sprintf('ε=%.2f',epsilon);
    
    % Vector field
    u = b*yv - c*xv.^3/3 + d*xv + e;
    v = epsilon*(a - xv) + f;
    m = sqrt(u.^2+v.^2); u=u./m; v=v./m;
    set(h_quiver,'UData',u,'VData',v);
    
    % Nullclines
    xgrid = linspace(-3,3,400);
    set(h_xnull,'XData',xgrid,'YData',(c*xgrid.^3/3 - d*xgrid - e)/b);
    h_yline.Value = a + f/epsilon;
    
    % Equilibrium (optional)
    x_eq = a + f/epsilon;
    y_eq = (c*x_eq^3/3 - d*x_eq - e)/b;
    
    % --- Solve trajectory (full param-dependent ODE)
    vdp = @(t,z)[ b*z(2) - c*z(1)^3/3 + d*z(1) + e; epsilon*(a - z(1)) + f ];
    t = 0:dt:Tfinal;
    [~,Z] = ode45(vdp,t,x0);
    
    % Last half + resample
    Z_ss = Z(t >= Tfinal/2,:);

    if useArcLength
        Z_res = resampleArcLength(Z_ss,M);
    else
        Z_res = interp1(linspace(0,1,size(Z_ss,1)),Z_ss,linspace(0,1,M));
    end
    
    trajData = Z_res; trajIndex = 1;
    
    % Update plots
    set(h_traj,'XData',Z_res(:,1),'YData',Z_res(:,2));
    set(h_dot,'XData',Z_res(1,1),'YData',Z_res(1,2));
    
    stop(animTimer); start(animTimer);
end

%% --- Animate dot
function animateDot(~,~)
    if isempty(trajData), return; end
    if trajIndex>size(trajData,1), trajIndex=1; end
    set(h_dot,'XData',trajData(trajIndex,1),'YData',trajData(trajIndex,2));
    trajIndex = trajIndex+1;
end

%% --- Resampling
function curveOut = resampleArcLength(curveIn,M)
    diffs = diff(curveIn,1,1);
    segLen = sqrt(sum(diffs.^2,2));
    cumLen = [0;cumsum(segLen)];
    t_new = linspace(0,cumLen(end),M);
    x_new = interp1(cumLen,curveIn(:,1),t_new);
    y_new = interp1(cumLen,curveIn(:,2),t_new);
    curveOut = [x_new(:),y_new(:)];
end

%% --- Slider helper
function [slider,label] = make_slider(ypos,range,init_val,name)
    slider = uicontrol('Parent',fig,'Style','slider','Units','normalized',...
        'Position',[0.72 ypos 0.20 0.04],'Min',range(1),'Max',range(2),'Value',init_val);
    label = uicontrol('Parent',fig,'Style','text','Units','normalized',...
        'Position',[0.94 ypos 0.05 0.04],'String',sprintf('%s=%.2f',name,init_val));
end

%% --- Connect callbacks
sliders = [slider_a,slider_b,slider_c,slider_d,slider_e,slider_f,slider_eps];
set(sliders,'Callback',@updatePlot);

%% --- Initial plot
updatePlot();

end