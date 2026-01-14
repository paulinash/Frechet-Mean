function rm_slider_explorer
% Interactive slider explorer for Rosenzweig–MacArthur model
% Now includes moving tracker dot along trajectory

clear; close all; clc;

%% ========== DEFAULT PARAMETERS ==========
alpha0 = 0.1;
beta0  = 0.5;
gamma0 = 1.5;

Tend = 300;
dt   = 0.02;
tspan = 0:dt:Tend;

%% ========== FIGURE AND AXES ==========
fig = figure('Color','w','Position',[200 100 1100 650]);

% Main axes
ax = axes('Parent',fig,'Position',[0.08 0.12 0.65 0.80]);
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
xlabel(ax,'x (prey)');
ylabel(ax,'y (predator)');
title(ax,'Rosenzweig–MacArthur Phase Portrait');

% Fixed limits
xlim(ax,[0 6]);
ylim(ax,[0 6]);

%% ========== UI SLIDERS (RIGHT PANEL) ==========
sAlpha = uicontrol('Style','slider','Min',0.01,'Max',0.95,'Value',alpha0,...
    'Units','normalized','Position',[0.78 0.65 0.18 0.04],'Callback',@updatePlot);

sBeta = uicontrol('Style','slider','Min',0.1,'Max',5,'Value',beta0,...
    'Units','normalized','Position',[0.78 0.55 0.18 0.04],'Callback',@updatePlot);

sGamma = uicontrol('Style','slider','Min',0.2,'Max',10,'Value',gamma0,...
    'Units','normalized','Position',[0.78 0.75 0.18 0.04],'Callback',@updatePlot);

% Labels
uicontrol('Style','text','String','\alpha','Units','normalized',...
    'Position',[0.72 0.65 0.05 0.04],'FontSize',11);
uicontrol('Style','text','String','\beta','Units','normalized',...
    'Position',[0.72 0.55 0.05 0.04],'FontSize',11);
uicontrol('Style','text','String','\gamma','Units','normalized',...
    'Position',[0.72 0.75 0.05 0.04],'FontSize',11);

% Value readouts
txtAlpha = uicontrol('Style','text','Units','normalized',...
    'Position',[0.78 0.60 0.18 0.035],'FontSize',11);
txtBeta  = uicontrol('Style','text','Units','normalized',...
    'Position',[0.78 0.50 0.18 0.035],'FontSize',11);
txtGamma = uicontrol('Style','text','Units','normalized',...
    'Position',[0.78 0.70 0.18 0.035],'FontSize',11);

% Hopf info
txtHopf = uicontrol('Style','text','Units','normalized',...
    'Position',[0.75 0.40 0.22 0.07],'FontSize',11,...
    'HorizontalAlignment','left');

%% ========== GRAPHICS OBJECTS ==========
trajLine   = plot(ax,NaN,NaN,'b','LineWidth',1.5);
eqPoint    = plot(ax,NaN,NaN,'ko','MarkerFaceColor','k');
trackerDot = plot(ax,NaN,NaN,'ro','MarkerFaceColor','r','MarkerSize',8);

%% ========== TIMER FOR ANIMATION ==========
animTimer = timer('ExecutionMode','fixedRate','Period',0.03,'TimerFcn',@animate);

% Shared animation data
trajX = []; 
trajY = [];
idxTracker = 1;

%% ========== INITIAL UPDATE ==========
updatePlot();

%% ========== CALLBACK: UPDATE PLOT ==========
function updatePlot(~,~)
    % Stop animation while updating
    stop(animTimer);

    % Read parameters
    alpha = get(sAlpha,'Value');
    beta  = get(sBeta,'Value');
    gamma = get(sGamma,'Value');

    % Update text
    txtAlpha.String = sprintf('\\alpha = %.3f',alpha);
    txtBeta.String  = sprintf('\\beta = %.3f',beta);
    txtGamma.String = sprintf('\\gamma = %.3f',gamma);

    % Hopf threshold
    if alpha < 1
        gammaH = (1+alpha)/(1-alpha);
    else
        gammaH = NaN;
    end

    if ~isnan(gammaH)
        if gamma > gammaH
            txtHopf.String = sprintf('gamma_H = %.3f\nLimit cycle likely',gammaH);
            txtHopf.ForegroundColor = [0 0.4 0];
        else
            txtHopf.String = sprintf('gamma_H = %.3f\nStable equilibrium',gammaH);
            txtHopf.ForegroundColor = [0.6 0 0];
        end
    else
        txtHopf.String = 'No coexistence equilibrium';
        txtHopf.ForegroundColor = [0.6 0 0];
    end

    % Define system
    rm = @(t,z) [
        z(1)*(1 - z(1)/gamma) - z(1)*z(2)/(1+z(1));
        beta*(z(1)/(1+z(1)) - alpha)*z(2)
    ];

    % Solve ODE
    z0 = [0.5*gamma; 0.5];
    opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [~,z] = ode45(rm,tspan,z0,opts);

    trajX = z(:,1);
    trajY = z(:,2);
    idxTracker = 1;

    % Update graphics
    set(trajLine,'XData',trajX,'YData',trajY);

    % Plot equilibrium
    if alpha < 1
        xstar = alpha/(1-alpha);
        ystar = (1+xstar)*(1 - xstar/gamma);
        if ystar > 0
            set(eqPoint,'XData',xstar,'YData',ystar);
        else
            set(eqPoint,'XData',NaN,'YData',NaN);
        end
    else
        set(eqPoint,'XData',NaN,'YData',NaN);
    end

    % Reset tracker
    set(trackerDot,'XData',trajX(1),'YData',trajY(1));

    % Restart animation
    start(animTimer);
end

%% ========== TIMER CALLBACK ==========
function animate(~,~)
    if isempty(trajX)
        return
    end

    idxTracker = idxTracker + 2;
    if idxTracker > numel(trajX)
        idxTracker = 1;
    end

    set(trackerDot,'XData',trajX(idxTracker),'YData',trajY(idxTracker));
end

end
