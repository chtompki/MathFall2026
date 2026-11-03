%% Problem 8.1.6: nullclines, equilibrium stability, and phase portraits
% System: x' = y - 2*x, y' = mu + x^2 - y.
% Run this script to create three figures matching problemSet4.tex.
% The script uses only base MATLAB (ode45 and graphics).
% Solid/dotted nullcline segments are classification guides at intersections;
% non-equilibrium points on a nullcline have no equilibrium stability.

blue = [0 0 1];
red = [1 0 0];
green = [0 0.45 0];

%% 1. Nullclines for mu < 1, mu = 1, and mu > 1
figure('Name','8.1.6: Nullclines and stability','Color','w', ...
    'Position',[100 100 1200 470]);
muCases = [0 1 2];
caseTitles = {'\mu < 1 (shown: \mu = 0)', '\mu = 1', ...
    '\mu > 1 (shown: \mu = 2)'};
for k = 1:3
    mu = muCases(k);
    ax = subplot(1,3,k); hold(ax,'on');
    if mu <= 1
        left = linspace(-0.65,1,250);
        right = linspace(1,3.05,250);
        plot(ax,left,2*left,'-','Color',blue,'LineWidth',1.6);
        plot(ax,right,2*right,':','Color',blue,'LineWidth',1.8);
        plot(ax,left,left.^2+mu,'-','Color',red,'LineWidth',1.6);
        plot(ax,right,right.^2+mu,':','Color',red,'LineWidth',1.8);
    else
        x = linspace(-0.65,3.05,500);
        plot(ax,x,2*x,'-','Color',blue,'LineWidth',0.8);
        plot(ax,x,x.^2+mu,'-','Color',red,'LineWidth',0.8);
    end
    markEquilibria(ax,mu);
    if mu < 1
        text(ax,0.15,0.5,'stable','BackgroundColor','w');
        text(ax,2.12,3.8,'saddle','BackgroundColor','w');
    elseif mu == 1
        text(ax,1.12,1.7,'saddle-node','BackgroundColor','w');
    else
        text(ax,0.5,1.4,'no equilibria','BackgroundColor','w');
    end
    setupAxes(ax,[-0.65 3.05 -0.65 6.6]);
    title(ax,caseTitles{k});
end
annotation('textbox',[0.07 0.005 0.9 0.075], ...
    'String',{'Blue: y = 2x; red: y = x^2 + \mu. Solid/dotted: attracting/saddle branch guides.', ...
    'At \mu = 1: one-sided attraction, not Lyapunov stable. At \mu > 1: no equilibrium to classify.'}, ...
    'EdgeColor','none','HorizontalAlignment','center');

%% 2. Bifurcation diagram: x = 1 +/- sqrt(1-mu)
figure('Name','8.1.6: Saddle-node bifurcation','Color','w');
ax = axes; hold(ax,'on');
mu = linspace(-2,1,700);
h1 = plot(ax,mu,1-sqrt(1-mu),'-','Color',blue,'LineWidth',1.8);
h2 = plot(ax,mu,1+sqrt(1-mu),'--','Color',red,'LineWidth',1.8);
plot(ax,1,1,'ko','MarkerFaceColor','k','MarkerSize',6);
text(ax,1.04,1.25,'saddle-node');
axis(ax,[-2 1.6 -0.9 2.9]);
xlabel(ax,'\mu'); ylabel(ax,'$1\pm\sqrt{1-\mu}$','Interpreter','latex');
title(ax,'Equilibrium branches: saddle-node at \mu = 1');
legend(ax,[h1 h2],{'Attracting: 1 - sqrt(1-mu)', ...
    'Saddle: 1 + sqrt(1-mu)'},'Location','southoutside');
grid(ax,'on'); box(ax,'on');

%% 3. Phase portraits, including saddle separatrices
figure('Name','8.1.6: Phase portraits','Color','w', ...
    'Position',[120 80 1050 850]);
muCases = [-1 0 1 1.5];
caseTitles = {'\mu = -1: stable spiral and saddle', ...
    '\mu = 0: stable node and saddle', ...
    '\mu = 1: saddle-node', '\mu = 3/2: no equilibria'};
limits = [-0.8 3.1 -1 6.5];
seeds = [-0.6 2.4; 0.3 4.8; 1.4 0.1; 2.7 1.5; 2.5 5.8; -0.5 0.3];
for k = 1:4
    mu = muCases(k);
    ax = subplot(2,2,k); hold(ax,'on');
    [X,Y] = meshgrid(-0.6:0.43:2.9,-0.7:0.62:6.2);
    U = Y-2*X; V = mu+X.^2-Y;
    % Normalize arrows with the same anisotropic scale used in the TikZ plots.
    scale = hypot(U,(0.61/1.12)*V);
    scale(scale < 1e-12) = 1;
    quiver(ax,X,Y,0.15*U./scale,0.15*V./scale,0, ...
        'Color',[0.7 0.7 0.7],'LineWidth',0.5);
    x = linspace(limits(1),limits(2),500);
    plot(ax,x,2*x,'-','Color',blue,'LineWidth',1.3);
    plot(ax,x,x.^2+mu,'--','Color',red,'LineWidth',1.3);
    starts = seeds;
    if mu == -1
        starts = [-0.1 -0.8; -0.65 -0.4; starts];
    end
    for j = 1:size(starts,1)
        Z = integrateOrbit(starts(j,:)',mu,1,limits);
        drawOrbit(ax,Z,mu,[0.15 0.15 0.15],0.8);
    end
    if mu < 1
        saddle = [1+sqrt(1-mu); 2+2*sqrt(1-mu)];
        lambdas = [(-3+sqrt(9+8*sqrt(1-mu)))/2, ...
                   (-3-sqrt(9+8*sqrt(1-mu)))/2];
        for lambda = lambdas
            v = [1; 2+lambda]; v = v/norm(v);
            for branchSign = [-1 1]
                z0 = saddle + branchSign*0.008*v;
                % Integrate stable branches backward to trace them outward.
                Z = integrateOrbit(z0,mu,sign(lambda),limits);
                drawOrbit(ax,Z,mu,green,1.4);
            end
        end
    end
    markEquilibria(ax,mu);
    setupAxes(ax,limits);
    title(ax,caseTitles{k});
end

%% Local functions
function dz = vectorField(~,z,mu)
    dz = [z(2)-2*z(1); mu+z(1)^2-z(2)];
end

function Z = integrateOrbit(z0,mu,timeSign,limits)
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.025, ...
        'Events',@(t,z) stopOrbit(t,z,mu,limits));
    [~,Z] = ode45(@(t,z) vectorField(t,z,mu),[0 timeSign*30],z0,opts);
end

function [value,isterminal,direction] = stopOrbit(~,z,mu,limits)
    value = [z(1)-limits(1); limits(2)-z(1); ...
             z(2)-limits(3); limits(4)-z(2); ...
             norm(vectorField(0,z,mu))-1e-5];
    isterminal = ones(5,1);
    direction = zeros(5,1);
end

function drawOrbit(ax,Z,mu,color,lineWidth)
    if size(Z,1) < 2, return; end
    plot(ax,Z(:,1),Z(:,2),'-','Color',color,'LineWidth',lineWidth);
    % Place the arrow by arc length. Its direction always uses forward time,
    % including when the stable separatrix was integrated backward.
    lengths = [0; cumsum(hypot(diff(Z(:,1)),(0.61/1.12)*diff(Z(:,2))))];
    if lengths(end) < 0.02, return; end
    idx = find(lengths >= 0.45*lengths(end),1);
    z = Z(idx,:)'; v = vectorField(0,z,mu);
    nv = hypot(v(1),(0.61/1.12)*v(2));
    if nv > 1e-10
        v = 0.12*v/nv;
        quiver(ax,z(1),z(2),v(1),v(2),0,'Color',color, ...
            'LineWidth',lineWidth,'MaxHeadSize',2);
    end
end

function markEquilibria(ax,mu)
    if mu < 1
        a = 1-sqrt(1-mu); b = 1+sqrt(1-mu);
        plot(ax,a,2*a,'ko','MarkerFaceColor','k','MarkerSize',6);
        plot(ax,b,2*b,'ko','MarkerFaceColor','w','MarkerSize',7);
    elseif mu == 1
        % Half-filled circle indicates the one-sided saddle-node behavior.
        plot(ax,1,2,'ko','MarkerFaceColor','w','MarkerSize',8);
        theta = linspace(pi/2,3*pi/2,60);
        patch(ax,1+0.035*cos(theta),2+0.064*sin(theta),'k', ...
            'EdgeColor','none');
    end
end

function setupAxes(ax,limits)
    axis(ax,limits);
    xlabel(ax,'x'); ylabel(ax,'y');
    set(ax,'XAxisLocation','origin','YAxisLocation','origin', ...
        'FontSize',10,'Layer','top');
    % Match the x/y unit proportions in the TikZ phase portraits.
    daspect(ax,[1 1.12/0.61 1]);
    box(ax,'off');
end
