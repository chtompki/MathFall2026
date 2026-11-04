%% Problem 8.4.4: companion to the TikZ figures in problemSet4.tex
% Run this script in MATLAB. No toolboxes or external data files are needed.
% Numerical solutions use ode45; the TeX coordinates were independently
% computed using fourth-order Runge-Kutta integration.
% Solid curves indicate stable branches and dotted curves unstable branches
% where a stability classification applies. A nullcline is not itself a
% trajectory and its non-equilibrium points have no equilibrium stability.

muHom=fzero(@(mu) sepVelocity(mu,false)-sepVelocity(mu,true),[3.5 4]);
fprintf('Homoclinic threshold, epsilon=1e-6: %.9f\n',muHom);
% Check sensitivity to the saddle starting offset.
muRefined=fzero(@(mu) sepVelocity(mu,false,1e-7)- ...
    sepVelocity(mu,true,1e-7),[3.5 4]);
fprintf('With epsilon=1e-7: %.9f\n',muRefined);
figure('Color','w','Name','8.4.4: phase-locked loop');
mus=[.8 2 3.7 3.9];limits=[-pi pi -7 7];
for j=1:4
    mu=mus(j);f=@(t,z) pllField(t,z,mu);
    ax=subplot(2,2,j);hold(ax,'on');axis(ax,limits);directionField(ax,f,limits);
    opts=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',.025);
    [t,Z]=ode45(f,[0 350],[.1;0],opts);
    if mu>1,Z=Z(t>200,:);end
    theta=mod(Z(:,1)+pi,2*pi)-pi;
    cut=[false;abs(diff(theta))>pi];theta(cut)=NaN;
    plot(ax,theta,Z(:,2),'k-','LineWidth',1.1);
    if mu<1,face='k';else,face='w';end
    plot(ax,0,0,'ko','MarkerFaceColor',face);
    plot(ax,[-pi pi],[0 0],'ko','MarkerFaceColor','w');
    finishAxes(ax,limits);xlabel(ax,'theta');ylabel(ax,'v');
    title(ax,sprintf('mu = %g',mu));
end
% On the plane the connection is heteroclinic between adjacent saddle copies;
% on the phase cylinder they are the same saddle: a homoclinic connection.
figure('Color','w','Name','8.4.4: critical saddle connection');hold on;
mu=muRefined;f=@(t,z) pllField(t,z,mu);ep=1e-7;
lu=(-(1+mu)+sqrt((1+mu)^2+4))/2;
opts=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',.01, ...
    'Events',@(t,z) rightSaddle(t,z));
[~,Z]=ode45(f,[0 350],[-pi+ep;lu*ep],opts);
plot(Z(:,1),Z(:,2),'g-','LineWidth',1.6);
plot(-Z(:,1),-Z(:,2),'g-','LineWidth',1.6);
plot([-pi pi],[0 0],'ko');xlabel('theta');ylabel('v');grid on;
title(sprintf('Near saddle connection, mu = %.6f',mu));

function dz=pllField(~,z,mu)
    dz=[z(2);-(1-mu*cos(z(1)))*z(2)-sin(z(1))];
end
function v=sepVelocity(mu,stable,epsilon)
    if nargin<3,epsilon=1e-6;end
    a=1+mu;
    if stable
        lambda=(-a-sqrt(a*a+4))/2;z0=[pi-epsilon;-lambda*epsilon];T=-300;
    else
        lambda=(-a+sqrt(a*a+4))/2;z0=[-pi+epsilon;lambda*epsilon];T=300;
    end
    opts=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',.02,'Events',@sectionZero);
    [~,~,~,Ze]=ode45(@(t,z) pllField(t,z,mu),[0 T],z0,opts);
    if isempty(Ze),error('Separatrix did not reach theta=0.');end
    v=Ze(end,2);
end
function [value,isterminal,direction]=sectionZero(~,z)
    value=z(1);isterminal=1;direction=0;
end
function [value,isterminal,direction]=rightSaddle(~,z)
    value=z(1)-(pi-1e-4);isterminal=1;direction=1;
end

%% Local plotting and integration helpers
function Z = trajectory(f,z0,tspan,limits)
    opts=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',0.04, ...
        'Events',@(t,z) windowExit(t,z,limits));
    [~,Z]=ode45(f,tspan,z0,opts);
end
function [value,isterminal,direction]=windowExit(~,z,limits)
    value=[z(1)-limits(1);limits(2)-z(1);z(2)-limits(3);limits(4)-z(2)];
    isterminal=ones(4,1); direction=zeros(4,1);
end
function directionField(ax,f,limits)
    [X,Y]=meshgrid(linspace(limits(1),limits(2),15), ...
        linspace(limits(3),limits(4),13));
    U=zeros(size(X));V=U;
    for j=1:numel(X)
        w=f(0,[X(j);Y(j)]); U(j)=w(1);V(j)=w(2);
    end
    dx=(limits(2)-limits(1))/15;dy=(limits(4)-limits(3))/13;
    n=hypot(U/dx,V/dy);n(n<1e-12)=1;
    quiver(ax,X,Y,0.35*U./n,0.35*V./n,0,'Color',[.75 .75 .75]);
end
function orbitLine(ax,Z,f,color,width)
    plot(ax,Z(:,1),Z(:,2),'-','Color',color,'LineWidth',width);
    if size(Z,1)<3,return;end
    j=max(1,round(size(Z,1)*0.45));w=f(0,Z(j,:)');
    xl=xlim(ax);yl=ylim(ax);n=hypot(w(1)/diff(xl),w(2)/diff(yl));
    if n>1e-12
        w=0.025*w/n;
        quiver(ax,Z(j,1),Z(j,2),w(1),w(2),0,'Color',color,'MaxHeadSize',2);
    end
end
function finishAxes(ax,limits)
    axis(ax,limits);xlabel(ax,'x');ylabel(ax,'y');grid(ax,'on');box(ax,'on');
end
