%% Problem 8.2.4: companion to the TikZ figures in problemSet4.tex
% Run this script in MATLAB. No toolboxes or external data files are needed.
% Numerical solutions use ode45; the TeX coordinates were independently
% computed using fourth-order Runge-Kutta integration.
% Solid curves indicate stable branches and dotted curves unstable branches
% where a stability classification applies. A nullcline is not itself a
% trajectory and its non-equilibrium points have no equilibrium stability.

mus=[-.01 -.0025 -.000625];radii=zeros(size(mus));
fprintf('mu          sqrt(-8*mu)    shooting section radius\n');
for j=1:numel(mus)
    radii(j)=shootCycle(mus(j));
    fprintf('% .7f    %.8f       %.8f\n',mus(j),sqrt(-8*mus(j)),radii(j));
end
figure('Color','w','Name','8.2.4: small unstable cycle');
mu=linspace(-.01,0,300);
plot(mu,sqrt(-8*mu),'b-','LineWidth',1.6);hold on;
plot(mus,radii,'ro-','LineWidth',1.3);grid on;
xlabel('mu');ylabel('Section radius R');
legend('Leading prediction sqrt(-8 mu)','Numerical shooting','Location','best');
title('Radius at upward crossing of y=0, x>0');
hopfPortraits();

function dz=hopfField(~,z,mu)
    x=z(1);y=z(2);dz=[-y+mu*x+x*y^2;x+mu*y-x^2];
end
function R=shootCycle(mu)
    estimate=sqrt(-8*mu);
    R=fzero(@(r) returnRadius(r,mu)-r,[.5*estimate 1.5*estimate]);
end
function R=returnRadius(r,mu)
    f=@(t,z) hopfField(t,z,mu);
    opts=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',.02, ...
        'Events',@returnEvent);
    [~,~,~,Ze,Ie]=ode45(f,[0 16],[r;0],opts);
    if isempty(Ie),error('No return to the section was found.');end
    if Ie(end)==2,R=2;else,R=Ze(end,1);end
end
function [value,isterminal,direction]=returnEvent(t,z)
    if t<1,v=1;else,v=z(2);end
    value=[v;2-norm(z)];isterminal=[1;1];direction=[1;-1];
end
function hopfPortraits()
    mus=[-.01 -.0025 0 .01];limits=[-.5 .5 -.5 .5];
    figure('Color','w','Name','Local subcritical Hopf portraits');
    for j=1:4
        mu=mus(j);ax=subplot(2,2,j);hold(ax,'on');axis(ax,limits);
        f=@(t,z) hopfField(t,z,mu);directionField(ax,f,limits);
        for r=[.07 .20 .38]
            Z=trajectory(f,[r;0],[0 160],limits);
            orbitLine(ax,Z,f,[.15 .15 .15],.8);
        end
        if mu<0
            R=shootCycle(mu);opts=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',.01);
            [~,Z]=ode45(f,[0 6.6],[R;0],opts);
            plot(ax,Z(:,1),Z(:,2),'r:','LineWidth',1.8);
            plot(ax,0,0,'ko','MarkerFaceColor','k');
        else
            plot(ax,0,0,'ko','MarkerFaceColor','w');
        end
        finishAxes(ax,limits);title(ax,sprintf('mu = %g',mu));
    end
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
