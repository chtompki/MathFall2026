%% Problem 8.3.1: companion to the TikZ figures in problemSet4.tex
% Run this script in MATLAB. No toolboxes or external data files are needed.
% Numerical solutions use ode45; the TeX coordinates were independently
% computed using fourth-order Runge-Kutta integration.
% Solid curves indicate stable branches and dotted curves unstable branches
% where a stability classification applies. A nullcline is not itself a
% trajectory and its non-equilibrium points have no equilibrium stability.

% Manual: Exercise 8.3.1, printed pages 237--239.
a=1;bc=1+a;
fprintf('bc=%.8f, limiting period=%.8f\n',bc,2*pi/sqrt(a));
fprintf('Cubic coefficient at Hopf = %.8f (supercritical)\n',-(a+2)/8);
figure('Color','w','Name','8.3.1: Brusselator');
bs=[1.5 2.5];
for j=1:2
    b=bs(j);f=@(t,z) [1-(b+1)*z(1)+a*z(1)^2*z(2); ...
        b*z(1)-a*z(1)^2*z(2)];
    J=[b-1 a;-b -a];fprintf('b=%g, eigenvalues:\n',b);disp(eig(J));
    ax=subplot(1,2,j);hold(ax,'on');limits=[0 4 0 5];axis(ax,limits);
    directionField(ax,f,limits);x=linspace(.15,4,600);
    plot(ax,x,((b+1)*x-1)./(a*x.^2),'b-','LineWidth',1.3);
    plot(ax,x,b./(a*x),'r--','LineWidth',1.3);
    opts=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',.025);
    [t,Z]=ode45(f,[0 250],[1.1;b/a],opts);
    if b>bc,Z=Z(t>170,:);end
    orbitLine(ax,Z,f,[0 0 0],1.1);
    if b<bc,face='k';else,face='w';end
    plot(ax,1,b/a,'ko','MarkerFaceColor',face);
    finishAxes(ax,limits);title(ax,sprintf('a=%g, b=%g',a,b));
end
b=2.5;L=1/(b+1);H=b*(b+1)/a;R=2-L;
figure('Color','w','Name','8.3.1: manual trapping polygon');
plot([L R+H R L L],[0 0 H H 0],'k-','LineWidth',1.8);hold on;
plot(1,b/a,'ko','MarkerFaceColor','w');xlabel('x');ylabel('y');grid on;
title('Inward boundary: x=L, y=0, y=H, x+y=R+H');

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
