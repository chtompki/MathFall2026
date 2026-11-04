%% Problem 8.2.10: companion to the TikZ figures in problemSet4.tex
% Run this script in MATLAB. No toolboxes or external data files are needed.
% Numerical solutions use ode45; the TeX coordinates were independently
% computed using fourth-order Runge-Kutta integration.
% Solid curves indicate stable branches and dotted curves unstable branches
% where a stability classification applies. A nullcline is not itself a
% trajectory and its non-equilibrium points have no equilibrium stability.

A=2;B=3;q=10;d=B-A;
if d<=0,error('Choose B>A for a physical equilibrium.');end
c=d/(1+q*d^2);k=A*(1-q*d^2)/(d*(1+q*d^2));
J=[-1-k -c;-k -c];
fprintf('Equilibrium: (%.8f, %.8f)\n',d,A*(1+q*d^2)/d);
fprintf('Trace=%.8f, determinant=%.8f\n',trace(J),det(J));disp(eig(J));
if q*d^2>1
    AH=d*(1+q*d^2+d)/(q*d^2-1);
    fprintf('Trace-zero threshold at fixed d,q: AH=%.8f, BH=%.8f\n',AH,AH+d);
end
f=@(t,z) [B-z(1)-z(1)*z(2)/(1+q*z(1)^2); ...
    A-z(1)*z(2)/(1+q*z(1)^2)];
opts=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',.04);
[t,Z]=ode45(f,[0 700],[1.05;22],opts);
figure('Color','w','Name','8.2.10: bacterial respiration');
ax=subplot(1,2,1);hold(ax,'on');limits=[0 2.7 12 29];axis(ax,limits);
directionField(ax,f,limits);x=linspace(.1,2.7,600);
plot(ax,x,(B-x).*(1+q*x.^2)./x,'b-','LineWidth',1.4);
plot(ax,x,A*(1+q*x.^2)./x,'r--','LineWidth',1.4);
plot(ax,Z(t>550,1),Z(t>550,2),'k-','LineWidth',1.3);
plot(ax,d,A*(1+q*d^2)/d,'ko','MarkerFaceColor','w');
finishAxes(ax,limits);title(ax,'Nullclines and stable cycle');
ax=subplot(1,2,2);hold(ax,'on');
plot(ax,t,Z(:,1),'b-',t,Z(:,2),'r-','LineWidth',1.2);
xlim(ax,[500 700]);xlabel(ax,'t');ylabel(ax,'Concentration');
legend(ax,{'x','y'});grid(ax,'on');title(ax,'Late-time oscillations');
% Show the invariant polygon from the proof separately.
M=B+1;epsilon=min(d/2,M/2);
x=linspace(epsilon,M,1000);Y=max([epsilon A*(1+q*x.^2)./x])+1;
figure('Color','w','Name','8.2.10: trapping polygon');
plot([0 M M epsilon 0 0],[0 0 Y Y Y-epsilon 0],'k-','LineWidth',1.6);
hold on;plot(d,A*(1+q*d^2)/d,'ko');xlabel('x');ylabel('y');grid on;
title('Compact trapping polygon (drawn on its full scale)');

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
