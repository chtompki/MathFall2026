%% Problem 8.5.4: companion to the TikZ figures in problemSet4.tex
% Run this script in MATLAB. No toolboxes or external data files are needed.
% Numerical solutions use ode45; the TeX coordinates were independently
% computed using fourth-order Runge-Kutta integration.
% Solid curves indicate stable branches and dotted curves unstable branches
% where a stability classification applies. A nullcline is not itself a
% trajectory and its non-equilibrium points have no equilibrium stability.

r=1;alpha=.5;T=2*pi;
low=r/(4*(1+alpha));high=r/4;
hc=fzero(@(h) hillTrace(h,r,alpha)-2,[low high]);
fprintf('Lower sufficient bound = %.9f\n',low);
fprintf('Computed critical h = %.9f; upper bound = %.9f\n',hc,high);
figure('Color','w','Name','8.5.4: seasonal harvesting');
hs=[.15 .23 .247];
for j=1:3
    h=hs(j);f=@(t,x) r*x*(1-x)-h*(1+alpha*sin(t));
    opts=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',.04);
    ax=subplot(1,3,j);hold(ax,'on');
    if h<hc
        [~,X]=ode45(f,[0 100*T],.9,opts);upper=X(end);
        [~,X]=ode45(f,[0 -100*T],.15,opts);lower=X(end);
        [t,U]=ode45(f,[0 2*T],upper,opts);
        [~,L]=ode45(f,t,lower,opts);
        plot(ax,t,U,'b-','LineWidth',1.6);
        plot(ax,t,L,'r:','LineWidth',1.8);
        multiplierUpper=exp(trapz(t,r*(1-2*U))/2);
        multiplierLower=exp(trapz(t,r*(1-2*L))/2);
        fprintf('h=%.3f: cycle x(0)=[%.8f,%.8f], multipliers=[%.6f,%.6f]\n', ...
            h,lower,upper,multiplierLower,multiplierUpper);
        xlim(ax,[0 2*T]);
    else
        opts=odeset(opts,'Events',@extinction);
        for x0=[.9 .6 .4]
            [t,X]=ode45(f,[0 8*T],x0,opts);
            plot(ax,t,X,'k-','LineWidth',1.1);
        end
        xlim(ax,[0 8*T]);
    end
    ylim(ax,[0 1]);xlabel(ax,'t');ylabel(ax,'x');grid(ax,'on');
    title(ax,sprintf('r=1, alpha=0.5, h=%.3f',h));
end
figure('Color','w','Name','8.5.4: Floquet threshold');
h=linspace(low,high,100);tr=zeros(size(h));
for j=1:numel(h),tr(j)=hillTrace(h(j),r,alpha);end
plot(h,tr,'b-','LineWidth',1.5);hold on;plot([low high],[2 2],'k:');
plot(hc,2,'ko','MarkerFaceColor','k');xlabel('h');ylabel('trace M(h)');grid on;
title('First trace=2 crossing is the fold of periodic populations');

function tr=hillTrace(h,r,alpha)
    f=@(t,z) [z(2);(r*r/4-r*h*(1+alpha*sin(t)))*z(1); ...
        z(4);(r*r/4-r*h*(1+alpha*sin(t)))*z(3)];
    opts=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',.025);
    [~,Z]=ode45(f,[0 2*pi],[1;0;0;1],opts);
    M=[Z(end,1) Z(end,3);Z(end,2) Z(end,4)];
    tr=trace(M);
end
function [value,isterminal,direction]=extinction(~,x)
    value=x;isterminal=1;direction=-1;
end
