%% Problem 8.1.15: companion to the TikZ figures in problemSet4.tex
% Run this script in MATLAB. No toolboxes or external data files are needed.
% Numerical solutions use ode45; the TeX coordinates were independently
% computed using fourth-order Runge-Kutta integration.
% Solid curves indicate stable branches and dotted curves unstable branches
% where a stability classification applies. A nullcline is not itself a
% trajectory and its non-equilibrium points have no equilibrium stability.

% The manual's corrected model is used: the final loss in nB' is
% -(p+nA)*nB, not -(p+nA)*nAB. Manual: printed pages 224--226.
pc=1-sqrt(3)/2;
f=@(p) @(t,z) [(p+z(1))*(1-p-z(1)-z(2))-z(1)*z(2); ...
    z(2)*(1-2*p-2*z(1)-z(2))];
figure('Color','w','Name','8.1.15: opinion dynamics');
ax=subplot(1,2,1);hold(ax,'on');
ps=[.10 .13 .14 .15];colors=lines(numel(ps));
for j=1:numel(ps)
    p=ps(j);opts=odeset('RelTol',1e-10,'AbsTol',1e-12);
    [t,Z]=ode45(f(p),[0 160],[0;1-p],opts);
    plot(ax,t,Z(:,2),'Color',colors(j,:),'LineWidth',1.4, ...
        'DisplayName',sprintf('p = %.2f',p));
    fprintf('p=%.2f: nA(160)=%.7f, nB(160)=%.7f\n',p,Z(end,1),Z(end,2));
end
xlabel(ax,'t');ylabel(ax,'n_B');legend(ax,'show');grid(ax,'on');
title(ax,'Initial state: n_A=0, n_B=1-p');
ax=subplot(1,2,2);hold(ax,'on');
p=linspace(0,pc,500);D=sqrt(max(0,4*p.^2-8*p+1));
uLower=(1-4*p-D)/6;uUpper=(1-4*p+D)/6;
plot(ax,p,1-2*p-2*uLower,'b-','LineWidth',1.8);
plot(ax,p,1-2*p-2*uUpper,'r:','LineWidth',1.8);
plot(ax,[0 .2],[0 0],'b-','LineWidth',1.8);
plot(ax,pc,1/sqrt(3),'ko','MarkerFaceColor','k');
xlabel(ax,'p');ylabel(ax,'n_B');axis(ax,[0 .2 0 1]);grid(ax,'on');
title(ax,'Saddle-node at p_c = 1 - sqrt(3)/2');
legend(ax,{'Stable B-dominated','Saddle','Stable A-consensus'},'Location','best');
fprintf('pc = %.10f\n',pc);
% Direct differentiation at consensus gives [-1,p-2;0,-1].
% Both consensus eigenvalues are -1; this check corrects a manual misprint.
