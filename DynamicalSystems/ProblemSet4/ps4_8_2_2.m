%% Problem 8.2.2: companion to the TikZ figures in problemSet4.tex
% Run this script in MATLAB. No toolboxes or external data files are needed.
% Numerical solutions use ode45; the TeX coordinates were independently
% computed using fourth-order Runge-Kutta integration.
% Solid curves indicate stable branches and dotted curves unstable branches
% where a stability classification applies. A nullcline is not itself a
% trajectory and its non-equilibrium points have no equilibrium stability.

J=@(mu) [mu -1;1 mu];
for mu=[-.1 0 .1]
    fprintf('mu=%g, eigenvalues:\n',mu);disp(eig(J(mu)));
end
mu=linspace(-.1,.1,300);
figure('Color','w','Name','8.2.2: eigenvalue crossing');
subplot(1,2,1);plot(mu,mu,'b-','LineWidth',1.8);hold on;
plot(mu,zeros(size(mu)),'k:');xlabel('mu');ylabel('Re(lambda)');grid on;
title('Both real parts equal mu');
subplot(1,2,2);plot(mu,ones(size(mu)),'b-',mu,-ones(size(mu)),'r-','LineWidth',1.8);
xlabel('mu');ylabel('Im(lambda)');grid on;title('Imaginary parts are +1 and -1');
% At mu=0, a linear center does not imply a nonlinear center.
