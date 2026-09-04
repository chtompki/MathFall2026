% Exercise 2.8.6: x' = x + exp(-x), x(0) = 0
clear; close all; clc;

f = @(t,x) x + exp(-x);

% High-accuracy reference solution.
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
[tRef,xRef] = ode45(f,[0 1],0,opts);
xExact = xRef(end);

% Explicit Euler approximations at t=1.
Nvals = [10 20 50 100 200 500 1000 2000 2500 5000 10000];
hvals = 1./Nvals;
xEuler = zeros(size(Nvals));
for j = 1:numel(Nvals)
    N = Nvals(j);
    h = 1/N;
    x = 0;
    for n = 1:N
        x = x + h*(x + exp(-x));
    end
    xEuler(j) = x;
end
errors = abs(xEuler-xExact);

fprintf('High-accuracy reference: x(1) = %.12f\n',xExact);
fprintf('       N            h          Euler x(1)       error\n');
fprintf('%8d   %11.4g   %15.10f   %11.4g\n', ...
    [Nvals; hvals; xEuler; errors]);

% A fine Euler trajectory for the solution plot.
Nplot = 5000;
hplot = 1/Nplot;
tE = (0:Nplot)*hplot;
xE = zeros(size(tE));
for n = 1:Nplot
    xE(n+1) = xE(n) + hplot*(xE(n) + exp(-xE(n)));
end

fig = figure('Color','w','Position',[100 100 1100 450]);

% Solution together with a slope field.
subplot(1,2,1);
tGrid = linspace(0,1,17);
xGrid = linspace(0,1.5,17);
[T,X] = meshgrid(tGrid,xGrid);
S = X + exp(-X);
L = sqrt(1+S.^2);
quiver(T,X,1./L,S./L,0.42,'Color',[0.68 0.68 0.68], ...
    'LineWidth',0.7,'MaxHeadSize',0.25);
hold on;
plot(tRef,xRef,'b-','LineWidth',2.2);
plot(tE,xE,'r--','LineWidth',1.3);
plot(1,xExact,'ko','MarkerFaceColor','k','MarkerSize',5);
xlim([0 1]); ylim([0 1.5]); grid on; box on;
xlabel('$t$','Interpreter','latex');
ylabel('$x(t)$','Interpreter','latex');
title('$\dot{x}=x+e^{-x},\quad x(0)=0$','Interpreter','latex');
legend('slope field','ode45','Euler, $h=2\times10^{-4}$','$x(1)$', ...
    'Interpreter','latex','Location','northwest');

% First-order convergence of Euler's method.
subplot(1,2,2);
loglog(hvals,errors,'o-','LineWidth',1.7,'MarkerFaceColor',[0.2 0.45 0.8]);
hold on;
loglog(hvals,errors(end)*(hvals/hvals(end)),'k--','LineWidth',1.2);
yline(5e-4,'r:','LineWidth',1.2);
grid on; box on;
xlabel('stepsize $h$','Interpreter','latex');
ylabel('$|x_{\rm Euler}(1)-x_{\rm ref}(1)|$','Interpreter','latex');
title('Euler global error','Interpreter','latex');
legend('computed error','$O(h)$ reference','$5\times10^{-4}$ tolerance', ...
    'Interpreter','latex','Location','southeast');

if exist('exportgraphics','file') == 2
    exportgraphics(fig,'exercise_2_8_6_matlab.png','Resolution',300);
    exportgraphics(fig,'exercise_2_8_6_matlab.pdf','ContentType','vector');
else
    print(fig,'exercise_2_8_6_matlab.png','-dpng','-r300');
    print(fig,'exercise_2_8_6_matlab.pdf','-dpdf','-bestfit');
end


% ----------------- OUTPUT ------------------------
%High-accuracy reference: x(1) = 1.153639002732
%       N            h          Euler x(1)       error
%      10           0.1      1.1274527973       0.02619
%      20          0.05      1.1400463413       0.01359
%      50          0.02      1.1480752578      0.005564
%     100          0.01      1.1508354119      0.002804
%     200         0.005      1.1522317190      0.001407
%     500         0.002      1.1530747644     0.0005642
%    1000         0.001      1.1533566621     0.0002823
%    2000        0.0005      1.1534977770     0.0001412
%    2500        0.0004      1.1535260133      0.000113
%    5000        0.0002      1.1535824991      5.65e-05
%   10000        0.0001      1.1536107487     2.825e-05
