clear all; close all;
plotDir = fullfile(fileparts(mfilename('fullpath')), 'plots');
if ~isfolder(plotDir)
    mkdir(plotDir);
end
%% circle map
Om = 0.2; % omega
K = 1;
th = linspace(0,1,1000);
figure;
hold on
plot(th, mod(th + Om - (K/(2*pi))*sin(2*pi*th), 1),'.');
xlabel('\theta_n');
ylabel('\theta_{n+1}');
axis square
exportgraphics(gcf, fullfile(plotDir, 'matlab_circle_map.png'), 'Resolution', 300);

%% staircase

Omega = linspace(0,1,20000);
N = 14000; % iterations
theta = zeros(size(Omega));
for k = 1:N
    theta = theta + Omega - (K/(2*pi))*sin(2*pi*theta);
end
W = theta / N;

figure;
plot(Omega, W, '.')
xlabel('\Omega');
ylabel('W');
axis square
exportgraphics(gcf, fullfile(plotDir, 'matlab_staircase.png'), 'Resolution', 300);
%% dimension calculation


% finding flat segments
tol = 1e-5;
flat = [0, abs(diff(W)) < tol, 0];
starts = find(diff(flat) == 1);
ends   = find(diff(flat) == -1);
widths = Omega(ends) - Omega(starts);


r  = linspace(2e-3, 0.1, 101);
sr = arrayfun(@(rk) sum(widths(widths >= rk)), r); % iterating over r(max width) values
Nr = (1 - sr) ./ r; % amount of empty space

x = log(1./r);
y = log(Nr);
p = polyfit(x, y, 1);
fprintf('slope = %.3f', p(1));

figure;
plot(x, y, 'b.','DisplayName', '$N(r)$ from plateau widths');
hold on
plot(x, polyval(p, x), 'r-', 'LineWidth', 2, 'DisplayName', sprintf('$y$ = %.4f $x$ + %.4f', p(1), p(2)));
grid on;
xlabel('log 1/r');
ylabel('log N(r)');
legend('show', 'Location', 'best', 'Interpreter', 'latex');
exportgraphics(gcf, fullfile(plotDir, 'matlab_dimension_fit.png'), 'Resolution', 300);
