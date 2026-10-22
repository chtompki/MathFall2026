% Problem 7.5.7: stiff integration of oscillations and excitable responses.
clear; close all;
b=100; B=1e-5;
u1=sqrt(2*B*(1+b*B)/(1-2*b*B+sqrt(1-8*b*B))); % stable fold formula
u2=sqrt((1-2*b*B+sqrt(1-8*b*B))/(2*b));
h=@(u)u+u./(b*(B+u.^2));
f=@(c) @(t,z)[b*(z(2)-z(1))*(B+z(1)^2)-z(1);c-z(1)];
opts=odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',.1);
fprintf('b=%g, B=%g, 8bB=%g; folds c1=%g, c2=%g\n',b,B,8*b*B,u1,u2);
c=(u1+u2)/2;
[t,z]=ode15s(f(c),[0 500],[c;h(c)*1.01],opts);
figure('Color','w');
subplot(1,2,1); u=linspace(0,2,6000); plot(u,h(u),'k'); hold on;
xline(c,'--'); plot(z(:,1),z(:,2),'b');
xlim([0 2]); ylim([0 2]); xlabel('u'); ylabel('v'); title('Relaxation orbit and nullclines');
subplot(1,2,2); plot(t,z(:,1)); xlim([400 500]); xlabel('t'); ylabel('u'); title('Late-time oscillations');
c=.95*u1; eq=[c;h(c)];
figure('Color','w'); hold on;
for delta=[1e-4,.05]
    [t,z]=ode15s(f(c),[0 200],eq+[delta;0],opts);
    plot(t,z(:,1),'DisplayName',sprintf('Initial u perturbation = %g',delta));
end
xlabel('t'); ylabel('u'); title('Below lower fold: small response versus excitable pulse'); legend show;
