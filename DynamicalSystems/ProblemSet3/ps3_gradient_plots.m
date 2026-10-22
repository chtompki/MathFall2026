% Problem 7.2.9: separate phase and equipotential graphs.
clear; close all;
[X,Y]=meshgrid(linspace(-2,2,25));
figure('Color','w');
U=2*X; W=8*Y; N=hypot(U,W); N(N==0)=1;
quiver(X,Y,U./N,W./N); hold on;
t=linspace(-2,.2,500);
for x0=[-.4,.4]
    for y0=[-.1,.1]
        plot(x0*exp(2*t),y0*exp(8*t));
    end
end
axis([-2 2 -2 2]); axis square; title('7.2.9(b): phase portrait'); xlabel('x'); ylabel('y');
figure('Color','w'); contour(X,Y,-X.^2-4*Y.^2,12); axis equal;
title('7.2.9(b): equipotentials'); xlabel('x'); ylabel('y');
figure('Color','w'); N=hypot(X,Y); N(N==0)=1;
quiver(X,Y,-X./N,-Y./N); axis equal;
title('7.2.9(c): radial inward flow (normalized arrows)'); xlabel('x'); ylabel('y');
figure('Color','w'); contour(X,Y,exp(X.^2+Y.^2),exp([.2,.5,1,2,3,4])); axis equal;
title('7.2.9(c): equipotentials'); xlabel('x'); ylabel('y');
