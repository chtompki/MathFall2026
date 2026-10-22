% Problem 6.6.4. Run this script to compare phase portraits.
clear; close all;
fields = {@(t,z)[z(2);3-z(1)-z(2)^2], ...
          @(t,z)[z(2)-z(2)^3;z(1)*cos(z(2))], ...
          @(t,z)[sin(z(2));z(2)^2-z(1)]};
windows = {[-1 5 -3 3],[-3 3 -4 4],[-2 5 -3 3],[-5 48 -8 8]};
seeds = {[3,.2;3,.7;3,1.2;1,0;4,0], ...
 [0,.3;0,.9;0,1.3;.2,0;1,2;1,-2;0,3.5], ...
 [.2,0;.8,0;1.5,0;3,0;4,2;4,-2], ...
 [.5,0;3,0;10,3;10,-3;39,6;39,-6;20,4;20,-4]};
for j=1:4
    f=fields{min(j,3)}; w=windows{j};
    figure('Color','w'); hold on;
    [X,Y]=meshgrid(linspace(w(1),w(2),25),linspace(w(3),w(4),25));
    U=zeros(size(X)); V=U;
    for k=1:numel(X)
        d=f(0,[X(k);Y(k)]); U(k)=d(1); V(k)=d(2);
    end
    N=hypot(U,V); N(N==0)=1;
    quiver(X,Y,U./N,V./N,.5,'Color',[.75 .75 .75]);
    opts=odeset('RelTol',1e-8,'AbsTol',1e-10,...
        'Events',@(t,z)leaveWindow(t,z,w));
    for k=1:size(seeds{j},1)
        for direction=[-1,1]
            [~,z]=ode45(f,[0,direction*40],seeds{j}(k,:)',opts);
            plot(z(:,1),z(:,2),'b-');
        end
    end
    axis(w); xlabel('x'); ylabel('y'); grid on;
    title(sprintf('6.6.4: panel %d (fourth panel is a wider view of c)',j));
end
function [value,isterminal,direction]=leaveWindow(~,z,w)
value=min([z(1)-w(1),w(2)-z(1),z(2)-w(3),w(4)-z(2)]);
isterminal=1; direction=0;
end
