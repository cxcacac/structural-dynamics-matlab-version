% second order function: van der pol
% y1''-u(1-y1^2)y1'+y1 =0
% let y1'= y2;
% y2' = u(1-y1^2)y2-y1
clear;
clc;
clf;
[t,y] = ode45(@vdp1,[0 20],[2; 0]);
plot(t,y(:,1), '-o', t, y(:,2),'-o')
title("Solution of van der pol equation with ODE45");
xlabel('Time t');
ylabel('Solution y');
legend('y_1', 'y_2');

function dydt = vdp1(t,y)
dydt = [y(2); (1-y(1)^2)*y(2)-y(1)];
end