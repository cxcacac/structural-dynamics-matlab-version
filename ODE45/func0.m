% y' = 2*t
% [t,y] = ode45(odefun, tspan, y0)
% tspan = [t0 tf] integrate the system of differential equations 
% from t0 to tf with initial condition y0, 
% each row in y corresponds to value return in column vector t;

% example:
% y1' = y1 + 2*y2;
% y2' = 3*y1 + 2*y2;
% function dydt = odefun(t,y)
% dydt = zeros(2,1)
% dydt(1) = y(1) + 2*y(2);
% dydt(2) = 3*y(1) + 2*y(2);
% end

clear;
clc;
clf;
tspan = [0 5];
y0 = 0;
[t,y] = ode45(@(t,y) 2*t, tspan, y0);
plot(t,y,'-o');