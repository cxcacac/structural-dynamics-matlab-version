% pass extra parameters to ode function
% y'' = A/B * t * y;
% y1' = y2;
% y2' = A/B * t * y;
clear;
clc;
clf;

A = 1;
B = 2;
tspan = [0 5];
y0 = [0 0.01];
[t,y] = ode45(@(t,y) odefcn(t,y,A,B), tspan, y0);
plot(t,y(:,1),'-o',t,y(:,2),'-.')

function dydt = odefcn(t,y,A,B)
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = (A/B)*t.*y(1);
end
