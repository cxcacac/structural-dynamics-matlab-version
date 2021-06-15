
clear;
clc;
x = 0;
dotx = 0;

[t,y] = ode45(@(t,y) boucWen(x, dotx, 0, y),[0 0.02],[2; 0]);
% inner displacement
innerd = y(length(y),1);

function dydt = boucWen(x,dotx,voltage,y)
gamma = 1;
n = 2;
beta = 1;
A = 1;
c0 = 1;
c1 = 1;
alpha = 1;
k0 = 1;
% y(1) -> y, y(2) -> z;
dydt = zeros(2,1);
dydt(1) = (1/(c0+c1))*(alpha*y(2) + c0*dotx + k0*(x-y(1)));
dydt(2) = -gamma*abs(dotx-(1/(c0+c1))*(alpha*y(2) + c0*dotx + k0*(x-y(1))))...
            *abs(y(2))^(n-1)*y(2) - beta*(dotx-(1/(c0+c1))*(alpha*y(2) + c0*dotx + k0*(x-y(1))))...
            *abs(y(2))^n + A*(dotx - (1/(c0+c1))*(alpha*y(2) + c0*dotx + k0*(x-y(1))));
end