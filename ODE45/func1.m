% \dot{y1} = y2*y3;
% \dot{y2} = y1*y3;
% \dot{y3} = -2*y1*y2;
clear;
clc;
clf;
[t,y] = ode45(@func, [0,20], [0,0.5,0.5]);
plot(t,y(:,1),t,y(:,2),t,y(:,3));
legend('y1', 'y2', 'y3');

function dy = func(t,y)
dy = zeros(3,1);
dy(1) = y(2)*y(3);
dy(2) = -y(1)*y(3);
dy(3) = -2*y(1)*y(2);
end