clear;clc;
% tic, toc is used to compute the computing time
tic
options = odeset('reltol',1e-13);
tspan = [0,60];
[tspan,x]=ode45(@vibration,tspan,[0.025 0],options);
toc

figure
plot(tspan,x(:,1),tspan,x(:,2));
legend('Displacement','Acceleration')
%Time Domain
% fig = figure('Name','Position');
% plot(x(:,1))

%Purpose Function

function dxdt = vibration(~ ,x)
%for t= 0:0.01:40
%f = cos(10*t);
f=0;
c = 0.;
k = 2.5;
m = 1;
dxdt = [0;0];
dxdt(1) = x(2);
dxdt(2) =(1/m)*( f - c*x(2) - k*x(1));
end
