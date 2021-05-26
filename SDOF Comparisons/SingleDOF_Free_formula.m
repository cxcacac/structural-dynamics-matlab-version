clc;
close all;

x0= 0.025; % m
v0= 0.000; % m/s
t = linspace(0,60,1000);
m = 1;     % kg
c = 0.1;   % N*sec/meter
k = 2.5;   % N/m

omega_n = sqrt(k/m);
zeta    = 0.5*c/m/omega_n;

x = exp(-zeta*omega_n*t).*(x0*cos(sqrt(1-zeta^2)*omega_n*t) +...
       (v0+zeta*omega_n*x0)/sqrt(1-zeta^2)/omega_n*...
       (sin(sqrt(1-zeta^2)*omega_n*t)));

figure(1)
plot(t,x.*1000,'r-','linewidth',1.5)
grid on;
title('\fontsize{10}\fontname{Times New Roman}Damped Response of a Single Degree of Freedom System');
xlabel('\fontsize{10}\fontname{Times New Roman}\it Time \rm / s');
ylabel('\fontsize{10}\fontname{Times New Roman}\it Displacement \rm / mm');
legend('Time domain response of SDF system')
set(gcf,'unit','centimeters','position',[0 10 13.53 9.03],'color','white');