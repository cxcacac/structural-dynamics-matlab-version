% dynamic equation is ordinary differential equations;
% initial state: X(0)=0, X'(0) = 0;
% use runge-kutta method to solve ode;

function [d,v] = rk4(odefile, tspan, dt, X0, A, B, F)
% fourth order rk;
% tspan, [t0, tf]; dt, interval; X0: initial state;
% A, B, coefficient matrix;
cn = length(X0)./2;
F0 = B*F;

t0 = tspan(1);
tf = tspan(2);
X = X0;
i = 1;
for t = t0:dt:tf-dt
    k1 = dt.*feval(odefile, X0, F0(:i), A);
    k2 = dt.*feval(odefile, X0+k1./2, F0(:i), A);
    k3 = dt.*feval(odefile, X0+k2./2, F0(:i), A);
    k4 = dt.*feval(odefile, X0+k3, F0(:i), A);
    X0 = X0 + (k1+2.*k2+ 2.*k3 + k4)./6;
    X = [X, X0];
    i = i+1;
end
d = X(1:cn, :);
v = X(cn+1:2*cn, :);