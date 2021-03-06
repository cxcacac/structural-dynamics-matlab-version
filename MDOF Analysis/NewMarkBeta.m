function [x,v,a] = NewMarkBeta(m, k, c, dt, F, X0, ks, R)
% use beta and gamma to replace linear acceleration formula;
% initial state: 3n*1;
% for linear model, ks = k, R = zeros(n,1)

cn = length(m);
dF = diff(F,1,2);
x(:,1) = X0(1:cn);
v(:,1) = X0(cn+1:2*cn);
a(:,1) = X0(2*cn+1:3*cn);

beta = 1/4;

for i = 1:length(F)-1
	k_ = k + (1/dt/dt/beta).*m +(1/2/dt/beta).*c;
	dF_ = dF(:, i) + (m.*1/beta/dt + 1/2/beta.*c)*v(:,i) + (1/2/beta.*m - dt*(1-1/4/beta).*c)*a(:,i);
	dx = inv(k_)*dF_;
	dv = (1/2/beta/dt).*dx - (1/2/beta).*v(:, i) + ((1-1/4/beta)*dt).*a(:,i);
	
    x(:, i+1) = x(:, i)+dx;
    v(:, i+1) = v(:, i)+dv;
    Fs(:, i+1) = ks*x(:, i+1) + R;
    a(:, i+1) = inv(m)*(F(:,i+1)-c*v(:,i+1)-Fs(:,i+1));
end