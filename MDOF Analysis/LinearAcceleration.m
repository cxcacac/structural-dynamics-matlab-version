function [x,v,a] = linearAcceleration(m,k,c,dt,F,X0)
cn = length(m);
% diff function
dF = diff(F,1,2);

x(:,1) = X0(1:cn);
v(:,1) = X0(cn+1:2*cn);
a(:,1) = X0(2*cn+1:3*cn);

for i=1:length(F)-1
	k_ = k+(6/dt/dt).*m+(3/dt).*c;
	dF_ = dF(:,i)+(6/dt.*m+3.*c)*v(:,i)+(3.*m+dt/2.*c)*a(:,i);
	dx = inv(k_)*dF_;
	dv = (3/dt).*dx-3.*v(:,i)-dt/2.*a(:,i);
	x(:,i+1) = x(:,i)+dx;
	v(:,i+1) = v(:,i)+dv;
    % use dynamic equation to get acceleration;
	a(:,i+1) = inv(m)*(F(:,i+1)-c*v(:,i+1)-k*x(:,i+1));
end