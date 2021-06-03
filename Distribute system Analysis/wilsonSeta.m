function [d,v,a] = wilsonSeta(m,k,c,dt,F,X0,ks,R)
% wilson-seta is modification of linear acceleration;
% assume [t, t+seta*delta t] is linear change;
% m,k,c: mass, stiffness, damping matrix;

cn = length(m);
nb = length(F);
beta = 1/6; % coefficient;
seta = 1.4; % seta means extension range, tau = 1.4*dt;
tao = seta*dt;
dF1 = diff(F,1,2);
dF = dF1(:,1:nb-2)+(seta-1).*dF1(:, 2:nb-1);

% 初始状态
d(:, 1) = X0(1:cn);
v(:, 1) = X0(cn+1: 2*cn);
a(:, 1) = X0(2*cn+1:3*cn);

for i=1:nb-2
	k_ = k+(1/tao/tao/beta).*m+(1/2/tao/beta).*c;
	dF_ = dF(:,i) + (m.*1/beta/tao + 1/2/beta.*c)*v(:,i) + (1/2/beta.*m-tao*(1-1/4/beta).*c)*a(:,i);
	dx_ = inv(k_)*dF_;
	da_ = (6/tao/tao).*dx_-(6/tao).*v(:,i)-3.*a(:,i);
	
	da = da_./seta;
	dv = a(:,i).*dt + da.*dt./2;
	dx = v(:,i).*dt + a(:,i).*dt.*dt./2 +da.*dt.*dt/6;
	
	d(:,i+1) = d(:, i) + dx;
	v(:, i+1) = v(:, i) + dv;
	Fs(:, i+1) = ks*d(:, i+1) + R;
	a(:, i+1) = inv(m)*(F(:, i+1)-c*v(:, i+1)-Fs(:,i+1));
end