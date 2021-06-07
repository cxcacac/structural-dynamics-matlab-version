function [x,v,a] = NewMarkBeta(M,K,C,F,dt,X0)
% Variable Description :
% INPUT :
%       M - Mass Matrix (in modal coordinates)
%       K - Stiffness Matrix (in modal coordinates)
%       C - Damping Matrix (in modal coordinates)
%       P - Force Matrix (in modal coordinates)
% OUTPUT :
%        x - modal displacement's 
%        v - modal velocities
%        a - modal accelerations 
%        U - system's displacement

cn = length(M); % dof;
n = length(F); % number of dt;

% when beta = 1/6, it is linear acceleration;
% when beta >= 1/4, it is unconditional stable;
beta = 1/4;
gamma = 1/2; 
a0=1/beta/dt^2; a1=gamma/beta/dt; a2=1/beta/dt; a3=1/2/beta-1;
a4=gamma/beta-1; a5=dt*(gamma/beta-2)/2; a6=dt*(1-gamma); a7=gamma*dt;

% initial state
x = zeros(cn, n); 
v = zeros(cn, n); 
a = zeros(cn, n); 
x(:,1) = X0(1:cn); 
v(:,1) = X0(cn+1:2*cn);
a(:,1) = X0(2*cn+1:3*cn);

K_ = K + a0*M + a1*C;
iK = inv(K_);

for i = 1:length(F)-1
    F_(:,i+1)=F(:,i+1)+M*(a0*x(:,i)+a2*v(:,i)+a3*a(:,i))+C*(a1*x(:,i)+a4*v(:,i)+a5*a(:,i));
    x(:,i+1)=iK*F_(:,i+1);
    a(:,i+1)=a0*(x(:,i+1)-x(:,i))-a2*v(:,i)-a3*a(:,i);
    v(:,i+1)=v(:,i)+a6*a(:,i)+a7*a(:,i+1);
end