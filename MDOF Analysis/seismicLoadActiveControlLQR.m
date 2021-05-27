clc;
clf;

m = [4,4,4].*1e5;
k = [2,2,2].*1e8;
cn = length(m);
location = ones(1,cn); % location of actuator
es = 0.05;

wavefile = char('elcentro_NS.dat');
ugmax = 0.7;
[ug,t,tf,dt] = wave(wavefile, ugmax);

% \dot{U} = AU + Bu + HF

[M] = lumpMass(m);
[E,F] = waveForce(ug, M);
[K] = stiffnessShear(k);
[C,T,z] = dampR(K,M,E,es,1);
[D] = relativeK(location);

[A,H,Dd,Ld] = ssLinear(M,K,C);
B = [zeros(cn); inv(M)*D];

% weight matrix;
a1 = 100; 
a2 = 0.01e-6;
Q = a1.*[K,zeros(cn); zeros(cn), M];
R = a2.*eye(cn);

G = lqr(A,B,Q,R);
Au = A - B*G;

X0 = zeros(2*cn, 1);
% [d,v] = rk4('stateSpaceEq', [0, tf], dt, X0, Au, H, F);
% [d0, v0] = rk4('stateSpaceEq', [0, tf], dt, X0, A, H, F);

[d,v] = stateSpaceMatrixTransfer(Au, H, F, X0, dt);
[d0, v0] = stateSpaceMatrixTransfer(A, H, F, X0, dt);
u = -G*[d;v];

figure(1)
plot(t, d(1,:),'r-', t, d0(1,:), 'k--');
xlabel('time/s');
ylabel('displacement/m');

figure(2)
plot(t, v(1,:),'r-', t, v0(1,:), 'k--');
xlabel('time/s');
ylabel('velocity/(ms-1)');

figure(3)
plot(t, u(1,:));
xlabel('time/s');
ylabel('Force/N');