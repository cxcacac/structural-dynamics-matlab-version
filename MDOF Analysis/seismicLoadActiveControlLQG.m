clc;
clf;

m = [2762,2760,2300];
k = [2.485,1.921,1.522].*1e4;
% m = [4,4,4].*1e5;
% k = [2,2,2].*1e8;
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
% when a2 = 1e-8, ruge-kutta method fails;
a1 = 100; 
a2 = 0.01e-6;
Q = a1.*[K,zeros(cn); zeros(cn), M];
R = a2.*eye(cn);

G = lqr(A,B,Q,R);
Au = A - B*G;

X0 = zeros(2*cn, 1);

% [d,v] = stateSpaceMatrixTransfer(Au, H, F, X0, dt);
[d0, v0] = stateSpaceMatrixTransfer(A, H, F, X0, dt);
% u = -G*[d;v];

w_cn = 3; % number of parameters observed;
yv = ones(cn,1)*ug; % absolute acceleration on each floors;
Cv = [zeros(cn), eye(cn)]*A;
Dv = [zeros(cn), eye(cn)]*B;
Hv = 2.*inv(M);
sys1 = ss(A, B, Cv, Dv);

Qe = 10e-4;
Re = 10e-2.*eye(w_cn);
[kest, Ge, p] = kalman(sys1, Qe, Re);

Ag = A - B*G - Ge*Cv + Ge*Dv*G;
Bg = [Ge, H-Ge*Hv];
Fg = [yv; F];
[dg, vg] = stateSpaceMatrixTransfer(Ag, Bg, Fg, X0, dt);
u_lqr = -G*[dg; vg];

figure(1)
plot(t, dg(1,:),'r-', t, d0(1,:), 'k--');
xlabel('time/s');
ylabel('displacement/m');

figure(2)
plot(t, vg(1,:),'r-', t, v0(1,:), 'k--');
xlabel('time/s');
ylabel('velocity/(ms-1)');

figure(3)
plot(t, u(1,:));
xlabel('time/s');
ylabel('Force/N');