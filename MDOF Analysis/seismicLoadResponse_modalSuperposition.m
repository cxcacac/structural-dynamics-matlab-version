clear;
m = [2762,2760,2300];
k = [2.485,1.921,1.522].*1e4;
es = 0.05;

wavefile = char('elcentro_NS.dat');
ugmax = 0.7;
[ug, t, tf,dt] = wave(wavefile, ugmax);

[M] = lumpMass(m);
[K] = stiffnessShear(k);
[E,F] = waveForce(ug, M);
flag = 1;
[C,T,z] = dampR(K,M,E,es,flag);

cn = length(m);
X0 = zeros(3*cn, 1);
[d1,v1,a1] = LinearAcceleration(M,K,C,dt,F,X0);
[d,v,a] = modalSuperposition(M,K,C,F,dt);

figure(1);
plot(t,d(1,:),'r-', t, d1(1,:),'k:');

xlabel('time/s');
ylabel('displacement/m');