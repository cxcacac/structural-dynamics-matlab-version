clear;
m = [2762,2760,2300];
k = [2.485,1.921,1.522].*1e4;

[M] = lumpMass(m);
[K] = stiffnessShear(k);
E = ones(length(m),1);
es = 0.05;
flag = 0;
[C,T,z] = dampR(K,M,E,es,flag);


