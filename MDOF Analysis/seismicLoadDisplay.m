clc;
clf;
m = 96.6;
wavefile = char('elcentro_NS.dat');
ugmax = 0.7;
[ug, t, tf,dt] = wave(wavefile, ugmax);

figure(1);
plot(t, ug, 'r-');
xlabel('time/s');
ylabel('acc/(ms-2)');
