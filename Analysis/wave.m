function [ug, t, tf, dt] = wave(wavefile, ugmax)
% input: 
% wavefile, filename; ugmax: maximum acceleration
% output: 
% ug, modified acceleration; t, time series; tf, time range; dt, interval

Dat = dlmread(wavefile,'');
ug1 = Dat(:,2);
t1 = Dat(:,1);
tf = Dat(length(t1),1);
dt = t1(3)-t1(2);
ug = ugmax./max(abs(ug1)).*ug1';
t = t1';
end