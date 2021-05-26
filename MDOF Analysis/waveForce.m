function [E, F] = waveForce(ug, M)
% ln: dimension of seismic load;
ln = size(ug,1);
cn = length(M);

E = ones(cn, ln);

F = -M*E*ug;