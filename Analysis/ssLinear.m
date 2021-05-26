function [A,B,D,L] = ssLinear(M,K,C);
% cn, degree of freedom;
cn =length(m);
A = [zeros(cn), eye(cn); -inv(M)*K, -inv(M)*C];
B = [zeros(cn), inv(M)];
D = [eye(cn), zeros(cn); zeros(cn), eye(cn); -inv(M)*K, -inv(M)*C];
L = [zeros(cn); zeros(cn); inv(M)];