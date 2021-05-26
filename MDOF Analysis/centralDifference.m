% step by step integration
% from x_i, x_{i-1} to get x_{i+1};
function [d,v,a] = centralDifference(m,k,c,dt,F,X0)

cn = length(m);
d(:,2) = X0(1:cn);
v(:,1) = X0(cn+1:2*cn);
a(:,1) = X0(2*cn+1:3*cn);

% x_{-1}, initial state;
d(:,1) = (eye(cn)-0.5.*dt.*dt.*inv(m)*k)*X0(1:cn)-(dt.*eye(cn) + 0.5.*dt.*dt.*inv(m)*c)*v(:,1)+...
        0.5.*dt.*dt.*inv(m)*F(:,1);

for i = 2:length(F)
    d(:, i+1) = inv(1./dt./dt.*m+0.5./dt.*c)*(F(:,i)-(k-2./dt./dt.*m)*d(:,i)-...
                (1./dt./dt.*m-0.5./dt.*c)*d(:,i-1));
    v(:,i) = 0.5./dt.*(d(:,i+1)-d(:,i-1));
    a(:,i) = 1./dt./dt.*(d(:,i+1)-2.*d(:,i)+d(:,i-1));
end
% delete the first row;
d(:, 1) = [];