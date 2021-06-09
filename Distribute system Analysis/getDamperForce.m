function [force] = getDamperForce(x, v, a, choice)
% non-dimensional in johnson's paper(2007), c = m*w0*L*c_d
% c_d = 220;
c_d = 5.07*85*0.6575*220; 
% in sapinski' paper(2002)
Fc = 580.31; 
c0 = 56.77;
f0 = -176.62;
if(strcmp(choice, 'passive'))
    force = -c_d*v;
elseif(strcmp(choice, 'free-vibration'))
    force = 0;
elseif(strcmp(choice, 'semi-active'))
    force = -(Fc*sign(v) + c0*v + f0);
end
