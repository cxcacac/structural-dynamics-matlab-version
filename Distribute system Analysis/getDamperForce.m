function [force] = getDamperForce(x, v, a, choice)
% non-dimensional in johnson's paper(2007), c = m*w0*L*c_d
% c_d = 220;
c_d = 5.07*85*0.6575*220; 
if(strcmp(choice, 'passive'))
    force = -c_d*v;
elseif(strcmp(choice, 'freeVibration'))
    force = 0;
else 
    force = 1;
end
