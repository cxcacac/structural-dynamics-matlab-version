function [force] = getDamperForce(x, v, a, choice)
c_d = 5.07;
if(strcmp(choice, 'passive'))
    force = -c_d*v;
elseif(strcmp(choice, 'freeVibration'))
    force = 0;
else 
    force = 1;
end
