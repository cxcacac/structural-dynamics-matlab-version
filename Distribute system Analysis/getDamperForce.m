function [force] = getDamperForce(x, v, a, choice)
% non-dimensional in johnson's paper(2007), c = m*w0*L*c_d
% c_d = 220;
c_d = 5.07*85*0.6575*220; 
%%  in Spencer's paper 1997 
% use MRX-135GD, Bingham model, 2.5Hz sinusoidal load, constant 1.5V;
Fc = 670; 
c0 = 5000;
f0 = -95;
%% parameters in different mode;
if(strcmp(choice, 'passive'))
    force = -c_d*v;
elseif(strcmp(choice, 'free-vibration'))
    force = 0;
elseif(strcmp(choice, 'semi-active'))
    force = -(Fc*sign(v) + c0*v + f0);
end
