function y = shapeFunction(x, num, dL, cL)
% x is coordinate along the cable
% num is rank of modes;
if(num==1)
    y = piecewise(x<=dL,x./dL,x>dL, (cL-x)./(cL-dL));
else
    y = sin(num*x);
end