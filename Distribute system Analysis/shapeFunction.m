function y = shapeFunction(x, num, dL, cL)
% x is coordinate along the cable, dL: damper location, cL: cable length;
% num is rank of modes;
if(num==0)
    y = piecewise(x<=dL, x./dL, x>dL,(cL-x)./(cL-dL));
else
    y = sin(num*pi*x./cL);
end