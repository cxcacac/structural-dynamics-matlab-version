function y = shapeFunction(x, num, damperLocation, cableLength)
% x is coordinate along the cable
% num is rank of modes;
if(num==0)
    y = (x./damperLocation).*(x<damperLocation) + ((cableLength-x)./(cableLength-damperLocation)).*(x>=damperLocation);
else
    y = sin(num*x);
end