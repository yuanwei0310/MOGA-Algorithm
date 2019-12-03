function [c ceq]=TPTNKcons(x)
if( x(2) == 0)
    c = -x(1)^2 - x(2)^2 + 1 + 0.1 * cos( 16 * atan(Inf) );
else
    c = -x(1)^2 - x(2)^2 + 1 + 0.1 * cos( 16 * atan(x(1)/x(2)) );
end
ceq=[];