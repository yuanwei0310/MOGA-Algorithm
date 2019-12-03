function [c ceq]=RTNKcons(x)
P1=4*rand(500,1)-ones(500,1);
P2=4*rand(500,1)-ones(500,1);
for i=1:500
if( x(2) == 0)
    c(i)=-x(1)^2-x(2)^2+1+0.1*cos(16*atan(Inf))+0.2*sin(P1(i))*cos(P2(i));
else
    c(i)=-x(1)^2-x(2)^2+1+0.1*cos(16*atan(x(1)/x(2)))+0.2*sin(P1(i))*cos(P2(i));
end
end
ceq=[];
end