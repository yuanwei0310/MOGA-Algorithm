function [c ceq]=TPCTP_objfun(x)
th=-0.2*pi; a=0.2; b=10; d=6; e=1;
g=abs(1+(sum(x(2:10))).^(0.25));
y(1)=x(1);
y(2)=g*(1-sqrt(x(1)/g));
c=-cos(th)*(y(2)-e)+sin(th)*y(1)+a*(abs(sin(b*pi*(sin(th)*(y(2)-e)+cos(th)*y(1)).^1))).^6;
ceq=[];