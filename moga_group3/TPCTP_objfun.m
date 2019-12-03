function y=TPCTP_objfun(x)
th=-0.2*pi; a=0.2; b=10; d=6; e=1;
y(1)=x(1);
g=abs(1+(sum(x(2:10))).^(0.25));
y(2)=g*(1-sqrt(x(1)/g));