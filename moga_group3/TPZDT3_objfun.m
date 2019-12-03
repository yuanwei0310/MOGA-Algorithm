function y= TPZDT3_objfun(x)
% Objective function : Test problem 'ZDT3'.
%*************************************************************************


g = 1 + 9*x(2)/29;


y(1) = x(1);
y(2) = g*(1-sqrt(x(1)/g)-x(1)/g*sin(10*pi*x(1)));
