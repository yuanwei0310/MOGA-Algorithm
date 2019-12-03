function y=ZDT1_objfun(x)
numVar = 30;
g = 1 + 9*x(2)/(numVar-1);


y(1) = x(1);
y(2) = g*(1-sqrt(x(1)/g));