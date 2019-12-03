A=load('FPP.txt');
x=A(:,1);
y=A(:,2);
a=plot(x,y,'o');
hold on
B=load('FPP9.txt');
x=B(:,1);
y=B(:,2);
b=plot(x,y,'*');
hold on
C=load('FPP64.txt');
x=C(:,1);
y=C(:,2);
b=plot(x,y,'.');
