function y=TPOSY_objfun(x)
y(1)=-(25*(x(1)-2)^2+(x(2)-2)^2+(x(3)-1)^2+(x(4)-4)^2+(x(5)-1)^2);
%y(2)=0;
%for i=1:6
    y(2)=x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(5)^2+x(6)^2;
end