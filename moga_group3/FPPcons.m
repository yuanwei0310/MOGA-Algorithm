function [c ceq]=FPPcons(x)
sizeX = 50;
sizeY = 25;
xWayPoints=[0;x(1);x(2);x(3);x(4);x(5);50];
yWayPoints=[12.5;x(6);x(7);x(8);x(9);x(10);12.5];
PathPoints = WayPoints_To_Path([xWayPoints,yWayPoints],'linear',sizeX,sizeY,501);
R1=5*rand(501,1)+3*ones(501,1);

    for j=1:501
    c(j)=-((PathPoints(j,1)-35)^2+(PathPoints(j,2)-8)^2)+R1(j,1)^2;%+(x(2)-35)^2+(x(3)-35)^2+(x(4)-35)^2+(x(5)-35)^2+...
    %(x(6)-8)^2+(x(7)-8)^2+(x(8)-8)^2+(x(9)-8)^2+(x(10)-8)^2)+R1(1,i)^2;
    end
%end
ceq=[];
end