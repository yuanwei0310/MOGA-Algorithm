function y=FPP_objfun(x)
AirSpeed = 500;
sizeX = 50;
sizeY = 25;
numWayPoints = 5;
rng(50);
W_x = makeWindFun(sizeX,sizeY);
W_y = makeWindFun(sizeX,sizeY);
xWayPoints=[0;x(1);x(2);x(3);x(4);x(5);50];
yWayPoints=[12.5;x(6);x(7);x(8);x(9);x(10);12.5];
PathPoints = WayPoints_To_Path([xWayPoints,yWayPoints],'linear',sizeX,sizeY,501);
y(1)=getTimeFromPath(PathPoints,W_x,W_y,AirSpeed,sizeX,sizeY,501);
%for i=1:101
    %y()
y(2)=-((x(1)-35)^2+(x(2)-35)^2+(x(3)-35)^2+(x(4)-35)^2+(x(5)-35)^2+...
    (x(6)-8)^2+(x(7)-8)^2+(x(8)-8)^2+(x(9)-8)^2+(x(10)-8)^2);