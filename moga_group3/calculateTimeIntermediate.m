function totalTime = calculateTimeIntermediate(x,W_x,W_y,x0,y0,AirSpeed,sizeX,sizeY,METHOD)
% This calculates the time starting from a specified start point (x0,y0)
% It is used in the time-dependent version of this optimization problem.
%
% Copyright (c) 2012, MathWorks, Inc. 
%

L = numel(x);
X_points = linspace(x0,sizeX,L+2)';

Points = [X_points [y0; x; sizeY/2]];

PointList = [interp1(X_points,Points(:,1),linspace(x0,sizeX,101)',METHOD,'extrap')...
    interp1(X_points,Points(:,2),linspace(x0,sizeX,101)',METHOD,'extrap')];

PointList(PointList < 0) = 0;
PointList(:,1) = min(PointList(:,1),sizeX);
PointList(:,2) = min(PointList(:,2),sizeY);

dP = diff(PointList);

V_WindPath = [interp2(W_x,PointList(1:end-1,1)+1,PointList(1:end-1,2)+1,'*linear')...
    interp2(W_y,PointList(1:end-1,1)+1,PointList(1:end-1,2)+1,'*linear')];


% Dot product the wind (Vwind_path) with the direction vector (dP)
V_add = (sum(V_WindPath.*dP,2))./sqrt(sum(dP.^2,2));

dx = sqrt(sum(dP.^2,2)); %dx is the length of each subinterval in plist
dT = dx./(AirSpeed+V_add);  %dT = dP/dV
totalTime = 100*sum(dT);
