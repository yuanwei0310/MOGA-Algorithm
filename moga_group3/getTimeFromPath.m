function TravelTime = getTimeFromPath(PathPoints,W_x,W_y,AirSpeed,sizeX,sizeY,METHOD)
% This is the main function that actually calculates the line integral
% along the path. This is the objective function for the optimizer.
%
% Copyright (c) 2012, MathWorks, Inc. 
%

% If we are called from the optimization routine (caller = 'optimizer')
% then we need to interpolate the fine path from the input control points.
if isvector(PathPoints)
    PathPoints = [0 sizeY/2; reshape(PathPoints,2,[])'; sizeX sizeY/2];
    PathPoints = WayPoints_To_Path(PathPoints,METHOD,sizeX,sizeY,101);
end

dP = diff(PathPoints);

% Interpolate the wind vector field at all the points in PathPoints.
V_wind = [interp2(W_x,PathPoints(1:end-1,1)+1,PathPoints(1:end-1,2)+1,'linear') interp2(W_y,PathPoints(1:end-1,1)+1,PathPoints(1:end-1,2)+1,'*linear')];

% Dot product the wind (Vwind_path) with the direction vector (dP) to get
% the tailwind/headwind contribution
V_add = (sum(V_wind.*dP,2))./sqrt(sum(dP.^2,2));
dx = sqrt(sum(dP.^2,2))*100; %dx is the length of each subinterval in PathPoints
dt = dx./(AirSpeed+V_add);  %dT = dP/dV
TravelTime = sum(dt);