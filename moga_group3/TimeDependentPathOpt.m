%  This demo shows continuous optimization through a time dependent
%  vector field. This is an extension of the time-invariant optimization
%  demo.
%
% Copyright (c) 2012, MathWorks, Inc. 
%

%% Initialize Parameters

clear; clc; clf;

rng(50);
AirSpeed = 500;
sizeX = 50;
sizeY = 25;
numWayPoints = 5;
METHOD = 'spline';

W_x0 = makeWindFun(sizeX,sizeY);
W_y0 = makeWindFun(sizeX,sizeY);
W_x1 = makeWindFun(sizeX,sizeY);
W_y1 = makeWindFun(sizeX,sizeY);

[Xgrid,Ygrid] = meshgrid(0:sizeX,0:sizeY);
hq = quiver(Xgrid,Ygrid,W_x0,W_y0,'k');
ax = gca;
xlabel('Units = 100 [km]');
axis equal tight

%% Add some color to make it more visible
L = (sqrt((Xgrid-sizeX).^2 + (Ygrid-sizeY/2).^2));
Favorability = ((sizeX-Xgrid).*W_x0 +  (sizeY/2-Ygrid).*W_y0)./L;
Favorability(~isfinite(Favorability)) = 0;
Favorability0 = Favorability; %Record this so we can reuse it later

Favorability = ((sizeX-Xgrid).*W_x1 +  (sizeY/2-Ygrid).*W_y1)./L;
Favorability(~isfinite(Favorability)) = 0;
Favorability1 = Favorability; %Record this so we can reuse it later

hold on;
h_im = imagesc(Favorability0); % This will be the background for the vector field
set(h_im,'Xdata',[0 sizeX],'Ydata',[0 sizeY]);
uistack(h_im,'bottom');

title('Time-varying Path Optimization')

% Change the colormap...
colormap(interp1([0,1,2],[1 0 0; 1 1 1; 0 1 0],0:0.01:2));
caxis(max(abs([Favorability0(:); Favorability1(:)]))*[-1 1]);

h_colorbar = colorbar;
title(h_colorbar,'Tailwind (km/h)')
xlabel(h_colorbar,'Headwind')
axis([-1 sizeX+1 -1 sizeY+1]);
axis equal tight
drawnow;

%% Re-calculate the solution at each "time step"
% In this optimization, we will assume that the x location of each waypoint
% is equally spaced, and we only optimize over the y values.

opts = optimset('fmincon');
opts.Display = 'off';
opts.Algorithm = 'active-set';

X_points = linspace(0,sizeX,numWayPoints+2)';
x_opt = sizeY/2*ones(numWayPoints,1);
ystart = sizeY/2;

kctr = 1;
timeStep = 1/40;
for time = 0: timeStep : 1-timeStep
    fprintf('%02.f%% Percent Complete\n',kctr*timeStep*100);
    W_x = (1-time)*W_x0 + time*W_x1;
    W_y = (1-time)*W_y0 + time*W_y1;
    
    % Make initial conditions based on the previous solution
    xstart = sizeX*time;
    Points = [X_points [ystart; x_opt; sizeY/2]];
    X_points_old = X_points;
    X_points = linspace(xstart,sizeX,numWayPoints+2)';
    xy0 = interp1(X_points_old,Points(:,2),X_points(2:end-1),METHOD,'extrap');
    ystart = interp1(X_points_old,Points(:,2),xstart,METHOD,'extrap');
    
    lb = zeros(size(xy0(:)));
    ub = sizeY*ones(numWayPoints,1);
    
    % Calculate the new optimization problem
    objectiveFun = @(P) ...
        calculateTimeIntermediate(P,W_x,W_y,xstart,ystart,AirSpeed,sizeX,sizeY,'spline');
    x_opt = fmincon(objectiveFun, xy0(:), [],[],[],[],lb,ub,[],opts);
    hold on;
    
    Points = [X_points [ystart; x_opt; sizeY/2]];
    PointList = WayPoints_To_Path(Points,METHOD,sizeX,sizeY,101);

    K_PointListX{kctr} = PointList(:,1);
    K_PointListY{kctr} = PointList(:,2);
    K_k{kctr} = time;
    kctr = kctr+1;
    if time == 0, plot(PointList(:,1), PointList(:,2),'b','linewidth',2);
        drawnow;
    end
end
%% Plot the results
pause;

for kk = 1:kctr-1
    time = K_k{kk};
    W_x = (1-time)*W_x0 + time*W_x1;
    W_y = (1-time)*W_y0 + time*W_y1;
    
    set(hq,'Udata',W_x,'Vdata',W_y);
    try, delete(hthis), end
    hthis = plot(K_PointListX{kk}, K_PointListY{kk},'k','linewidth',2);
    if kk == 1, plot(K_PointListX{kk}, K_PointListY{kk},'b','linewidth',2);
        pause(1);
    end
    
    %Update the background image, as well as the path line
    plot(K_PointListX{kk}(1), K_PointListY{kk}(1),'k.','markersize',16);
    Favorability = ((sizeX-Xgrid).*W_x +  (sizeY/2-Ygrid).*W_y)./L;
    Favorability(~isfinite(Favorability)) = 0;
    set(h_im,'Cdata',Favorability);
    drawnow
end

