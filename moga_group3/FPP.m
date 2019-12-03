%% First make some random vector field of wind, and set parameters

clear; clf;
close all

AirSpeed = 500;
sizeX = 50;
sizeY = 25;
numWayPoints = 5;
rng(50);

W_x = makeWindFun(sizeX,sizeY);
W_y = makeWindFun(sizeX,sizeY);


nvars = 10;
lb = [0 0 0 0 0 0 0 0 0 0];
ub = [50 50 50 50 50 25 25 25 25 25];
options = gaoptimset;
options = gaoptimset(options,'PopulationSize', 50);
%options = gaoptimset(options,'Generations', 500);
options = gaoptimset(options,'Display', 'iter');
options = gaoptimset(options,'PlotFcns',@gaplotpareto);


options = gaoptimset(options, 'ParetoFraction', 0.8);
objfunction=@(x)FPP_objfun(x);
[x,fval]=gamultiobj(objfunction,nvars,[],[],[],[],lb,ub,@FPPcons,options)
fval1=min(fval(:,1));
fval2=min(fval(:,2));
row1=find(fval(:,1)==fval1);
row2=find(fval(:,2)==fval2);
xWayPoints=zeros(7,1);
xWayPoints(1,1)=0;
xWayPoints(2,1)=x(row2,1);
xWayPoints(3,1)=x(row2,2);
xWayPoints(4,1)=x(row2,3);
xWayPoints(5,1)=x(row2,4);
xWayPoints(6,1)=x(row2,5);
xWayPoints(7,1)=50;
yWayPoints(1,1)=12.5;
yWayPoints(2,1)=x(row2,6);
yWayPoints(3,1)=x(row2,7);
yWayPoints(4,1)=x(row2,8);
yWayPoints(5,1)=x(row2,9);
yWayPoints(6,1)=x(row2,10);
yWayPoints(7,1)=12.5;

 
%% Plot the optimal solution:
[Xgrid,Ygrid] = meshgrid(0:sizeX,0:sizeY);
hq = quiver(Xgrid,Ygrid,W_x,W_y,'k');
hold on;
xlabel('Units = 100 [km]');
axis equal tight
plot([0 sizeX],[sizeY sizeY]/2,'k.','markersize',16)
%% Add some color to make it more visible
L = (sqrt((Xgrid-sizeX).^2 + (Ygrid-sizeY/2).^2));
Favorability =((sizeX-Xgrid).*W_x+(sizeY/2-Ygrid).*W_y)./L; 
Favorability(~isfinite(Favorability)) = 0; 

hold on;
h_im = imagesc(Favorability); % This will be the background for the vector field

set(h_im,'Xdata',[0 sizeX],'Ydata',[0 sizeY]);
uistack(h_im,'bottom');

% Change the colormap...
colormap(interp1([0,1,2],[1 0 0; 1 1 1; 0 1 0],0:0.01:2));
caxis(max(abs(Favorability(:)))*[-1 1]);

h_colorbar = colorbar;
title(h_colorbar,'Tailwind (km/h)')
xlabel(h_colorbar,'Headwind')

h_wp = plot(xWayPoints,yWayPoints,'color','k','linestyle','none','marker','.','markersize',16);

PathPoints = WayPoints_To_Path([xWayPoints,yWayPoints],'cubic',sizeX,sizeY,501);
h_path = plot(PathPoints(:,1),PathPoints(:,2),'k','linewidth',2);
m=rectangle('Position',[27,0,16,16],'Curvature',[1,1],'LineWidth',2);
u=rectangle('Position',[32,5,6,6],'Curvature',[1,1],'LineWidth',2);
LineTime = getTimeFromPath(PathPoints,W_x,W_y,AirSpeed);
fprintf('Optimal Travel Time: %d hours, %.1f minutes\n',floor(LineTime),rem(LineTime,1)*60);



