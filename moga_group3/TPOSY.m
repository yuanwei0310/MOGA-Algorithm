clear all
close all
nvars = 6;
lb = [0 0 1 0 1 0];
ub = [10 10 5 6 5 10];
options = gaoptimset;
options = gaoptimset(options,'PopulationSize', 100);
%options = gaoptimset(options,'Generations', 500);
options = gaoptimset(options,'Display', 'iter');
options = gaoptimset(options,'PlotFcns',@gaplotpareto);


options = gaoptimset(options, 'ParetoFraction', 0.8);
[~,fval] = gamultiobj(@TPOSY_objfun,nvars,[],[],[],[],lb,ub,@TPOSYcons,options);

filename='TPOSY.txt';	
obj=get(gca,'children');
 x=get(obj,'xdata');
y=get(obj,'ydata');
 x=x(:);
 y=y(:);
M=[x y];
save(filename,'M','-ascii');

%calculate coverage difference
A=load('TPOSY.txt'); %load pareto point as a metric
B=sortrows(A,1);   %rearrange according to f1
% bad point 
pb_f1=0;
pb_f2=100;
pg_f1=-300;
cd1=0;
num1=size(B,1);
num=num1-1;
for i=1:num
    cd1=cd1+(B(i+1,1)-B(i,1))*B(i,2);
end
cd2=cd1+(pb_f1-B(num1,1))*B(num1,2)+(B(1,1)-pg_f1)*pb_f2;
cd_TPOSY=cd2/((pb_f1-pg_f1)*pb_f2);
cd_TPOSY

%calculat accuracy of observed pareto
C=sortrows(A,-2);
ap1=0;
for i=1:num
    ap1=ap1+(C(i,2)-C(i+1,2))*(pb_f1-C(i,1));
end
ap2=ap1+(pb_f1-pg_f1)*(pb_f2-C(1,2))+(pb_f1-C(num1,1))*C(num1,2);
ap=((ap2+cd2)-((pb_f1-pg_f1)*pb_f2))/((pb_f1-pg_f1)*pb_f2);
ac_TPOSY=1/ap;
ac_TPOSY