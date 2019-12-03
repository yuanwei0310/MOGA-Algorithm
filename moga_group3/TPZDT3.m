clear all
close all
nvars = 2;
lb = zeros(1,30);
ub = ones(1,30);
options = gaoptimset;
options = gaoptimset(options,'PopulationSize', 50);
%options = gaoptimset(options,'Generations', 500);
options = gaoptimset(options,'Display', 'iter');
options = gaoptimset(options,'PlotFcns',@gaplotpareto);


options = gaoptimset(options, 'ParetoFraction', 0.8);
[~,fval] = gamultiobj(@TPZDT3_objfun,nvars,[],[],[],[],lb,ub,options);

filename='TPZDT3.txt';	
obj=get(gca,'children');
 x=get(obj,'xdata');
y=get(obj,'ydata');
 x=x(:);
 y=y(:);
M=[x y];
save(filename,'M','-ascii');

%calculate coverage difference
A=load('TPZDT3.txt'); %load pareto point as a metric
B=sortrows(A,1);   %rearrange according to f1
% bad point 
pb_f1=1;
pb_f2=1.2;
pg_f1=0;
pg_f2=-1;
cd1=0;
num1=size(B,1);
num=num1-1;
for i=1:num
    cd1=cd1+(B(i+1,1)-B(i,1))*(B(i,2)-pg_f2);
end
cd2=cd1+(pb_f1-B(num1,1))*(B(num1,2)-pg_f2)+B(1,1)*(pb_f2-pg_f2);
cd_TPZDT3=cd2/(pb_f1*(pb_f2-pg_f2));
cd_TPZDT3

%calculat accuracy of observed pareto
C=sortrows(A,-2);
ap1=0;
for i=1:num
    ap1=ap1+(C(i,2)-C(i+1,2))*(pb_f1-C(i,1));
end
ap2=ap1+pb_f1*(pb_f2-C(1,2))+(pb_f1-C(num1,1))*(C(num1,2)-pg_f2);
ap=((ap2+cd2)-(pb_f1*(pb_f2-pg_f2)))/(pb_f1*(pb_f2-pg_f2));
ac_TPZDT3=1/ap;
ac_TPZDT3

