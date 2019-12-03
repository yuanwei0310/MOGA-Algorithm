%**************************************************
%Test Problem 'ZDT2'
%**************************************************
clear
close all
n_pop=50; %population size
lb=[0 0]; %lower bound 
ub=[1 1]; %upper bound
n_var=2; %number of variables
maxGen=500; %max generation
ShareNum=0.5;
a=1;
epsi=0.25;
X=zeros(n_pop,n_var);
for i=1:n_var
    X(:,i)=(ub(1,i)-lb(1,i))*rand(n_pop,1)+lb(1,i)*ones(n_pop,1);
end
Obj=zeros(n_pop,2);
ngen=0;
while(ngen<maxGen)
    ngen=ngen+1;
for j=1:n_pop
    Obj(j,1)=X(j,1);
    %g = 1 + 9*X(j,2)/29;
    Obj(j,2)=(1+9*X(j,2)/29)*(1-(X(j,1)/(1+9*X(j,2)/29))^2);
end
    SharedFitness=CalSharedFitness(Obj,n_var,ShareNum,a,epsi);
    population=NewGA(SharedFitness,n_var,n_pop,lb,ub,X);
    X=population;
end
Obj
R=CalLayerRank(Obj)
P=find(R==1);
n=size(P,2);
for i=1:n
    pareto(i,:)=Obj(P(1,i),:);
end
x=pareto(:,1);
y=pareto(:,2);
plot(x,y,'*')

filename='ZDT2.txt';	
obj=get(gca,'children');
 x=get(obj,'xdata');
y=get(obj,'ydata');
 x=x(:);
 y=y(:);
M=[x y];
save(filename,'M','-ascii');

%calculate coverage difference
A=load('ZDT2.txt'); %load pareto point as a metric
B=sortrows(A,1);   %rearrange according to f1
% bad point 
pb_f1=1.3;
pb_f2=1.3;
cd1=0;
num1=size(B,1);
num=num1-1;
for i=1:num
    cd1=cd1+(B(i+1,1)-B(i,1))*B(i,2);
end
cd2=cd1+(pb_f1-B(num1,1))*B(num1,2)+B(1,1)*pb_f2;
cd_ZDT2=cd2/(pb_f1*pb_f2);
cd_ZDT2

%calculat accuracy of observed pareto
C=sortrows(A,-2);
ap1=0;
for i=1:num
    ap1=ap1+(C(i,2)-C(i+1,2))*(pb_f1-C(i,1));
end
ap2=ap1+pb_f1*(pb_f2-C(1,2))+(pb_f1-C(num1,1))*C(num1,2);
ap=((ap2+cd2)-(pb_f1*pb_f2))/(pb_f1*pb_f2);
ac_ZDT2=1/ap;
ac_ZDT2

