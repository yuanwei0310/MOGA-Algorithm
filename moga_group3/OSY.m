%**************************************************
%Test Problem 'OSY'
%**************************************************
clear
close all
n_pop=200; %population size
lb=[0 0 1 0 1 0]; %lower bound 
ub=[10 10 5 6 5 10]; %upper bound
n_var=6; %number of variables
maxGen=500; %max generation
ShareNum=0.5;
a=1;
epsi=0.25;
X=zeros(n_pop,n_var);
R=1000000; %penalty number
for i=1:n_var
    X(:,i)=(ub(1,i)-lb(1,i))*rand(n_pop,1)+lb(1,i)*ones(n_pop,1); %randomly generate initial population
end
%G=zeros(n_pop,6); %generate constraint matrix
%Obj=zeros(n_pop,2);
ngen=0;
while(ngen<maxGen)
    ngen=ngen+1;
    G=zeros(n_pop,6); %generate constraint matrix
Obj=zeros(n_pop,2);
for j=1:n_pop
    G(j,1)=1-(X(j,1)+X(j,2))/2;
    G(j,2)=(X(j,1)+X(j,2))/6-1;
    G(j,3)=(X(j,2)-X(j,1))/2-1;
    G(j,4)=(X(j,1)-3*X(j,2))/2-1;
    G(j,5)=((X(j,3)-3)^2+X(j,4))/4-1;
    G(j,6)=1-((X(j,5)-3)^2+X(j,6))/4;
    Obj(j,1)=-(25*(X(j,1)-2)^2+(X(j,2)-2)^2+(X(j,3)-1)^2+(X(j,4)-4)^2+(X(j,5)-1)^2)+R*(max(0,G(j,1))+max(0,G(j,2))+max(0,G(j,3))+max(0,G(j,4))+max(0,G(j,5))+max(0,G(j,6)));
    Obj(j,2)=X(j,1)^2+X(j,2)^2+X(j,3)^2+X(j,4)^2+X(j,5)^2+X(j,6)^2+R*(max(0,G(j,1))+max(0,G(j,2))+max(0,G(j,3))+max(0,G(j,4))+max(0,G(j,5))+max(0,G(j,6)));
end
    SharedFitness=CalSharedFitness(Obj,n_var,ShareNum,a,epsi);
    %for i=1:n_pop
   %SharedFitness(i,1)=SharedFitnessold(i,1);%+R*(max(0,G(i,1))+max(0,G(i,2))+max(0,G(i,3))+max(0,G(i,4))+max(0,G(i,5))+max(0,G(i,6)));
    %end
    population=NewGA(SharedFitness,n_var,n_pop,lb,ub,X);
    X=population;
end
D=CalLayerRank(Obj);
P=find(D==1);
n=size(P,2);
for i=1:n
    pareto(i,:)=Obj(P(1,i),:);
end
x=pareto(:,1);
y=pareto(:,2);
plot(x,y,'o');
filename='OSY.txt';	
obj=get(gca,'children');
 x=get(obj,'xdata');
y=get(obj,'ydata');
 x=x(:);
 y=y(:);
M=[x y];
save(filename,'M','-ascii');

%calculate coverage difference
A=load('OSY.txt'); %load pareto point as a metric
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
cd_OSY=cd2/((pb_f1-pg_f1)*pb_f2);
cd_OSY

%calculat accuracy of observed pareto
C=sortrows(A,-2);
ap1=0;
for i=1:num
    ap1=ap1+(C(i,2)-C(i+1,2))*(pb_f1-C(i,1));
end
ap2=ap1+(pb_f1-pg_f1)*(pb_f2-C(1,2))+(pb_f1-C(num1,1))*C(num1,2);
ap=((ap2+cd2)-((pb_f1-pg_f1)*pb_f2))/((pb_f1-pg_f1)*pb_f2);
ac_OSY=1/ap;
ac_OSY
