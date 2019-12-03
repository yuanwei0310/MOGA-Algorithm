%**************************************************
%Test Problem 'TNK'
%**************************************************
clear
close all
n_pop=300; %population size
lb=[0 0]; %lower bound 
ub=[pi pi]; %upper bound
n_var=2; %number of variables
maxGen=500; %max generation
ShareNum=0.5;
a=1;
epsi=0.25;
X=zeros(n_pop,n_var);
R=10^6; %penalty number
for i=1:n_var
    X(:,i)=(ub(1,i)-lb(1,i))*rand(n_pop,1)+lb(1,i)*ones(n_pop,1); %randomly generate initial population
end
G=zeros(n_pop,2); %generate constraint matrix
Obj=zeros(n_pop,2);
ngen=0;
while(ngen<maxGen)
    ngen=ngen+1;
    G=zeros(n_pop,2); %generate constraint matrix
Obj=zeros(n_pop,2);
for j=1:n_pop
    if( X(j,2) == 0)
    G(j,1)=-X(j,1)^2-X(j,2)^2+1+0.1*cos(16*atan(Inf));
else
    G(j,1)=-X(j,1)^2-X(j,2)^2+1 + 0.1 * cos(16*atan(X(j,1)/X(j,2)));
    end
G(j,2)=(X(j,1)-0.5)^2 + (X(j,2)-0.5)^2 - 0.5;
    Obj(j,1)=X(j,1)+R*(max(0,G(j,1))+max(0,G(j,2)));
    Obj(j,2)=X(j,2)+R*(max(0,G(j,1))+max(0,G(j,2)));
end
    SharedFitness=CalSharedFitness(Obj,n_var,ShareNum,a,epsi);
    %for i=1:n_pop
    %SharedFitness(i,1)=SharedFitnessold(i,1)+R*(max(0,G(i,1))+max(0,G(i,2)));
    %end
    population=NewGA(SharedFitness,n_var,n_pop,lb,ub,X);
    X=population;
end
Obj
M=CalLayerRank(Obj)
P=find(M==1);
n=size(P,2);
for i=1:n
    pareto(i,:)=Obj(P(1,i),:);
end
x=pareto(:,1);
y=pareto(:,2);
plot(x,y,'*');

filename='TNK.txt';	
obj=get(gca,'children');
 x=get(obj,'xdata');
y=get(obj,'ydata');
 x=x(:);
 y=y(:);
M=[x y];
save(filename,'M','-ascii');

%calculate coverage difference
A=load('TNK.txt'); %load pareto point as a metric
B=sortrows(A,1);   %rearrange according to f1
% bad point 
pb_f1=1.5;
pb_f2=1.5;
cd1=0;
num1=size(B,1);
num=num1-1;
for i=1:num
    cd1=cd1+(B(i+1,1)-B(i,1))*B(i,2);
end
cd2=cd1+(pb_f1-B(num1,1))*B(num1,2)+B(1,1)*pb_f2;
cd_TNK=cd2/(pb_f1*pb_f2);
cd_TNK

%calculat accuracy of observed pareto
C=sortrows(A,-2);
ap1=0;
for i=1:num
    ap1=ap1+(C(i,2)-C(i+1,2))*(pb_f1-C(i,1));
end
ap2=ap1+pb_f1*(pb_f2-C(1,2))+(pb_f1-C(num1,1))*C(num1,2);
ap=((ap2+cd2)-(pb_f1*pb_f2))/(pb_f1*pb_f2);
ac_TNK=1/ap;
ac_TNK