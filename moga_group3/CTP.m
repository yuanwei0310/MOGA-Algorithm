%**************************************************
%Test Problem 'CTP'
%**************************************************
clear
close all
n_pop=300; %population size
lb=[0 -5 -5 -5 -5 -5 -5 -5 -5 -5]; %lower bound 
ub=[1 5 5 5 5 5 5 5 5 5]; %upper bound
n_var=10; %number of variables
maxGen=500; %max generation
ShareNum=0.5;
a=1;
epsi=0.25;
X=zeros(n_pop,n_var);
R=10^5; %penalty number
for i=1:n_var
    X(:,i)=(ub(1,i)-lb(1,i))*rand(n_pop,1)+lb(1,i)*ones(n_pop,1); %randomly generate initial population
end
G=zeros(n_pop,1); %generate constraint matrix
Obj=zeros(n_pop,2);
ngen=0;
while(ngen<maxGen)
    ngen=ngen+1;
    G=zeros(n_pop,1); %generate constraint matrix
Obj=zeros(n_pop,2);
Y=zeros(n_pop,2);
for j=1:n_pop
    th=-0.2*pi; a=0.2; b=10; d=6; e=1;
    g=abs(1+(sum(X(j,2:10))).^(0.25));
    Y(j,1)=X(j,1);
    Y(j,2)=g*(1-sqrt(X(j,1)/g));
    G(j,1)=-cos(th)*(Y(j,2)-e)+sin(th)*Y(j,1)+a*(abs(sin(b*pi*(sin(th)*(Y(j,2)-e)+cos(th)*Y(j,1)).^1))).^6;
    Obj(j,1)=Y(j,1)+R*max(0,G(j,1));
    Obj(j,2)=Y(j,2)+R*max(0,G(j,1));
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

filename='CTP.txt';	
obj=get(gca,'children');
 x=get(obj,'xdata');
y=get(obj,'ydata');
 x=x(:);
 y=y(:);
M=[x y];
save(filename,'M','-ascii');

%calculate coverage difference
A=load('CTP.txt'); %load pareto point as a metric
B=sortrows(A,1);   %rearrange according to f1
% bad point 
pb_f1=1;
pb_f2=3;
cd1=0;
num1=size(B,1);
num=num1-1;
for i=1:num
    cd1=cd1+(B(i+1,1)-B(i,1))*B(i,2);
end
cd2=cd1+(pb_f1-B(num1,1))*B(num1,2)+B(1,1)*pb_f2;
cd_CTP=cd2/(pb_f1*pb_f2);
cd_CTP

%calculat accuracy of observed pareto
C=sortrows(A,-2);
ap1=0;
for i=1:num
    ap1=ap1+(C(i,2)-C(i+1,2))*(pb_f1-C(i,1));
end
ap2=ap1+pb_f1*(pb_f2-C(1,2))+(pb_f1-C(num1,1))*C(num1,2);
ap=((ap2+cd2)-(pb_f1*pb_f2))/(pb_f1*pb_f2);
ac_CTP=1/ap;
ac_CTP