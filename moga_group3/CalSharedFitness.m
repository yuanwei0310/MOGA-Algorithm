% main nondominated sorting process, complete fitness assignment process
function [NicheCount,SharedFitness]=CalSharedFitness(Obj,n_var,ShareNum,a,epsi)
G=size(Obj,1); %popsize

%Run CalLayerRank.m to get the layer rank matrix
R=CalLayerRank(Obj);
R=R';
ObjDiff1=max(Obj(:,1))-min(Obj(:,1)); %Calculate f1max-f1min
ObjDiff2=max(Obj(:,2))-max(Obj(:,2)); %Calculate f2max-f2min
Fmin=n_var+epsi; %initialize fitness value
NicheCount=zeros(G,1); %initialize niche count
for n=1:max(R) % n is rank of each layer
    %if any(R==n)~=0
    m=sum(R==n); %calculate number of individuals in layer n
    lc=find(R==n); %find location of individuals in layer n
    for i=1:m
        for j=1:m
            %calculate distance dij
            d(lc(i),lc(j))=sqrt(((Obj(lc(i),1)-Obj(lc(j),1))/ObjDiff1)^2+((Obj(lc(i),2)-Obj(lc(j),2))/ObjDiff2)^2);
            if d(lc(i),lc(j))<=ShareNum
                sh(i,j)=1-(d(lc(i),lc(j))/ShareNum)^a;
            else
                sh(i,j)=0;
            end
            NicheCount(lc(i))=NicheCount(lc(i))+sh(i,j); %niche count is the sum of all j
        end
        Fit=Fmin-epsi; %fitness value of specific layer
        M(i)=Fit/NicheCount(lc(i)); %obtain the minimum fitness value for calculating higher layer's fitness value
        SharedFitness(lc(i),1)=Fit/NicheCount(lc(i));
    end
    %end
        Fmin=min(M); %fitness value of next layer
end

end
        